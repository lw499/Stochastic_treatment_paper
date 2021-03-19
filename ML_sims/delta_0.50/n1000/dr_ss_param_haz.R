library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(detectCores()-1)                                              
registerDoParallel(cl)                                                                          
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  library(geepack);library(MASS);library(ResourceSelection);library(SuperLearner)
  library(dplyr); library(glm2);
  library(data.table)
  #library(reshape2)  #do not use for data frame only
  setDTthreads(1)
  
  logit <- function(term) {
    return( ifelse(!is.na(term),log(term/(1-term)),NA) )
  }
  
  EXPIT <- function(term) {
    return( ifelse(!is.na(term),exp(term)/(1+exp(term)),NA) )
  }
  
  source("datagen.R")
  set.seed(1129)
  seeds = floor(runif(1000)*10^8);
  set.seed(seeds[m])
  
  n <- 1000
  K <- 5
  delta=0.5
  beta0=2; beta1=1; beta2=-1; beta3=0.5; beta4=1;
  gamma0=1.5; gamma1=-1; gamma2=-0.5; gamma3=1; gamma4=0.25; gamma5=1;
  alpha0=-3; alpha1=1; alpha2=-0.5; alpha3=0.25; alpha4=0.5; alpha5=0.25;
  theta0=-1; theta1=2; theta2=-2; theta3=0.5; theta4=-0.25; theta5=0.5; theta6=0.75;
  eta0=-4; eta1=-1; eta2=-1; eta3=-0.5; eta4=0;
  sigma=1
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen(ind, K=K,
            alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3, alpha4=alpha4, alpha5=alpha5, 
            beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4,
            theta0=theta0, theta1=theta1, theta2=theta2, theta3=theta3, theta4=theta4, theta5=theta5, theta6=theta6, 
            eta0=eta0, eta1=eta1, eta2=eta2, eta3=eta3, eta4=eta4, 
            gamma0=gamma0, gamma1=gamma1, gamma2=gamma2, gamma3=gamma3, gamma4=gamma4, gamma5=gamma5, sigma=sigma)
  })
  
  dffull <- rbindlist(df)
  
  dffull[, paste("lag_A") := shift(A, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_C") := shift(C, 1, NA, type='lag'), by=id]
  dffull$lag_A = ifelse(dffull$t0==0, 0 , dffull$lag_A)
  dffull$lag_C = ifelse(dffull$t0==0, 0 , dffull$lag_C)
  
  #mtry_seq    <- seq(2, 4, 1)
  #n.trees_seq <- 500#seq(500, 1000, 500) 
  learners <- create.Learner("SL.ranger", tune = list(num.trees = 3000, min.mode.size=50))  
  #deg.gam = c(4,5,6)
  #learners <- create.Learner("SL.gam", tune = list(deg.gam = c(1,3,4)))
  sl.lib = c("SL.gam","SL.glm", "SL.bayesglm", "SL.glm.interaction","SL.ranger")
  
  #sl.lib = c("SL.gam","SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.glm.interaction","SL.step.forward","SL.ranger")
  
  set.seed(100388)
  train_ind = sample(unique(dffull$id),size = n/2)
  
  #fit afit and cfit:
  tmpdffull = dffull[!is.na(dffull$A) & dffull$lag_A==0,]
  tmpdffull  = tmpdffull[tmpdffull$id %in% train_ind,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  y = tmpdffull$A
  set.seed(1003)
  afit = glm2(A ~ L1 + L2 + L1*L2 + W1 + W2 + I(abs(W2)), family=binomial(), data = tmpdffull )

  tmpdffull = dffull[!is.na(dffull$C) & dffull$lag_C==0,]
  tmpdffull  = tmpdffull[tmpdffull$id %in% train_ind,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "A" = tmpdffull$A, "t0" = tmpdffull$t0))
  y = tmpdffull$C
  set.seed(1003)
  cfit = glm2(C ~ A + L1 + I(sqrt(abs(W2*L2))) + I(abs(L2)/(1+exp(W2))), family=binomial(), data = tmpdffull )
  
  #fit afit and cfit (swap):
  tmpdffull = dffull[!is.na(dffull$A) & dffull$lag_A==0,]
  tmpdffull  = tmpdffull[!tmpdffull$id %in% train_ind,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  y = tmpdffull$A
  set.seed(1003)
  afitb = glm2(A ~ L1 + L2 + L1*L2 + W1 + W2 + I(abs(W2)), family=binomial(), data = tmpdffull )
  
  tmpdffull = dffull[!is.na(dffull$C) & dffull$lag_C==0,]
  tmpdffull  = tmpdffull[!tmpdffull$id %in% train_ind,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "A" = tmpdffull$A, "t0" = tmpdffull$t0))
  y = tmpdffull$C
  set.seed(1003)
  cfitb = glm2(C ~ A + L1 + I(sqrt(abs(W2*L2))) + I(abs(L2)/(1+exp(W2))), family=binomial(), data = tmpdffull )
  
  #fit hazard
  tmpdffull = dffull[!is.na(dffull$Y) & dffull$C==0,]
  tmpdffull  = tmpdffull[tmpdffull$id %in% train_ind,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  A = tmpdffull$A; y = tmpdffull$Y
  set.seed(1003)
  yfitog = glm2(Y ~ A + L1 + I(L1*L2) + W1 + I(abs(L2+W2)^1.5), family=binomial(), data = tmpdffull)
  #fit hazard (swap train and test)
  tmpdffull = dffull[!is.na(dffull$Y) & dffull$C==0,]
  tmpdffull  = tmpdffull[!tmpdffull$id %in% train_ind,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  A = tmpdffull$A; y = tmpdffull$Y
  set.seed(1003)
  yfitogb = glm2(Y ~ A + L1 + I(L1*L2) + W1 + I(abs(L2+W2)^1.5), family=binomial(), data = tmpdffull)

  #iteration 1 (this is only for computationl simplicity. 
  #Even though fit for everyone, we will eventually only use those in training):
  dffulltmp = dffull[,c("L1","L2","W1","W2","t0")]
  dffull$pred_obs = predict(afit, newdata = dffulltmp, type="response")
  dffull$pred_obs = ifelse(dffull$A==1, dffull$pred_obs, 1-dffull$pred_obs)
  dffull$pred_obs = ifelse(!is.na(dffull$lag_A) & dffull$lag_A==1, 1, dffull$pred_obs)
  dffulltmp = dffull[,c("L1","L2","W1","W2","A","t0")]
  dffull$pred_obsc = predict(cfit, newdata = dffulltmp, type="response")
  dffull$pred_obsc = ifelse(dffull$C==1, NA, 1-dffull$pred_obsc)
  dffull$pred_obsc = ifelse(!is.na(dffull$lag_C) & dffull$lag_C==1, 1, dffull$pred_obsc)
  
  dffull$fint = (dffull$L1*delta + 1 - dffull$L1)*dffull$pred_obs + dffull$L1*dffull$A*(1-delta)
  
  #dffull$lag_A =  dffull$lag_C = NULL
  dffullwide = dcast(dffull, id + W1 + W2 ~ t0, value.var = c("L1","L2","A","C","Y","pred_obs","pred_obsc","fint"))
  
  tmpdata = dffullwide
  tmpdata$Y_5 = tmpdata$Y_4
  tmpdata$Y_4 = tmpdata$Y_3
  tmpdata$Y_3 = tmpdata$Y_2
  tmpdata$Y_2 = tmpdata$Y_1
  tmpdata$Y_1 = tmpdata$Y_0
  tmpdata$Y_0 = NULL

  tmpdata$pi0 = tmpdata$pred_obs_0
  tmpdata$pi1 = tmpdata$pi0*tmpdata$pred_obs_1
  tmpdata$pi2 = tmpdata$pi1*tmpdata$pred_obs_2
  tmpdata$pi3 = tmpdata$pi2*tmpdata$pred_obs_3
  tmpdata$pi4 = tmpdata$pi3*tmpdata$pred_obs_4
  
  tmpdata$fintc0 = tmpdata$fint_0
  tmpdata$fintc1 = tmpdata$fintc0*tmpdata$fint_1
  tmpdata$fintc2 = tmpdata$fintc1*tmpdata$fint_2
  tmpdata$fintc3 = tmpdata$fintc2*tmpdata$fint_3
  tmpdata$fintc4 = tmpdata$fintc3*tmpdata$fint_4
  
  tmpdata$pi0c = tmpdata$pred_obsc_0
  tmpdata$pi1c = tmpdata$pi0c*tmpdata$pred_obsc_1
  tmpdata$pi2c = tmpdata$pi1c*tmpdata$pred_obsc_2
  tmpdata$pi3c = tmpdata$pi2c*tmpdata$pred_obsc_3
  tmpdata$pi4c = tmpdata$pi3c*tmpdata$pred_obsc_4
  
  ##################
  ######time 5######
  ##################
  tmp_train =  tmpdata[tmpdata$id %in% train_ind,]
  #ydat = tmp_train[!is.na(tmp_train$C_4) & tmp_train$C_4==0,]; 
  #l = as.data.frame(cbind("L1_4" = ydat$L1_4, "L2_4" = ydat$L2_4, "W1" = ydat$W1, "W2" = ydat$W2))
  #A_4 = ydat$A_4; y = ydat$Y_5
  #set.seed(1003)
  #yfitog = SuperLearner(Y=y, X=cbind(l,A_4), SL.library=sl.lib, family=binomial)
  y4dat = tmp_train[!is.na(tmp_train$C_3) & tmp_train$C_3==0,]; 
  ydat = y4dat; ydat = as.data.table(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4,  "W1" = ydat$W1, "W2" = ydat$W2, "t0" = 4))
  ydat[is.na(ydat)] = 0;  ydat$A = 1;
  ypred = predict(yfitog, newdata = ydat, type="response"); 
  ydat2 = y4dat; ydat2 = as.data.table(cbind("L1" = ydat2$L1_4, "L2" = ydat2$L2_4, "W1" = ydat2$W1, "W2" = ydat2$W2, "t0" = 4, "A" = ydat2$A_4))
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict(yfitog, newdata = ydat2, type="response")*(y4dat$L1_4*delta + 1-y4dat$L1_4) + ypred*(1-delta)*(y4dat$L1_4)
  y4dat$y4pred = ifelse(y4dat$Y_4==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_3" = y4dat$L1_3, "L2_3" = y4dat$L2_3, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_3 = y4dat$A_3; y = y4dat$y4pred
  set.seed(1003); yfitextra1 = SuperLearner(Y=y, X=cbind(l,A_3), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmp_train[!is.na(tmp_train$C_2) & tmp_train$C_2==0,]; 
  ydat = y4dat; ydat = ydat[,c("L1_3","L2_3","W1","W2","A_3")]
  ydat[is.na(ydat)] = 0; ydat$A_3 = 1;
  ypred = predict.SuperLearner(yfitextra1, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_3","L2_3","W1","W2","A_3")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra1, newdata = ydat2, type="response", onlySL = T)$pred*(y4dat$L1_3*delta + 1-y4dat$L1_3) + ypred*(1-delta)*(y4dat$L1_3)
  y4dat$y4pred = ifelse(y4dat$Y_3==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_2" = y4dat$L1_2, "L2_2" = y4dat$L2_2, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_2 = y4dat$A_2; y = y4dat$y4pred
  set.seed(1003); yfitextra2 = SuperLearner(Y=y, X=cbind(l,A_2), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmp_train[!is.na(tmp_train$C_1) & tmp_train$C_1==0,]; 
  ydat = y4dat; ydat = ydat[,c("L1_2","L2_2","W1","W2","A_2")]
  ydat[is.na(ydat)] = 0; ydat$A_2 = 1;
  ypred = predict.SuperLearner(yfitextra2, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_2","L2_2","W1","W2","A_2")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra2, newdata = ydat2, type="response", onlySL = T)$pred*(y4dat$L1_2*delta + 1-y4dat$L1_2) + ypred*(1-delta)*(y4dat$L1_2)
  y4dat$y4pred = ifelse(y4dat$Y_2==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_1" = y4dat$L1_1, "L2_1" = y4dat$L2_1, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_1 = y4dat$A_1; y = y4dat$y4pred
  set.seed(1003); yfitextra3 = SuperLearner(Y=y, X=cbind(l,A_1), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmp_train[!is.na(tmp_train$C_0) & tmp_train$C_0==0,]; 
  ydat = y4dat; ydat = ydat[,c("L1_1","L2_1","W1","W2","A_1")]
  ydat[is.na(ydat)] = 0; ydat$A_1 = 1;
  ypred = predict.SuperLearner(yfitextra3, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_1","L2_1","W1","W2","A_1")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra3, newdata = ydat2, type="response", onlySL = T)$pred*(y4dat$L1_1*delta + 1-y4dat$L1_1) + ypred*(1-delta)*(y4dat$L1_1)
  y4dat$y4pred = ifelse(y4dat$Y_1==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_0" = y4dat$L1_0, "L2_0" = y4dat$L2_0, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_0 = y4dat$A_0; y = y4dat$y4pred
  set.seed(1003); yfitextra4 = SuperLearner(Y=y, X=cbind(l,A_0), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmp_train; 
  ydat = y4dat; ydat = ydat[,c("L1_0","L2_0","W1","W2","A_0")]
  ydat[is.na(ydat)] = 0; ydat$A_0 = 1;
  ypred = predict.SuperLearner(yfitextra4, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_0","L2_0","W1","W2","A_0")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra4, newdata = ydat2, type="response", onlySL = T)$pred*(y4dat$L1_0*delta + 1-y4dat$L1_0) + ypred*(1-delta)*(y4dat$L1_0)
  
  ### Predict and refit on test dataset
  ##################
  ######time 5######
  ##################
  tmp_test =  tmpdata[!tmpdata$id %in% train_ind,]
  ydat = tmp_test[!is.na(tmp_test$C_4) & tmp_test$C_4==0,]; 
  #l = as.data.frame(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4, "W1" = ydat$W1, "W2" = ydat$W2, "t0" = 4))
  #A = ydat$A_4; y = ydat$ypred
  tmp = as.data.frame(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4, "W1" = ydat$W1, "W2" = ydat$W2, "t0" = 4,  "A" = ydat$A_4)); tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict(yfitog, newdata = tmp, type="response"))
  yfit2 = glm2(Y_5 ~ 1, family = binomial(), offset = ypred0, weights = (fintc4)/(pi4*pi4c), data = ydat) 
  ydat = tmp_test[!is.na(tmp_test$C_3) & tmp_test$C_3==0,]; 
  ndat = ydat; ndat = as.data.frame(cbind("L1" = ndat$L1_4, "L2" = ndat$L2_4, "W1" = ndat$W1, "W2" = ndat$W2, "t0" = 4, "A" = ndat$A_4));
  ndat[is.na(ndat)] = 0;  ndat$A = 1;
  ndat$ypred0 = qlogis(predict(yfitog, newdata = ndat, type="response")); 
  ypred = predict(yfit2, newdata = ndat, type="response"); 
  ndat2 = ydat; ndat2 = as.data.frame(cbind("L1" = ndat2$L1_4, "L2" = ndat2$L2_4, "W1" = ndat2$W1, "W2" = ndat2$W2, "t0" = 4, "A" = ndat2$A_4)); ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict(yfitog, newdata = ndat2, type="response")); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_4*delta + 1-ydat$L1_4) + ypred*(1-delta)*(ydat$L1_4)
  ydat$ypred = ifelse(ydat$Y_4==0,0,ydat$ypred); 
  
  l = as.data.frame(cbind("L1_3" = ydat$L1_3, "L2_3" = ydat$L2_3, "W1" = ydat$W1, "W2" = ydat$W2))
  A_3 = ydat$A_3; y = ydat$ypred
  tmp  = ydat[,c("L1_3","L2_3","W1","W2","A_3")]; tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitextra1, newdata = tmp, onlySL = T)$pred)
  yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = (fintc3)/(pi3*pi3c), data = ydat) 
  ydat = tmp_test[!is.na(tmp_test$C_2) & tmp_test$C_2==0,]; 
  ndat = ydat; ndat = ndat[,c("L1_3","L2_3","W1","W2","A_3")]
  ndat[is.na(ndat)] = 0; ndat$A_3 = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitextra1, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response");   
  ndat2 = ydat; ndat2 = ndat2[,c("L1_3","L2_3","W1","W2","A_3")]; ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitextra1, newdata = ndat2, type="response", onlySL = T)$pred); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_3*delta + 1-ydat$L1_3) + ypred*(1-delta)*(ydat$L1_3)  
  ydat$ypred = ifelse(ydat$Y_3==0,0,ydat$ypred); 
  
  l = as.data.frame(cbind("L1_2" = ydat$L1_2, "L2_2" = ydat$L2_2, "W1" = ydat$W1, "W2" = ydat$W2))
  A_2 = ydat$A_2; y = ydat$ypred
  tmp = ydat[,c("L1_2","L2_2","W1","W2","A_2")]; tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitextra2, newdata = tmp, onlySL = T)$pred) 
  yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = (fintc2)/(pi2*pi2c), data = ydat) 
  ydat = tmp_test[!is.na(tmp_test$C_1) & tmp_test$C_1==0,]; 
  ndat = ydat; ndat = ndat[,c("L1_2","L2_2","W1","W2","A_2")]
  ndat[is.na(ndat)] = 0; ndat$A_2 = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitextra2, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response");   
  ndat2 = ydat; ndat2 = ndat2[,c("L1_2","L2_2","W1","W2","A_2")]; ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitextra2, newdata = ndat2, type="response", onlySL = T)$pred); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_2*delta + 1-ydat$L1_2) + ypred*(1-delta)*(ydat$L1_2)
  ydat$ypred = ifelse(ydat$Y_2==0,0,ydat$ypred); 
  
  l = as.data.frame(cbind("L1_1" = ydat$L1_1, "L2_1" = ydat$L2_1, "W1" = ydat$W1, "W2" = ydat$W2))
  A_1 = ydat$A_1; y = ydat$ypred
  tmp = ydat[,c("L1_1","L2_1","W1","W2","A_1")]; tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitextra3, newdata = tmp, onlySL = T)$pred) 
  yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = (fintc1)/(pi1*pi1c), data = ydat)
  ydat = tmp_test[!is.na(tmp_test$C_0) & tmp_test$C_0==0,]; 
  ndat = ydat; ndat = ndat[,c("L1_1","L2_1","W1","W2","A_1")]
  ndat[is.na(ndat)] = 0; ndat$A_1 = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitextra3, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response");   
  ndat2 = ydat; ndat2 = ndat2[,c("L1_1","L2_1","W1","W2","A_1")]; ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitextra3, newdata = ndat2, type="response", onlySL = T)$pred); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_1*delta + 1-ydat$L1_1) + ypred*(1-delta)*(ydat$L1_1)
  ydat$ypred = ifelse(ydat$Y_1==0,0,ydat$ypred); 
  
  l = as.data.frame(cbind("L1_0" = ydat$L1_0, "L2_0" = ydat$L2_0, "W1" = ydat$W1, "W2" = ydat$W2))
  A_0 = ydat$A_0; y = ydat$ypred
  tmp = ydat[,c("L1_0","L2_0","W1","W2","A_0")]; tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitextra4, newdata = tmp, onlySL = T)$pred) 
  yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = (fintc0)/(pi0*pi0c), data = ydat)
  ydat = tmp_test; 
  ndat = ydat; ndat = ndat[,c("L1_0","L2_0","W1","W2","A_0")]
  ndat[is.na(ndat)] = 0; ndat$A_0 = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitextra4, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response");   
  ndat2 = ydat; ndat2 = ndat2[,c("L1_0","L2_0","W1","W2","A_0")]; ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitextra4, newdata = ndat2, type="response", onlySL = T)$pred); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_0*delta + 1-ydat$L1_0) + ypred*(1-delta)*(ydat$L1_0)
  
  meany12_1 = mean(ydat$ypred)
  
  ## REPEAT (swap the train and test dataset)
  dffulltmp = dffull[,c("L1","L2","W1","W2","t0")]
  dffull$pred_obs = predict(afitb, newdata = dffulltmp, type="response")
  dffull$pred_obs = ifelse(dffull$A==1, dffull$pred_obs, 1-dffull$pred_obs)
  dffull$pred_obs = ifelse(!is.na(dffull$lag_A) & dffull$lag_A==1, 1, dffull$pred_obs)
  dffulltmp = dffull[,c("L1","L2","W1","W2","A","t0")]
  dffull$pred_obsc = predict(cfitb, newdata = dffulltmp, type="response")
  dffull$pred_obsc = ifelse(dffull$C==1, NA, 1-dffull$pred_obsc)
  dffull$pred_obsc = ifelse(!is.na(dffull$lag_C) & dffull$lag_C==1, 1, dffull$pred_obsc)
  
  dffull$fint = (dffull$L1*delta + 1 - dffull$L1)*dffull$pred_obs + dffull$L1*dffull$A*(1-delta)
  
  dffullwide = dcast(dffull, id + W1 + W2 ~ t0, value.var = c("L1","L2","A","C","Y","pred_obs","pred_obsc","fint"))
  
  tmpdata = dffullwide
  tmpdata$Y_5 = tmpdata$Y_4
  tmpdata$Y_4 = tmpdata$Y_3
  tmpdata$Y_3 = tmpdata$Y_2
  tmpdata$Y_2 = tmpdata$Y_1
  tmpdata$Y_1 = tmpdata$Y_0
  tmpdata$Y_0 = NULL
  
  tmpdata$pi0 = tmpdata$pred_obs_0
  tmpdata$pi1 = tmpdata$pi0*tmpdata$pred_obs_1
  tmpdata$pi2 = tmpdata$pi1*tmpdata$pred_obs_2
  tmpdata$pi3 = tmpdata$pi2*tmpdata$pred_obs_3
  tmpdata$pi4 = tmpdata$pi3*tmpdata$pred_obs_4
  
  tmpdata$fintc0 = tmpdata$fint_0
  tmpdata$fintc1 = tmpdata$fintc0*tmpdata$fint_1
  tmpdata$fintc2 = tmpdata$fintc1*tmpdata$fint_2
  tmpdata$fintc3 = tmpdata$fintc2*tmpdata$fint_3
  tmpdata$fintc4 = tmpdata$fintc3*tmpdata$fint_4
  
  tmpdata$pi0c = tmpdata$pred_obsc_0
  tmpdata$pi1c = tmpdata$pi0c*tmpdata$pred_obsc_1
  tmpdata$pi2c = tmpdata$pi1c*tmpdata$pred_obsc_2
  tmpdata$pi3c = tmpdata$pi2c*tmpdata$pred_obsc_3
  tmpdata$pi4c = tmpdata$pi3c*tmpdata$pred_obsc_4
  
  ##################
  ######time 5######
  ##################
  tmp_train =  tmpdata[!tmpdata$id %in% train_ind,]
  #ydat = tmp_train[!is.na(tmp_train$C_4) & tmp_train$C_4==0,]; 
  #l = as.data.frame(cbind("L1_4" = ydat$L1_4, "L2_4" = ydat$L2_4, "W1" = ydat$W1, "W2" = ydat$W2))
  #A_4 = ydat$A_4; y = ydat$Y_5
  #set.seed(1003)
  #yfitog = SuperLearner(Y=y, X=cbind(l,A_4), SL.library=sl.lib, family=binomial)
  y4dat = tmp_train[!is.na(tmp_train$C_3) & tmp_train$C_3==0,]; 
  ydat = y4dat; ydat = as.data.table(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4,  "W1" = ydat$W1, "W2" = ydat$W2, "t0" = 4))
  ydat[is.na(ydat)] = 0;  ydat$A = 1;
  ypred = predict(yfitogb, newdata = ydat, type="response"); 
  ydat2 = y4dat; ydat2 = as.data.table(cbind("L1" = ydat2$L1_4, "L2" = ydat2$L2_4, "W1" = ydat2$W1, "W2" = ydat2$W2, "t0" = 4, "A" = ydat2$A_4))
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict(yfitogb, newdata = ydat2, type="response")*(y4dat$L1_4*delta + 1-y4dat$L1_4) + ypred*(1-delta)*(y4dat$L1_4)
  y4dat$y4pred = ifelse(y4dat$Y_4==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_3" = y4dat$L1_3, "L2_3" = y4dat$L2_3, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_3 = y4dat$A_3; y = y4dat$y4pred
  set.seed(1003); yfitextra1 = SuperLearner(Y=y, X=cbind(l,A_3), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmp_train[!is.na(tmp_train$C_2) & tmp_train$C_2==0,]; 
  ydat = y4dat; ydat = ydat[,c("L1_3","L2_3","W1","W2","A_3")]
  ydat[is.na(ydat)] = 0; ydat$A_3 = 1;
  ypred = predict.SuperLearner(yfitextra1, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_3","L2_3","W1","W2","A_3")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra1, newdata = ydat2, type="response", onlySL = T)$pred*(y4dat$L1_3*delta + 1-y4dat$L1_3) + ypred*(1-delta)*(y4dat$L1_3)
  y4dat$y4pred = ifelse(y4dat$Y_3==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_2" = y4dat$L1_2, "L2_2" = y4dat$L2_2, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_2 = y4dat$A_2; y = y4dat$y4pred
  set.seed(1003); yfitextra2 = SuperLearner(Y=y, X=cbind(l,A_2), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmp_train[!is.na(tmp_train$C_1) & tmp_train$C_1==0,]; 
  ydat = y4dat; ydat = ydat[,c("L1_2","L2_2","W1","W2","A_2")]
  ydat[is.na(ydat)] = 0; ydat$A_2 = 1;
  ypred = predict.SuperLearner(yfitextra2, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_2","L2_2","W1","W2","A_2")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra2, newdata = ydat2, type="response", onlySL = T)$pred*(y4dat$L1_2*delta + 1-y4dat$L1_2) + ypred*(1-delta)*(y4dat$L1_2)
  y4dat$y4pred = ifelse(y4dat$Y_2==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_1" = y4dat$L1_1, "L2_1" = y4dat$L2_1, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_1 = y4dat$A_1; y = y4dat$y4pred
  set.seed(1003); yfitextra3 = SuperLearner(Y=y, X=cbind(l,A_1), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmp_train[!is.na(tmp_train$C_0) & tmp_train$C_0==0,]; 
  ydat = y4dat; ydat = ydat[,c("L1_1","L2_1","W1","W2","A_1")]
  ydat[is.na(ydat)] = 0; ydat$A_1 = 1;
  ypred = predict.SuperLearner(yfitextra3, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_1","L2_1","W1","W2","A_1")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra3, newdata = ydat2, type="response", onlySL = T)$pred*(y4dat$L1_1*delta + 1-y4dat$L1_1) + ypred*(1-delta)*(y4dat$L1_1)
  y4dat$y4pred = ifelse(y4dat$Y_1==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_0" = y4dat$L1_0, "L2_0" = y4dat$L2_0, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_0 = y4dat$A_0; y = y4dat$y4pred
  set.seed(1003); yfitextra4 = SuperLearner(Y=y, X=cbind(l,A_0), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmp_train; 
  ydat = y4dat; ydat = ydat[,c("L1_0","L2_0","W1","W2","A_0")]
  ydat[is.na(ydat)] = 0; ydat$A_0 = 1;
  ypred = predict.SuperLearner(yfitextra4, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_0","L2_0","W1","W2","A_0")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra4, newdata = ydat2, type="response", onlySL = T)$pred*(y4dat$L1_0*delta + 1-y4dat$L1_0) + ypred*(1-delta)*(y4dat$L1_0)
  
  ### Predict and refit on test dataset
  ##################
  ######time 5######
  ##################
  tmp_test =  tmpdata[tmpdata$id %in% train_ind,]
  ydat = tmp_test[!is.na(tmp_test$C_4) & tmp_test$C_4==0,]; 
  #l = as.data.frame(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4, "W1" = ydat$W1, "W2" = ydat$W2, "t0" = 4))
  #A = ydat$A_4; y = ydat$ypred
  tmp = as.data.frame(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4, "W1" = ydat$W1, "W2" = ydat$W2, "t0" = 4,  "A" = ydat$A_4)); tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict(yfitogb, newdata = tmp, type="response"))
  yfit2 = glm2(Y_5 ~ 1, family = binomial(), offset = ypred0, weights = (fintc4)/(pi4*pi4c), data = ydat) 
  ydat = tmp_test[!is.na(tmp_test$C_3) & tmp_test$C_3==0,]; 
  ndat = ydat; ndat = as.data.frame(cbind("L1" = ndat$L1_4, "L2" = ndat$L2_4, "W1" = ndat$W1, "W2" = ndat$W2, "t0" = 4, "A" = ndat$A_4));
  ndat[is.na(ndat)] = 0;  ndat$A = 1;
  ndat$ypred0 = qlogis(predict(yfitogb, newdata = ndat, type="response")); 
  ypred = predict(yfit2, newdata = ndat, type="response"); 
  ndat2 = ydat; ndat2 = as.data.frame(cbind("L1" = ndat2$L1_4, "L2" = ndat2$L2_4, "W1" = ndat2$W1, "W2" = ndat2$W2, "t0" = 4, "A" = ndat2$A_4)); ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict(yfitogb, newdata = ndat2, type="response")); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_4*delta + 1-ydat$L1_4) + ypred*(1-delta)*(ydat$L1_4)
  ydat$ypred = ifelse(ydat$Y_4==0,0,ydat$ypred); 
  
  l = as.data.frame(cbind("L1_3" = ydat$L1_3, "L2_3" = ydat$L2_3, "W1" = ydat$W1, "W2" = ydat$W2))
  A_3 = ydat$A_3; y = ydat$ypred
  tmp  = ydat[,c("L1_3","L2_3","W1","W2","A_3")]; tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitextra1, newdata = tmp, onlySL = T)$pred)
  yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = (fintc3)/(pi3*pi3c), data = ydat) 
  ydat = tmp_test[!is.na(tmp_test$C_2) & tmp_test$C_2==0,]; 
  ndat = ydat; ndat = ndat[,c("L1_3","L2_3","W1","W2","A_3")]
  ndat[is.na(ndat)] = 0; ndat$A_3 = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitextra1, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response");   
  ndat2 = ydat; ndat2 = ndat2[,c("L1_3","L2_3","W1","W2","A_3")]; ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitextra1, newdata = ndat2, type="response", onlySL = T)$pred); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_3*delta + 1-ydat$L1_3) + ypred*(1-delta)*(ydat$L1_3)  
  ydat$ypred = ifelse(ydat$Y_3==0,0,ydat$ypred); 
  
  l = as.data.frame(cbind("L1_2" = ydat$L1_2, "L2_2" = ydat$L2_2, "W1" = ydat$W1, "W2" = ydat$W2))
  A_2 = ydat$A_2; y = ydat$ypred
  tmp = ydat[,c("L1_2","L2_2","W1","W2","A_2")]; tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitextra2, newdata = tmp, onlySL = T)$pred) 
  yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = (fintc2)/(pi2*pi2c), data = ydat) 
  ydat = tmp_test[!is.na(tmp_test$C_1) & tmp_test$C_1==0,]; 
  ndat = ydat; ndat = ndat[,c("L1_2","L2_2","W1","W2","A_2")]
  ndat[is.na(ndat)] = 0; ndat$A_2 = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitextra2, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response");   
  ndat2 = ydat; ndat2 = ndat2[,c("L1_2","L2_2","W1","W2","A_2")]; ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitextra2, newdata = ndat2, type="response", onlySL = T)$pred); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_2*delta + 1-ydat$L1_2) + ypred*(1-delta)*(ydat$L1_2)
  ydat$ypred = ifelse(ydat$Y_2==0,0,ydat$ypred); 
  
  l = as.data.frame(cbind("L1_1" = ydat$L1_1, "L2_1" = ydat$L2_1, "W1" = ydat$W1, "W2" = ydat$W2))
  A_1 = ydat$A_1; y = ydat$ypred
  tmp = ydat[,c("L1_1","L2_1","W1","W2","A_1")]; tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitextra3, newdata = tmp, onlySL = T)$pred) 
  yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = (fintc1)/(pi1*pi1c), data = ydat)
  ydat = tmp_test[!is.na(tmp_test$C_0) & tmp_test$C_0==0,]; 
  ndat = ydat; ndat = ndat[,c("L1_1","L2_1","W1","W2","A_1")]
  ndat[is.na(ndat)] = 0; ndat$A_1 = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitextra3, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response");   
  ndat2 = ydat; ndat2 = ndat2[,c("L1_1","L2_1","W1","W2","A_1")]; ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitextra3, newdata = ndat2, type="response", onlySL = T)$pred); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_1*delta + 1-ydat$L1_1) + ypred*(1-delta)*(ydat$L1_1)
  ydat$ypred = ifelse(ydat$Y_1==0,0,ydat$ypred); 
  
  l = as.data.frame(cbind("L1_0" = ydat$L1_0, "L2_0" = ydat$L2_0, "W1" = ydat$W1, "W2" = ydat$W2))
  A_0 = ydat$A_0; y = ydat$ypred
  tmp = ydat[,c("L1_0","L2_0","W1","W2","A_0")]; tmp[is.na(tmp)] = 0
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitextra4, newdata = tmp, onlySL = T)$pred) 
  yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = (fintc0)/(pi0*pi0c), data = ydat)
  ydat = tmp_test; 
  ndat = ydat; ndat = ndat[,c("L1_0","L2_0","W1","W2","A_0")]
  ndat[is.na(ndat)] = 0; ndat$A_0 = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitextra4, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response");   
  ndat2 = ydat; ndat2 = ndat2[,c("L1_0","L2_0","W1","W2","A_0")]; ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitextra4, newdata = ndat2, type="response", onlySL = T)$pred); 
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(ydat$L1_0*delta + 1-ydat$L1_0) + ypred*(1-delta)*(ydat$L1_0)
  
  meany12_2 = mean(ydat$ypred)
  
  meany12 = (meany12_1 + meany12_2)/2
  
  myparam = cbind(meany12)
  gc()
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"dr_ss_param.csv")

stopCluster(cl)