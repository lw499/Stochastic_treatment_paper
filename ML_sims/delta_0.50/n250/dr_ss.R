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
  library(geepack);library(MASS);library(ResourceSelection);library(ltmle); library(SuperLearner)
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
  set.seed(987654)
  seeds = floor(runif(1000)*10^8);
  set.seed(seeds[m])
  
  n <- 250
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
  #learners <- create.Learner("SL.ranger", tune = list(mtry = 5, num.trees = 500))  
  #learners <- create.Learner("SL.rpart", tune = list(cp = c(0.025,0.05)))  
  #deg.gam = c(4,5,6)
  #learners <- create.Learner("SL.gam", tune = list(deg.gam = c(1,2,3,4,5,6,7,8)))
  sl.lib = c("SL.gam","SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.glm.interaction","SL.step","SL.nnet","SL.earth","SL.ranger")
  
  #sl.lib = c("SL.gam","SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.glm.interaction","SL.step.forward","SL.ranger")
  
  set.seed(100388)
  train_ind = sample(unique(dffull$id),size = n/2)
  
  #fit afit and cfit:
  tmpdffull  = dffull[dffull$id %in% train_ind,]
  tmpdffull = tmpdffull[!is.na(tmpdffull$A) & tmpdffull$lag_A==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  y = tmpdffull$A
  set.seed(1003)
  afit = SuperLearner(Y=y, X=cbind(l), SL.library=sl.lib, family=binomial)
  
  tmpdffull  = dffull[dffull$id %in% train_ind,]
  tmpdffull = tmpdffull[!is.na(tmpdffull$C) & tmpdffull$lag_C==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "A" = tmpdffull$A, "t0" = tmpdffull$t0))
  y = tmpdffull$C
  set.seed(1003)
  cfit = SuperLearner(Y=y, X=cbind(l), SL.library=sl.lib, family=binomial)
  
  #fit afit and cfit (and then swap):
  tmpdffull  = dffull[!dffull$id %in% train_ind,]
  tmpdffull = tmpdffull[!is.na(tmpdffull$A) & tmpdffull$lag_A==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  y = tmpdffull$A
  set.seed(1003)
  afitb = SuperLearner(Y=y, X=cbind(l), SL.library=sl.lib, family=binomial)
  
  tmpdffull  = dffull[!dffull$id %in% train_ind,]
  tmpdffull = tmpdffull[!is.na(tmpdffull$C) & tmpdffull$lag_C==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "A" = tmpdffull$A, "t0" = tmpdffull$t0))
  y = tmpdffull$C
  set.seed(1003)
  cfitb = SuperLearner(Y=y, X=cbind(l), SL.library=sl.lib, family=binomial)
  
  #fit hazard
  tmpdffull  = dffull[dffull$id %in% train_ind,]
  tmpdffull = tmpdffull[!is.na(tmpdffull$Y) & tmpdffull$C==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  A = tmpdffull$A; y = tmpdffull$Y
  set.seed(1003)
  yfitog = SuperLearner(Y=y, X=cbind(l,A), SL.library=sl.lib, family=binomial)
  #fit hazard (swap train and test)
  tmpdffull  = dffull[!dffull$id %in% train_ind,]
  tmpdffull = tmpdffull[!is.na(tmpdffull$Y) & tmpdffull$C==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  A = tmpdffull$A; y = tmpdffull$Y
  set.seed(1003)
  yfitogb = SuperLearner(Y=y, X=cbind(l,A), SL.library=sl.lib, family=binomial)
  
  
  #iteration 1 (this is only for computationl simplicity. 
  #Even though fit for everyone, we will eventually only use those in training):
  dffulltmp = dffull[,c("L1","L2","W1","W2","t0")]
  dffull$pred_obs = predict.SuperLearner(afit, newdata = dffulltmp, type="response", onlySL = T)$pred
  dffull$pred_obs = ifelse(dffull$A==1, dffull$pred_obs, 1-dffull$pred_obs)
  dffull$pred_obs = ifelse(!is.na(dffull$lag_A) & dffull$lag_A==1, 1, dffull$pred_obs)
  dffulltmp = dffull[,c("L1","L2","W1","W2","A","t0")]
  dffull$pred_obsc = predict.SuperLearner(cfit, newdata = dffulltmp, type="response", onlySL = T)$pred
  dffull$pred_obsc = ifelse(dffull$C==1, NA, 1-dffull$pred_obsc)
  dffull$pred_obsc = ifelse(!is.na(dffull$lag_C) & dffull$lag_C==1, 1, dffull$pred_obsc)
  # for q^g for all time points
  dffull$fint = (dffull$L1*delta + 1 - dffull$L1)*dffull$pred_obs + dffull$L1*dffull$A*(1-delta)
  
  dffullwide = dcast(dffull, id + W1 + W2 ~ t0, value.var = c("L1","L2","A","C","Y","pred_obs","pred_obsc","fint"))
  
  tmpdata = dffullwide
  
  tmpindC = NULL; tmpindA = NULL; tmpindY = NULL;
  tmppredA = NULL; tmppredC = NULL
  tmpfintc  = NULL;
  for(j in 1:(K)) #dwide[,c(paste("overallcen","_",4, sep=""))]
  {
    tmpindC = c(tmpindC, paste("C","_",j-1, sep=""))
    tmpindA = c(tmpindA, paste("A","_",j-1, sep=""))
    tmpindY = c(tmpindY, paste("Y","_",j-1, sep=""))
    tmppredA = c(tmppredA, paste("pi","_",j-1, sep=""))
    tmppredC = c(tmppredC, paste("pic","_",j-1, sep=""))
    tmpfintc = c(tmpfintc, paste("fintc","_",j-1, sep=""))
    tmpdata[ , tmppredA] <- NA
    tmpdata[ , tmppredC] <- NA
    tmpdata[ , tmpfintc] <- NA
  }
  
  tmpdata$pi_0 = tmpdata$pred_obs_0; tmpdata$pic_0 = tmpdata$pred_obsc_0;
  tmpdata$fintc_0 = tmpdata$fint_0
  for(i in 1:(K-1))
  {
    tmppi = select(tmpdata, c(paste("pi","_",i-1, sep="")))[[1]] * select(tmpdata, c(paste("pred_obs","_",i, sep="")))[[1]]
    set(tmpdata,NULL , paste0("pi_",i), tmppi)
    tmppic = select(tmpdata, c(paste("pic","_",i-1, sep="")))[[1]] * select(tmpdata, c(paste("pred_obsc","_",i, sep="")))[[1]]
    set(tmpdata,NULL , paste0("pic_",i), tmppic)
    
    lastfintc = select(tmpdata, c(paste("fintc","_",i-1, sep="")))[[1]]
    currentfint = select(tmpdata, c(paste("fint","_",i, sep="")))[[1]]
    tmpfint = currentfint * lastfintc
    set(tmpdata,NULL , paste0("fintc_",i), tmpfint)
  }
  
  #learners <- create.Learner("SL.gam", tune = list(deg.gam = c(0,1,3)))
  #learners <- create.Learner("SL.ranger", tune = list(mtry = 6, num.trees = 1000, min.node.size=10))  
  sl.lib = c("SL.gam","SL.glm", "SL.bayesglm", "SL.glm.interaction","SL.ranger")
  
  ##################
  ######time 5######
  ##################
  tmp_train =  tmpdata[tmpdata$id %in% train_ind,]
  mydata = tmp_train %>% filter_(paste0(c("C_"),c(3), "==",0)  ) 
  ydat = mydata; ydat = as.data.table(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4,  "W1" = ydat$W1, "W2" = ydat$W2, "t0" = K-1))
  ydat[is.na(ydat)] = 0;  ydat$A = 1;
  set.seed(1003); ypred = predict.SuperLearner(yfitog, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = mydata; ydat2 = as.data.table(cbind("L1" = ydat2$L1_4, "L2" = ydat2$L2_4, "W1" = ydat2$W1, "W2" = ydat2$W2, "t0" = K-1, "A" = ydat2$A_4))
  ydat2[is.na(ydat2)] = 0
  Ltmp = unlist(ydat2 %>% select_at(paste0(c("L1")) ))
  mydata$ypred = predict.SuperLearner(yfitog, newdata = ydat2, type="response", onlySL = T)$pred*(Ltmp*delta + 1-Ltmp) + ypred*(1-delta)*(Ltmp)
  Ytmp = unlist(mydata %>% select_at(paste0(c("Y_"),c(3)) ))
  mydata$ypred = ifelse(Ytmp==0,0,mydata$ypred); ## This is really Y_4
  
  res_list <- vector(mode="list", length=K-1)
  
  #m=K-2
  for(m in (K-2):0)
  {
    Lformulatmp = c(paste0(c("L1_", "L2_"),c(m,m)),"W1", "W2")
    l = mydata %>% select_at(Lformulatmp)
    a = unlist(mydata %>% select_at(paste0(c("A_"),c(m)))); y = mydata$ypred
    set.seed(1003); res_list[[m+1]] = SuperLearner(Y=y, X=cbind(l,a), SL.library=sl.lib, family=quasibinomial (link = "logit"))
    if(m>0) {mydata = tmp_train %>% filter_(paste0(c("C_"),c(m-1), "==",0)  )} else{mydata = tmp_train}  #C_2
    ydat = mydata %>% select_at(Lformulatmp); ydat$a=1;
    ydat[is.na(ydat)] = 0; 
    ypred = predict.SuperLearner(res_list[[m+1]], newdata = ydat, type="response", onlySL = T)$pred; 
    ydat2 = mydata %>% select_at(Lformulatmp); 
    ydat2$a = unlist(mydata %>% select_at(paste0(c("A_"),c(m)) )); ydat2[is.na(ydat2)] = 0; 
    Ltmp = unlist(mydata %>% select_at(paste0(c("L1_"),c(m)) ))
    mydata$ypred = predict.SuperLearner(res_list[[m+1]], newdata = ydat2, type="response", onlySL = T)$pred*(Ltmp*delta + 1-Ltmp) + ypred*(1-delta)*Ltmp
    if(m>0){Ytmp = as.vector(unlist(mydata %>% select_at(paste0(c("Y_"),c(m-1)) )));
    mydata$ypred = ifelse(Ytmp==0,0,mydata$ypred);}
  }
  
  ### Predict and refit on test dataset
  ##################
  ######time 5######
  ##################
  tmp_test =  tmpdata[!tmpdata$id %in% train_ind,]
  ydat = tmp_test %>% filter_(paste0(c("C_"), c(K-1), "==",0)  ) 
  tmp = as.data.table(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4, "W1" = ydat$W1, "W2" = ydat$W2, "t0" = K-1, "A" = ydat$A_4)); tmp[is.na(tmp)] = 0; 
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitog, newdata = tmp, onlySL = T)$pred)
  Ytmp = unlist(ydat %>% select_at(paste0(c("Y_"),c(K-1)) ))
  denom = unlist(ydat %>% select_at(paste0(c("pi_"),c(K-1)) )) * unlist(ydat %>% select_at(paste0(c("pic_"),c(K-1)) ))
  weight = unlist(ydat %>% select_at(paste0(c("fintc_"),c(K-1)) ))/denom
  yfit2 = glm2(Ytmp ~ 1, family = binomial(), offset = ypred0, weights = weight, data = ydat) 
  ###
  ydat = tmp_test %>% filter_(paste0(c("C_"),c(K-2), "==",0) )
  ndat = ydat; ndat = as.data.frame(cbind("L1" = ndat$L1_4, "L2" = ndat$L2_4, "W1" = ndat$W1, "W2" = ndat$W2, "t0" = K-1, "A" = ndat$A_4));
  ndat[is.na(ndat)] = 0;  ndat$A = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitog, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response"); 
  ndat2 = ydat; ndat2 = as.data.frame(cbind("L1" = ndat2$L1_4, "L2" = ndat2$L2_4, "W1" = ndat2$W1, "W2" = ndat2$W2, "t0" = K-1, "A" = ndat2$A_4)); ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitog, newdata = ndat2, type="response", onlySL = T)$pred); 
  Ltmp = unlist(ydat %>% select_at(paste0(c("L1_"),c(4)) ))
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(Ltmp*delta + 1-Ltmp) + ypred*(1-delta)*(Ltmp)
  Ytmp = unlist(ydat %>% select_at(paste0(c("Y_"),c(3)) ))
  ydat$ypred = ifelse(Ytmp==0,0,ydat$ypred); ## This is really Y_4
  
  #m=K-2
  for(m in (K-2):0)
  {
    formulatmp = c(paste0(c("L1_", "L2_"),c(m,m)) ,"W1", "W2")
    tmp = ydat %>% select_at(formulatmp)
    Aformulatmp = paste0(c("A_"),c(m)); tmp$a = unlist(ydat %>% select_at(Aformulatmp)); tmp[is.na(tmp)] = 0
    ydat$ypred0 = qlogis(predict.SuperLearner(res_list[[m+1]], newdata = tmp, onlySL = T)$pred)
    denom = unlist(ydat %>% select_at(paste0(c("pi_"),c(m)) )) * unlist(ydat %>% select_at(paste0(c("pic_"),c(m)) ))
    weight = unlist(ydat %>% select_at(paste0(c("fintc_"),c(m)) ))/denom
    yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = weight, data = ydat) 
    if(m>0) {ydat = tmp_test %>% filter_(paste0(c("C_"),c(m-1), "==",0) )} else{ydat = tmp_test} 
    ndat = ydat %>% select_at(formulatmp); ndat$a=1; ndat[is.na(ndat)] = 0; 
    ndat$ypred0 = qlogis(predict.SuperLearner(res_list[[m+1]], newdata = ndat, type="response", onlySL = T)$pred); 
    ypred = predict(yfit2, newdata = ndat, type="response");   
    ndat2 = ydat %>% select_at(formulatmp); 
    Aformulatmp = paste0(c("A_"),c(m)); ndat2$a = unlist(ydat %>% select_at(Aformulatmp)); ndat2[is.na(ndat2)] = 0
    ndat2$ypred0 = qlogis(predict.SuperLearner(res_list[[m+1]], newdata = ndat2, type="response", onlySL = T)$pred); 
    Ltmp = unlist(ydat %>% select_at(paste0(c("L1_"),c(m)) ))
    ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(Ltmp*delta + 1-Ltmp) + ypred*(1-delta)*(Ltmp)
    if(m>0){Ytmp = unlist(ydat %>% select_at(paste0(c("Y_"),c(m-1)) ));
    ydat$ypred = ifelse(Ytmp==0,0,ydat$ypred);}
  }
  
  meany12_1 = mean(ydat$ypred)
  
  ## REPEAT (swap the train and test dataset)
  dffulltmp = dffull[,c("L1","L2","W1","W2","t0")]
  dffull$pred_obs = predict.SuperLearner(afitb, newdata = dffulltmp, type="response", onlySL = T)$pred
  dffull$pred_obs = ifelse(dffull$A==1, dffull$pred_obs, 1-dffull$pred_obs)
  dffull$pred_obs = ifelse(!is.na(dffull$lag_A) & dffull$lag_A==1, 1, dffull$pred_obs)
  dffulltmp = dffull[,c("L1","L2","W1","W2","A","t0")]
  dffull$pred_obsc = predict.SuperLearner(cfitb, newdata = dffulltmp, type="response", onlySL = T)$pred
  dffull$pred_obsc = ifelse(dffull$C==1, NA, 1-dffull$pred_obsc)
  dffull$pred_obsc = ifelse(!is.na(dffull$lag_C) & dffull$lag_C==1, 1, dffull$pred_obsc)
  
  dffull$fint = (dffull$L1*delta + 1 - dffull$L1)*dffull$pred_obs + dffull$L1*dffull$A*(1-delta)
  
  dffullwide = dcast(dffull, id + W1 + W2 ~ t0, value.var = c("L1","L2","A","C","Y","pred_obs","pred_obsc","fint"))
  
  tmpdata = dffullwide
  
  tmpindC = NULL; tmpindA = NULL; tmpindY = NULL;
  tmppredA = NULL; tmppredC = NULL
  tmpfintc  = NULL;
  for(j in 1:(K)) #dwide[,c(paste("overallcen","_",4, sep=""))]
  {
    tmpindC = c(tmpindC, paste("C","_",j-1, sep=""))
    tmpindA = c(tmpindA, paste("A","_",j-1, sep=""))
    tmpindY = c(tmpindY, paste("Y","_",j-1, sep=""))
    tmppredA = c(tmppredA, paste("pi","_",j-1, sep=""))
    tmppredC = c(tmppredC, paste("pic","_",j-1, sep=""))
    tmpfintc = c(tmpfintc, paste("fintc","_",j-1, sep=""))
    tmpdata[ , tmppredA] <- NA
    tmpdata[ , tmppredC] <- NA
    tmpdata[ , tmpfintc] <- NA
  }
  
  tmpdata$pi_0 = tmpdata$pred_obs_0; tmpdata$pic_0 = tmpdata$pred_obsc_0;
  tmpdata$fintc_0 = tmpdata$fint_0
  for(i in 1:(K-1))
  {
    tmppi = select(tmpdata, c(paste("pi","_",i-1, sep="")))[[1]] * select(tmpdata, c(paste("pred_obs","_",i, sep="")))[[1]]
    set(tmpdata,NULL , paste0("pi_",i), tmppi)
    tmppic = select(tmpdata, c(paste("pic","_",i-1, sep="")))[[1]] * select(tmpdata, c(paste("pred_obsc","_",i, sep="")))[[1]]
    set(tmpdata,NULL , paste0("pic_",i), tmppic)
    
    lastfintc = select(tmpdata, c(paste("fintc","_",i-1, sep="")))[[1]]
    currentfint = select(tmpdata, c(paste("fint","_",i, sep="")))[[1]]
    tmpfint = currentfint * lastfintc
    set(tmpdata,NULL , paste0("fintc_",i), tmpfint)
  }
  
  ##################
  ######time 5######
  ##################
  tmp_train =  tmpdata[!tmpdata$id %in% train_ind,]
  mydata = tmp_train %>% filter_(paste0(c("C_"),c(3), "==",0)  ) 
  ydat = mydata; ydat = as.data.table(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4,  "W1" = ydat$W1, "W2" = ydat$W2, "t0" = K-1))
  ydat[is.na(ydat)] = 0;  ydat$A = 1;
  set.seed(1003); ypred = predict.SuperLearner(yfitogb, newdata = ydat, type="response", onlySL = T)$pred; 
  ydat2 = mydata; ydat2 = as.data.table(cbind("L1" = ydat2$L1_4, "L2" = ydat2$L2_4, "W1" = ydat2$W1, "W2" = ydat2$W2, "t0" = K-1, "A" = ydat2$A_4))
  ydat2[is.na(ydat2)] = 0
  Ltmp = unlist(ydat2 %>% select_at(paste0(c("L1")) ))
  mydata$ypred = predict.SuperLearner(yfitogb, newdata = ydat2, type="response", onlySL = T)$pred*(Ltmp*delta + 1-Ltmp) + ypred*(1-delta)*(Ltmp)
  Ytmp = unlist(mydata %>% select_at(paste0(c("Y_"),c(3)) ))
  mydata$ypred = ifelse(Ytmp==0,0,mydata$ypred); ## This is really Y_4
  
  res_list <- vector(mode="list", length=K-1)
  
  #m=K-2
  for(m in (K-2):0)
  {
    Lformulatmp = c(paste0(c("L1_", "L2_"),c(m,m)),"W1", "W2")
    l = mydata %>% select_at(Lformulatmp)
    a = unlist(mydata %>% select_at(paste0(c("A_"),c(m)))); y = mydata$ypred
    set.seed(1003); res_list[[m+1]] = SuperLearner(Y=y, X=cbind(l,a), SL.library=sl.lib, family=quasibinomial (link = "logit") )
    if(m>0) {mydata = tmp_train %>% filter_(paste0(c("C_"),c(m-1), "==",0)  )} else{mydata = tmp_train}  #C_2
    ydat = mydata %>% select_at(Lformulatmp); ydat$a=1;
    ydat[is.na(ydat)] = 0; 
    ypred = predict.SuperLearner(res_list[[m+1]], newdata = ydat, type="response", onlySL = T)$pred; 
    ydat2 = mydata %>% select_at(Lformulatmp); 
    ydat2$a = unlist(mydata %>% select_at(paste0(c("A_"),c(m)) )); ydat2[is.na(ydat2)] = 0; 
    Ltmp = unlist(mydata %>% select_at(paste0(c("L1_"),c(m)) ))
    mydata$ypred = predict.SuperLearner(res_list[[m+1]], newdata = ydat2, type="response", onlySL = T)$pred*(Ltmp*delta + 1-Ltmp) + ypred*(1-delta)*Ltmp
    if(m>0){Ytmp = as.vector(unlist(mydata %>% select_at(paste0(c("Y_"),c(m-1)) )));
    mydata$ypred = ifelse(Ytmp==0,0,mydata$ypred);}
  }
  
  ### Predict and refit on test dataset
  ##################
  ######time 5######
  ##################
  tmp_test =  tmpdata[tmpdata$id %in% train_ind,]
  ydat = tmp_test %>% filter_(paste0(c("C_"), c(K-1), "==",0)  ) 
  tmp = as.data.table(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4, "W1" = ydat$W1, "W2" = ydat$W2, "t0" = K-1, "A" = ydat$A_4)); tmp[is.na(tmp)] = 0; 
  ydat$ypred0 = qlogis(predict.SuperLearner(yfitogb, newdata = tmp, onlySL = T)$pred)
  Ytmp = unlist(ydat %>% select_at(paste0(c("Y_"),c(K-1)) ))
  denom = unlist(ydat %>% select_at(paste0(c("pi_"),c(K-1)) )) * unlist(ydat %>% select_at(paste0(c("pic_"),c(K-1)) ))
  weight = unlist(ydat %>% select_at(paste0(c("fintc_"),c(K-1)) ))/denom
  yfit2 = glm2(Ytmp ~ 1, family = binomial(), offset = ypred0, weights = weight, data = ydat) 
  ###
  ydat = tmp_test %>% filter_(paste0(c("C_"),c(K-2), "==",0) )
  ndat = ydat; ndat = as.data.frame(cbind("L1" = ndat$L1_4, "L2" = ndat$L2_4, "W1" = ndat$W1, "W2" = ndat$W2, "t0" = K-1, "A" = ndat$A_4));
  ndat[is.na(ndat)] = 0;  ndat$A = 1;
  ndat$ypred0 = qlogis(predict.SuperLearner(yfitogb, newdata = ndat, type="response", onlySL = T)$pred); 
  ypred = predict(yfit2, newdata = ndat, type="response"); 
  ndat2 = ydat; ndat2 = as.data.frame(cbind("L1" = ndat2$L1_4, "L2" = ndat2$L2_4, "W1" = ndat2$W1, "W2" = ndat2$W2, "t0" = K-1, "A" = ndat2$A_4)); ndat2[is.na(ndat2)] = 0
  ndat2$ypred0 = qlogis(predict.SuperLearner(yfitogb, newdata = ndat2, type="response", onlySL = T)$pred); 
  Ltmp = unlist(ydat %>% select_at(paste0(c("L1_"),c(4)) ))
  ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(Ltmp*delta + 1-Ltmp) + ypred*(1-delta)*(Ltmp)
  Ytmp = unlist(ydat %>% select_at(paste0(c("Y_"),c(3)) ))
  ydat$ypred = ifelse(Ytmp==0,0,ydat$ypred); ## This is really Y_4
  
  #m=K-2
  for(m in (K-2):0)
  {
    formulatmp = c(paste0(c("L1_", "L2_"),c(m,m)) ,"W1", "W2")
    tmp = ydat %>% select_at(formulatmp)
    Aformulatmp = paste0(c("A_"),c(m)); tmp$a = unlist(ydat %>% select_at(Aformulatmp)); tmp[is.na(tmp)] = 0
    ydat$ypred0 = qlogis(predict.SuperLearner(res_list[[m+1]], newdata = tmp, onlySL = T)$pred)
    denom = unlist(ydat %>% select_at(paste0(c("pi_"),c(m)) )) * unlist(ydat %>% select_at(paste0(c("pic_"),c(m)) ))
    weight = unlist(ydat %>% select_at(paste0(c("fintc_"),c(m)) ))/denom
    yfit2 = glm2(ypred ~ 1, family = binomial(), offset = ypred0, weights = weight, data = ydat) 
    if(m>0) {ydat = tmp_test %>% filter_(paste0(c("C_"),c(m-1), "==",0) )} else{ydat = tmp_test} 
    ndat = ydat %>% select_at(formulatmp); ndat$a=1; ndat[is.na(ndat)] = 0; 
    ndat$ypred0 = qlogis(predict.SuperLearner(res_list[[m+1]], newdata = ndat, type="response", onlySL = T)$pred); 
    ypred = predict(yfit2, newdata = ndat, type="response");   
    ndat2 = ydat %>% select_at(formulatmp); 
    Aformulatmp = paste0(c("A_"),c(m)); ndat2$a = unlist(ydat %>% select_at(Aformulatmp)); ndat2[is.na(ndat2)] = 0
    ndat2$ypred0 = qlogis(predict.SuperLearner(res_list[[m+1]], newdata = ndat2, type="response", onlySL = T)$pred); 
    Ltmp = unlist(ydat %>% select_at(paste0(c("L1_"),c(m)) ))
    ydat$ypred = predict(yfit2, newdata = ndat2, type="response")*(Ltmp*delta + 1-Ltmp) + ypred*(1-delta)*(Ltmp)
    if(m>0){Ytmp = unlist(ydat %>% select_at(paste0(c("Y_"),c(m-1)) ));
    ydat$ypred = ifelse(Ytmp==0,0,ydat$ypred);}
  }
  
  meany12_2 = mean(ydat$ypred)
  
  meany12 = (meany12_1 + meany12_2)/2
  
  myparam = cbind(meany12)
  gc()
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"dr_ss.csv")

stopCluster(cl)

