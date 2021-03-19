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
  
  sl.lib = c("SL.gam","SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.glm.interaction","SL.step","SL.nnet","SL.earth","SL.ranger")
  
  #fit hazard
  tmpdffull = dffull[!is.na(dffull$Y) & dffull$C==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  A = tmpdffull$A; y = tmpdffull$Y
  set.seed(1003)
  yfitog = SuperLearner(Y=y, X=cbind(l,A), SL.library=sl.lib, family=binomial)
  
  dffullwide = dcast(dffull, id + W1 + W2 ~ t0, value.var = c("L1","L2","A","C","Y"))
  
  tmpdata = dffullwide
  #subset data
  tmpdata$Y_5 = tmpdata$Y_4
  tmpdata$Y_4 = tmpdata$Y_3
  tmpdata$Y_3 = tmpdata$Y_2
  tmpdata$Y_2 = tmpdata$Y_1
  tmpdata$Y_1 = tmpdata$Y_0
  tmpdata$Y_0 = NULL
  
  tmpdata$id = seq(1,n,by=1)

  sl.lib = c("SL.gam","SL.glm", "SL.bayesglm", "SL.glm.interaction","SL.ranger")
  ##################
  ######time 5######
  ##################
  y4dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0,]; 
    ydat = y4dat; ydat = as.data.table(cbind("L1" = ydat$L1_4, "L2" = ydat$L2_4,  "W1" = ydat$W1, "W2" = ydat$W2, "t0" = 4))
    ydat[is.na(ydat)] = 0;  ydat$A = 1;
    ypred = predict.SuperLearner(yfitog, newdata = ydat, type="response",onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = as.data.table(cbind("L1" = ydat2$L1_4, "L2" = ydat2$L2_4, "A" = ydat2$A_4, "W1" = ydat2$W1, "W2" = ydat2$W2, "t0" = 4))
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitog, newdata = ydat2, type="response",onlySL = T)$pred*(y4dat$L1_4*delta + 1-y4dat$L1_4) + ypred*(1-delta)*(y4dat$L1_4)
  y4dat$y4pred = ifelse(y4dat$Y_4==0,0,y4dat$y4pred); 
  
  #sl.lib = c("SL.earth","SL.gam","SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.glm.interaction")
  l = as.data.frame(cbind("L1_3" = y4dat$L1_3, "L2_3" = y4dat$L2_3, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_3 = y4dat$A_3; y = y4dat$y4pred
  set.seed(1003); yfitextra = SuperLearner(Y=y, X=cbind(l,A_3), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0,]; 
    ydat = y4dat; ydat = ydat[,c("L1_3","L2_3","A_3","W1","W2")]
    ydat[is.na(ydat)] = 0; ydat$A_3 = 1;
    ypred = predict.SuperLearner(yfitextra, newdata = ydat, type="response",onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_3","L2_3","A_3","W1","W2")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra, newdata = ydat2, type="response",onlySL = T)$pred*(y4dat$L1_3*delta + 1-y4dat$L1_3) + ypred*(1-delta)*(y4dat$L1_3)
  y4dat$y4pred = ifelse(y4dat$Y_3==0,0,y4dat$y4pred); 

  l = as.data.frame(cbind("L1_2" = y4dat$L1_2, "L2_2" = y4dat$L2_2, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_2 = y4dat$A_2; y = y4dat$y4pred
  set.seed(1003); yfitextra = SuperLearner(Y=y, X=cbind(l,A_2), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0,]; 
    ydat = y4dat; ydat = ydat[,c("L1_2","L2_2","A_2","W1","W2")]
    ydat[is.na(ydat)] = 0; ydat$A_2 = 1;
    ypred = predict.SuperLearner(yfitextra, newdata = ydat, type="response",onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_2","L2_2","A_2","W1","W2")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra, newdata = ydat2, type="response",onlySL = T)$pred*(y4dat$L1_2*delta + 1-y4dat$L1_2) + ypred*(1-delta)*(y4dat$L1_2)
  y4dat$y4pred = ifelse(y4dat$Y_2==0,0,y4dat$y4pred); 

  
  l = as.data.frame(cbind("L1_1" = y4dat$L1_1, "L2_1" = y4dat$L2_1, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_1 = y4dat$A_1; y = y4dat$y4pred
  set.seed(1003); yfitextra = SuperLearner(Y=y, X=cbind(l,A_1), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmpdata[!is.na(tmpdata$C_0) & tmpdata$C_0==0,]; 
  ydat = y4dat; ydat = ydat[,c("L1_1","L2_1","A_1","W1","W2")]
  ydat[is.na(ydat)] = 0; ydat$A_1 = 1;
  ypred = predict.SuperLearner(yfitextra, newdata = ydat, type="response",onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_1","L2_1","A_1","W1","W2")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra, newdata = ydat2, type="response",onlySL = T)$pred*(y4dat$L1_1*delta + 1-y4dat$L1_1) + ypred*(1-delta)*(y4dat$L1_1)
  y4dat$y4pred = ifelse(y4dat$Y_1==0,0,y4dat$y4pred); 
  
  l = as.data.frame(cbind("L1_0" = y4dat$L1_0, "L2_0" = y4dat$L2_0, "W1" = y4dat$W1, "W2" = y4dat$W2))
  A_0 = y4dat$A_0; y = y4dat$y4pred
  set.seed(1003); yfitextra = SuperLearner(Y=y, X=cbind(l,A_0), SL.library=sl.lib, family=quasibinomial)
  y4dat = tmpdata; 
  ydat = y4dat; ydat = ydat[,c("L1_0","L2_0","A_0","W1","W2")]
  ydat[is.na(ydat)] = 0; ydat$A_0 = 1;
  ypred = predict.SuperLearner(yfitextra, newdata = ydat, type="response",onlySL = T)$pred; 
  ydat2 = y4dat; ydat2 = ydat2[,c("L1_0","L2_0","A_0","W1","W2")]
  ydat2[is.na(ydat2)] = 0
  y4dat$y4pred = predict.SuperLearner(yfitextra, newdata = ydat2, type="response",onlySL = T)$pred*(y4dat$L1_0*delta + 1-y4dat$L1_0) + ypred*(1-delta)*(y4dat$L1_0)
  
    
  meany4tmp = c(mean(y4dat$y4pred))
  
  meany5 = (meany4tmp)
  #meany5
  
  myparam = cbind(meany5)
  
  gc()
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"seq_reg.csv")

stopCluster(cl)