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
  options(warn=-1)
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
  set.seed(123456)
  seeds = floor(runif(1000)*10^8);
  set.seed(seeds[m])
  
  n <- 2500
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
  
  sl.lib = c("SL.gam","SL.glm", "SL.stepAIC", "SL.bayesglm", "SL.glm.interaction","SL.step","SL.nnet","SL.earth","SL.ranger")
  
  tmpdffull = dffull[!is.na(dffull$A) & dffull$lag_A==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "t0" = tmpdffull$t0))
  y = tmpdffull$A
  set.seed(1003)
  afit = SuperLearner(Y=y, X=cbind(l), SL.library=sl.lib, family=binomial)
  
  tmpdffull = dffull[!is.na(dffull$C) & dffull$lag_C==0,]
  l = as.data.frame(cbind("L1" = tmpdffull$L1, "L2" = tmpdffull$L2, "W1" = tmpdffull$W1, "W2" = tmpdffull$W2, "A" = tmpdffull$A, "t0" = tmpdffull$t0))
  y = tmpdffull$C
  set.seed(1003)
  cfit = SuperLearner(Y=y, X=cbind(l), SL.library=sl.lib, family=binomial)
  
  dffulltmp = dffull[,c("L1","L2","W1","W2","t0")]
  dffull$pred_obs = predict.SuperLearner(afit, newdata = dffulltmp, type="response", onlySL = T)$pred
  dffull$pred_obs = ifelse(dffull$A==1, dffull$pred_obs, 1-dffull$pred_obs)
  dffull$pred_obs = ifelse(!is.na(dffull$lag_A) & dffull$lag_A==1, 1, dffull$pred_obs)
  dffulltmp = dffull[,c("L1","L2","W1","W2","A","t0")]
  dffull$pred_obsc = predict.SuperLearner(cfit, newdata = dffulltmp, type="response", onlySL = T)$pred
  dffull$pred_obsc = ifelse(dffull$C==1, NA, 1-dffull$pred_obsc)
  dffull$pred_obsc = ifelse(!is.na(dffull$lag_C) & dffull$lag_C==1, 1, dffull$pred_obsc)
  
  dffull$lag_A =  dffull$lag_C = NULL
  dffullwide = dcast(dffull, id + W1 + W2 ~ t0, value.var = c("L1","L2","A","C","Y","pred_obs","pred_obsc"))
  
  tmpdata = dffullwide
  tmpdata$Y_5 = tmpdata$Y_4
  tmpdata$Y_4 = tmpdata$Y_3
  tmpdata$Y_3 = tmpdata$Y_2
  tmpdata$Y_2 = tmpdata$Y_1
  tmpdata$Y_1 = tmpdata$Y_0
  tmpdata$Y_0 = NULL
  
  tmpdata$id = seq(1,n,by=1)
  tmpdata$pi4 <- tmpdata$pi3 <- tmpdata$pi2 <- tmpdata$pi1 <- tmpdata$pi0 <- NA
  tmpdata$pi4c <- tmpdata$pi3c <- tmpdata$pi2c <- tmpdata$pi1c <- tmpdata$pi0c <- NA
  
  tmpdata$pi0 = tmpdata$pred_obs_0
  tmpdata$pi1 = tmpdata$pi0*tmpdata$pred_obs_1
  tmpdata$pi2 = tmpdata$pi1*tmpdata$pred_obs_2
  tmpdata$pi3 = tmpdata$pi2*tmpdata$pred_obs_3
  tmpdata$pi4 = tmpdata$pi3*tmpdata$pred_obs_4
  
  tmpdata$pi0c = tmpdata$pred_obsc_0
  tmpdata$pi1c = tmpdata$pi0c*tmpdata$pred_obsc_1
  tmpdata$pi2c = tmpdata$pi1c*tmpdata$pred_obsc_2
  tmpdata$pi3c = tmpdata$pi2c*tmpdata$pred_obsc_3
  tmpdata$pi4c = tmpdata$pi3c*tmpdata$pred_obsc_4
  
  ##calculate risk
  mean = NULL
  ind = NA;
  #time 1
  tmp = tmpdata[!is.na(tmpdata$Y_1),]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  fit1 = glm2(Y_1 ~ 1, family = binomial(), data = tmp, weights = fint0/(pi0*pi0c))
  param1 = plogis(summary(fit1)$coef[1,1])
  #time 2
  tmp  = tmpdata[tmpdata$Y_1==1 & !is.na(tmpdata$Y_2),]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  tmp$fint1 = (tmp$L1_1*delta + 1 - tmp$L1_1)*tmp$pred_obs_1 + tmp$L1_1*tmp$A_1*(1-delta)
  if(nrow(tmp)>0) {fit2 = glm2(Y_2 ~ 1, family = binomial(), data = tmp, weights = (fint0*fint1)/(pi1*pi1c))
  param2 = plogis(summary(fit2)$coef[1,1])} else{param2 = NA}
  #time 3
  tmp  = tmpdata[tmpdata$Y_2==1 & !is.na(tmpdata$Y_3),]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  tmp$fint1 = (tmp$L1_1*delta + 1 - tmp$L1_1)*tmp$pred_obs_1 + tmp$L1_1*tmp$A_1*(1-delta)
  tmp$fint2 = (tmp$L1_2*delta + 1 - tmp$L1_2)*tmp$pred_obs_2 + tmp$L1_2*tmp$A_2*(1-delta)
  if(nrow(tmp)>0) {fit3 = glm2(Y_3 ~ 1, family = binomial(), data = tmp, weights = (fint0*fint1*fint2)/(pi2*pi2c))
  param3 = plogis(summary(fit3)$coef[1,1])} else{param3 = NA}
  #time 4
  tmp  = tmpdata[tmpdata$Y_3==1 & !is.na(tmpdata$Y_4),]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  tmp$fint1 = (tmp$L1_1*delta + 1 - tmp$L1_1)*tmp$pred_obs_1 + tmp$L1_1*tmp$A_1*(1-delta)
  tmp$fint2 = (tmp$L1_2*delta + 1 - tmp$L1_2)*tmp$pred_obs_2 + tmp$L1_2*tmp$A_2*(1-delta)
  tmp$fint3 = (tmp$L1_3*delta + 1 - tmp$L1_3)*tmp$pred_obs_3 + tmp$L1_3*tmp$A_3*(1-delta)
  if(nrow(tmp)>0) {fit4 = glm2(Y_4 ~ 1, family = binomial(), data = tmp, weights = (fint0*fint1*fint2*fint3)/(pi3*pi3c))
  param4 = plogis(summary(fit4)$coef[1,1])} else{param4 = NA}
  #time 5
  tmp  = tmpdata[tmpdata$Y_4==1 & !is.na(tmpdata$Y_5),]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  tmp$fint1 = (tmp$L1_1*delta + 1 - tmp$L1_1)*tmp$pred_obs_1 + tmp$L1_1*tmp$A_1*(1-delta)
  tmp$fint2 = (tmp$L1_2*delta + 1 - tmp$L1_2)*tmp$pred_obs_2 + tmp$L1_2*tmp$A_2*(1-delta)
  tmp$fint3 = (tmp$L1_3*delta + 1 - tmp$L1_3)*tmp$pred_obs_3 + tmp$L1_3*tmp$A_3*(1-delta)
  tmp$fint4 = (tmp$L1_4*delta + 1 - tmp$L1_4)*tmp$pred_obs_4 + tmp$L1_4*tmp$A_4*(1-delta)
  if(nrow(tmp)>0) {fit5 = glm2(Y_5 ~ 1, family = binomial(), data = tmp, weights = (fint0*fint1*fint2*fint3*fint4)/(pi4*pi4c))
  param5 = plogis(summary(fit5)$coef[1,1])} else{param5 = NA}
  
  #ind = ifelse(comp1==0 | comp2==0 | comp3==0 | comp4==0, 1, ind)
  
  t1 = param1
  t2 = param2*(t1)
  t3 = param3*(t2)
  t4 = param4*(t3)
  t5 = param5*(t4)
  
  myparam = c(t1, t2, t3, t4, t5)
  gc()
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"ipw_haz.csv")

stopCluster(cl)
