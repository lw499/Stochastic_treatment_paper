library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(7)                                              
registerDoParallel(cl)                                                                          
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  options(warn=-1)
  #library(geepack);
  library(MASS);library(ResourceSelection);
  #library(ltmle); library(SuperLearner)
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
  
  n <- 500
  K <- 5
  delta=0.25
  alpha0=-1; alpha1=-2; alpha2=-1; alpha3=1; alpha4=2; 
  beta1_0=-1; beta1_1=0; beta1_2=-1; beta1_3=1; beta1_4=-1; beta1_5=1; beta1_6=0; #rectal STI indicator
  beta2_0=-1; beta2_1=0; beta2_2=0; beta2_3=1; beta2_4=-1; beta2_5=1; beta2_6=0; #cd4 count
  beta3_0=1; beta3_1=0; beta3_2=1; beta3_3=1; beta3_4=0; beta3_5=1; beta3_6=0;  #unprotected sexual activity (H/L)
  theta0=1; theta1=0; theta2=3; theta3=-2; theta4=1; theta5=-1; theta6=0; #theta1 always 0
  cens0=-2; cens1=0; cens2=1; cens3=-1
  
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen(ind, K=K, sigma=1,
            alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3, alpha4 = alpha4,
            beta1_0=beta1_0, beta1_1=beta1_1, beta1_2=beta1_2, beta1_3=beta1_3, beta1_4=beta1_4, beta1_5=beta1_5, beta1_6=beta1_6,
            beta2_0=beta2_0, beta2_1=beta2_1, beta2_2=beta2_2, beta2_3=beta2_3, beta2_4=beta2_4, beta2_5=beta2_5, beta2_6=beta2_6, 
            beta3_0=beta3_0, beta3_1=beta3_1, beta3_2=beta3_2, beta3_3=beta3_3, beta3_4=beta3_4, beta3_5=beta3_5, beta3_6=beta3_6, 
            theta0=theta0, theta1=theta1, theta2=theta2, theta3=theta3, theta4=theta4, theta5=theta5, theta6=theta6,
            cens0=cens0, cens1=cens1, cens2=cens2, cens3=cens3)
  })
  
  dffull <- rbindlist(df)

  dffull[, paste("lag_A") := shift(A, 1, NA, type='lag'), by=id]
  dffull$lag_A = ifelse(dffull$t0==0, 0, dffull$lag_A)
  
  afit = glm2(A ~ L1 + L2 + L3, family = binomial(), data = dffull) #This is from the data generation mechanism
  dffull$pred_obs = predict(afit, newdata = dffull, type="response")
  dffull$pred_obs = ifelse(dffull$A==1, dffull$pred_obs, 1-dffull$pred_obs)
  dffull$lag_A=NULL
  
  #cfit =  glm2(Cen~L2+L3,family=binomial,data=dffull)
  dffull$pred_obsc = 1#predict(cfit, newdata = dffull, type="response")
  #dffull$pred_obsc = ifelse(dffull$Cen==1, dffull$pred_obsc, 1-dffull$pred_obsc)
  
  dffullwide = dcast(dffull, id ~ t0, value.var = c("L1","L2","L3","A","Cen","Y","pred_obs","pred_obsc","U"))
  
  tmpdata = dffullwide
  #subset data
  tmpdata$Y_1 = ifelse(tmpdata$Y_0==0,0,tmpdata$Y_1)
  tmpdata$Y_2 = ifelse(!is.na(tmpdata$Y_1) & tmpdata$Y_1==0,0,tmpdata$Y_2)
  tmpdata$Y_3 = ifelse(!is.na(tmpdata$Y_2) & tmpdata$Y_2==0,0,tmpdata$Y_3)
  tmpdata$Y_4 = ifelse(!is.na(tmpdata$Y_3) & tmpdata$Y_3==0,0,tmpdata$Y_4)
  
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
  tmp = tmpdata[tmpdata$Y_1<2,]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  fit1 = glm2(Y_1 ~ 1, family = binomial(), data = tmp, weights = fint0/(pi0*pi0c))
  param1 = plogis(summary(fit1)$coef[1,1])
  #time 2
  tmp  = tmpdata[tmpdata$Y_1==1 & tmpdata$Y_2<2,]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  tmp$fint1 = (tmp$L1_1*delta + 1 - tmp$L1_1)*tmp$pred_obs_1 + tmp$L1_1*tmp$A_1*(1-delta)
  if(nrow(tmp)>0) {fit2 = glm2(Y_2 ~ 1, family = binomial(), data = tmp, weights = (fint0*fint1)/(pi1*pi1c))
  param2 = plogis(summary(fit2)$coef[1,1])} else{param2 = NA}
  #time 3
  tmp  = tmpdata[tmpdata$Y_2==1 & tmpdata$Y_3<2,]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  tmp$fint1 = (tmp$L1_1*delta + 1 - tmp$L1_1)*tmp$pred_obs_1 + tmp$L1_1*tmp$A_1*(1-delta)
  tmp$fint2 = (tmp$L1_2*delta + 1 - tmp$L1_2)*tmp$pred_obs_2 + tmp$L1_2*tmp$A_2*(1-delta)
  if(nrow(tmp)>0) {fit3 = glm2(Y_3 ~ 1, family = binomial(), data = tmp, weights = (fint0*fint1*fint2)/(pi2*pi2c))
  param3 = plogis(summary(fit3)$coef[1,1])} else{param3 = NA}
  #time 4
  tmp  = tmpdata[tmpdata$Y_3==1 & tmpdata$Y_4<2,]
  tmp$fint0 = (tmp$L1_0*delta + 1 - tmp$L1_0)*tmp$pi0 + tmp$L1_0*tmp$A_0*(1-delta)
  tmp$fint1 = (tmp$L1_1*delta + 1 - tmp$L1_1)*tmp$pred_obs_1 + tmp$L1_1*tmp$A_1*(1-delta)
  tmp$fint2 = (tmp$L1_2*delta + 1 - tmp$L1_2)*tmp$pred_obs_2 + tmp$L1_2*tmp$A_2*(1-delta)
  tmp$fint3 = (tmp$L1_3*delta + 1 - tmp$L1_3)*tmp$pred_obs_3 + tmp$L1_3*tmp$A_3*(1-delta)
  if(nrow(tmp)>0) {fit4 = glm2(Y_4 ~ 1, family = binomial(), data = tmp, weights = (fint0*fint1*fint2*fint3)/(pi3*pi3c))
  param4 = plogis(summary(fit4)$coef[1,1])} else{param4 = NA}
  #time 5
  tmp  = tmpdata[tmpdata$Y_4==1 & tmpdata$Y_5<2,]
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
  
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"ipw_haz_wrong.csv")

stopCluster(cl)
