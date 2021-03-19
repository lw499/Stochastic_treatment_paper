set.seed(12345)

rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
N <- 10000000
Ndat <- 5
delta=0.25
alpha0=-1; alpha1=-2; alpha2=-1; alpha3=1; alpha4=2; 
beta1_0=-1; beta1_1=0; beta1_2=-1; beta1_3=1; beta1_4=-1; beta1_5=1; beta1_6=0; #rectal STI indicator
beta2_0=-1; beta2_1=0; beta2_2=0; beta2_3=1; beta2_4=-1; beta2_5=1; beta2_6=0; #cd4 count
beta3_0=1; beta3_1=0; beta3_2=1; beta3_3=1; beta3_4=0; beta3_5=1; beta3_6=0;  #unprotected sexual activity (H/L)
theta0=1; theta1=0; theta2=3; theta3=-2; theta4=1; theta5=-1; theta6=0; #theta1 always 0
#cens0=-2; cens1=0; cens2=1; cens3=-0.1
sigma=1
#all time points: probability highest in ppl w/ no insurance, with rectal STI is less than 30%

#ua <- rep(TRUE, N)
L10 <- L20 <- A0 <- Y1 <- L11 <- L21 <- A1 <- Y2 <- L12 <- L22 <- A2 <- Y3 <- L13 <- L23 <- A3 <- Y4 <- L14 <- L24 <- A4 <- Y5 <- L15 <- L25 <- A5 <- Y6 <- L16 <- L26 <- A6 <- Y7 <- as.numeric(rep(NA, N))
L30 <- L31 <- L32 <- L33 <- L34 <- L35 <- as.numeric(rep(NA, N))

ids <- as.list(1:N)
U0 <- rbinom(N, 1, 0.7)
L10 <- rbinom(N, 1, plogis(beta1_0+beta1_1*U0))
L20 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U0+beta2_2*L10))
L30 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U0+beta3_3*L10))
A0 = rbinom(N, 1, plogis(alpha0+alpha1*L10+alpha2*L20+alpha3*L30+alpha4*0))
Y1 <- plogis(theta0+theta1*U0+theta2*A0+theta3*L10+theta4*L20+theta5*L30)

U1=rbinom(N, 1, 1/(1+exp(-0.6+A0-0.2*L10-0.4*L30))) 
L11 = rexpit(beta1_0+beta1_1*U1+beta1_2*A0+beta1_3*L10+beta1_4*L20+beta1_5*L30+beta1_6*1)
L21 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U1+beta2_2*L11+beta2_3*L20+beta2_4*L30+beta2_5*A0+beta2_6*1))
L31 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U1+beta3_2*A0+beta3_3*L11+beta3_4*L21+beta3_5*L30+beta3_6*1))
A1 = rbinom(N, 1, plogis(alpha0+alpha1*L11+alpha2*L21+alpha3*L31+alpha4*A0)) 
Y2 = plogis((theta0+theta1*U1+theta2*A1+theta3*L11+theta4*L21+theta5*L31+theta6*1))

U2=rbinom(N, 1, 1/(1+exp(-0.6+A1-0.2*L11-0.4*L31))) 
L12 = rexpit(beta1_0+beta1_1*U2+beta1_2*A1+beta1_3*L11+beta1_4*L21+beta1_5*L31+beta1_6*2)
L22 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U2+beta2_2*L12+beta2_3*L21+beta2_4*L31+beta2_5*A1+beta2_6*2))
L32 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U2+beta3_2*A1+beta3_3*L12+beta3_4*L22+beta3_5*L31+beta3_6*2))
A2 = rbinom(N, 1, plogis(alpha0+alpha1*L12+alpha2*L22+alpha3*L32+alpha4*A1)) 
Y3 = plogis((theta0+theta1*U2+theta2*A2+theta3*L12+theta4*L22+theta5*L32+theta6*2))

U3=rbinom(N, 1, 1/(1+exp(-0.6+A2-0.2*L12-0.4*L32))) 
L13 = rexpit(beta1_0+beta1_1*U3+beta1_2*A2+beta1_3*L12+beta1_4*L22+beta1_5*L32+beta1_6*3)
L23 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U3+beta2_2*L13+beta2_3*L22+beta2_4*L32+beta2_5*A2+beta2_6*3))
L33 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U3+beta3_2*A2+beta3_3*L13+beta3_4*L23+beta3_5*L32+beta3_6*3))
A3 = rbinom(N, 1, plogis(alpha0+alpha1*L13+alpha2*L23+alpha3*L33+alpha4*A2)) 
Y4 = plogis((theta0+theta1*U3+theta2*A3+theta3*L13+theta4*L23+theta5*L33+theta6*3))

U4=rbinom(N, 1, 1/(1+exp(-0.6+A3-0.2*L13-0.4*L33))) 
L14 = rexpit(beta1_0+beta1_1*U4+beta1_2*A3+beta1_3*L13+beta1_4*L23+beta1_5*L33+beta1_6*4)
L24 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U4+beta2_2*L14+beta2_3*L23+beta2_4*L33+beta2_5*A3+beta2_6*4))
L34 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U4+beta3_2*A3+beta3_3*L14+beta3_4*L24+beta3_5*L33+beta3_6*4))
A4 = rbinom(N, 1, plogis(alpha0+alpha1*L14+alpha2*L24+alpha3*L34+alpha4*A3)) 
Y5 = plogis((theta0+theta1*U4+theta2*A4+theta3*L14+theta4*L24+theta5*L34+theta6*4))

surv1 = mean(Y1)
surv2 = mean((Y1)*(Y2))
surv3 = mean((Y1)*(Y2)*(Y3))
surv4 = mean((Y1)*(Y2)*(Y3)*(Y4))
surv5 = mean((Y1)*(Y2)*(Y3)*(Y4)*(Y5))

c(surv1, surv2, surv3, surv4, surv5)
#0.6242012 0.4190589 0.2985796 0.2182048 0.1608798

#### Incremental
L10 <- L20 <- A0 <- Y1 <- L11 <- L21 <- A1 <- Y2 <- L12 <- L22 <- A2 <- Y3 <- L13 <- L23 <- A3 <- Y4 <- L14 <- L24 <- A4 <- Y5 <- L15 <- L25 <- A5 <- Y6 <- L16 <- L26 <- A6 <- Y7 <- as.numeric(rep(NA, N))
L30 <- L31 <- L32 <- L33 <- L34 <- L35 <- as.numeric(rep(NA, N))

EXPITT <- function(term) {
  return( ifelse(!is.na(term),1/(1+exp(term)),NA) )
}

ids <- as.list(1:N)
U0 <- rbinom(N, 1, 0.7)
L10 <- rbinom(N, 1, plogis(beta1_0+beta1_1*U0))
L20 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U0+beta2_2*L10))
L30 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U0+beta3_3*L10))
A0 = rbinom(N, 1, plogis(alpha0+alpha1*L10+alpha2*L20+alpha3*L30+alpha4*0))
  Q0 = I(L10==1)
  #tmp = rbinom(N,1, (1-delta) + delta*plogis(alpha0+alpha1*L10+alpha2*L20+alpha3*L30+alpha4*0))
  tmp = rbinom(N,1,1-delta*EXPITT(alpha0+alpha1*L10+alpha2*L20+alpha3*L30+alpha4*0))
  A0 = ifelse(Q0==1,tmp,A0)
Y1 <- plogis(theta0+theta1*U0+theta2*A0+theta3*L10+theta4*L20+theta5*L30)

U1=rbinom(N, 1, 1/(1+exp(-0.6+A0-0.2*L10-0.4*L30))) 
L11 = rexpit(beta1_0+beta1_1*U1+beta1_2*A0+beta1_3*L10+beta1_4*L20+beta1_5*L30+beta1_6*1)
L21 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U1+beta2_2*L11+beta2_3*L20+beta2_4*L30+beta2_5*A0+beta2_6*1))
L31 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U1+beta3_2*A0+beta3_3*L11+beta3_4*L21+beta3_5*L30+beta3_6*1))
A1 = rbinom(N, 1, plogis(alpha0+alpha1*L11+alpha2*L21+alpha3*L31+alpha4*A0)) 
  Q1 = I(L11==1)
  #tmp = rbinom(N,1, (1-delta) + delta*plogis(alpha0+alpha1*L11+alpha2*L21+alpha3*L31+alpha4*A0))
  tmp = rbinom(N,1, 1-delta*EXPITT(alpha0+alpha1*L11+alpha2*L21+alpha3*L31+alpha4*A0))
  A1 = ifelse(Q1==1,tmp,A1)
Y2 = plogis((theta0+theta1*U1+theta2*A1+theta3*L11+theta4*L21+theta5*L31+theta6*1))

U2=rbinom(N, 1, 1/(1+exp(-0.6+A1-0.2*L11-0.4*L31))) 
L12 = rexpit(beta1_0+beta1_1*U2+beta1_2*A1+beta1_3*L11+beta1_4*L21+beta1_5*L31+beta1_6*2)
L22 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U2+beta2_2*L12+beta2_3*L21+beta2_4*L31+beta2_5*A1+beta2_6*2))
L32 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U2+beta3_2*A1+beta3_3*L12+beta3_4*L22+beta3_5*L31+beta3_6*2))
A2 = rbinom(N, 1, plogis(alpha0+alpha1*L12+alpha2*L22+alpha3*L32+alpha4*A1)) 
  Q2 = I(L12==1)
  #tmp = rbinom(N,1, (1-delta) + delta*plogis(alpha0+alpha1*L12+alpha2*L22+alpha3*L32+alpha4*A1))
  tmp = rbinom(N,1, 1-delta*EXPITT(alpha0+alpha1*L12+alpha2*L22+alpha3*L32+alpha4*A1))
  A2 = ifelse(Q2==1,tmp,A2)
Y3 = plogis((theta0+theta1*U2+theta2*A2+theta3*L12+theta4*L22+theta5*L32+theta6*2))

U3=rbinom(N, 1, 1/(1+exp(-0.6+A2-0.2*L12-0.4*L32))) 
L13 = rexpit(beta1_0+beta1_1*U3+beta1_2*A2+beta1_3*L12+beta1_4*L22+beta1_5*L32+beta1_6*3)
L23 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U3+beta2_2*L13+beta2_3*L22+beta2_4*L32+beta2_5*A2+beta2_6*3))
L33 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U3+beta3_2*A2+beta3_3*L13+beta3_4*L23+beta3_5*L32+beta3_6*3))
A3 = rbinom(N, 1, plogis(alpha0+alpha1*L13+alpha2*L23+alpha3*L33+alpha4*A2)) 
  Q3 = I(L13==1)
  #tmp = rbinom(N,1, (1-delta) + delta*plogis(alpha0+alpha1*L13+alpha2*L23+alpha3*L33+alpha4*A2))
  tmp = rbinom(N,1, 1-delta*EXPITT(alpha0+alpha1*L13+alpha2*L23+alpha3*L33+alpha4*A2))
  A3 = ifelse(Q3==1,tmp,A3)
Y4 = plogis((theta0+theta1*U3+theta2*A3+theta3*L13+theta4*L23+theta5*L33+theta6*3))

#beta3_0=6.5; beta3_1=-0.5; beta3_2=0.5; beta3_3=0; beta3_4=0; beta3_5=0.1; beta3_6=0; 
#alpha0=1; alpha1=1; alpha2=-3; alpha3=-0.1; alpha4=0; #probability highest in ppl w/ no insurance, with rectal STI is less than 15%
U4=rbinom(N, 1, 1/(1+exp(-0.6+A3-0.2*L13-0.4*L33))) 
L14 = rexpit(beta1_0+beta1_1*U4+beta1_2*A3+beta1_3*L13+beta1_4*L23+beta1_5*L33+beta1_6*4)
L24 <- rbinom(N, 1, plogis(beta2_0+beta2_1*U4+beta2_2*L14+beta2_3*L23+beta2_4*L33+beta2_5*A3+beta2_6*4))
L34 <- rbinom(N, 1, plogis(beta3_0+beta3_1*U4+beta3_2*A3+beta3_3*L14+beta3_4*L24+beta3_5*L33+beta3_6*4))
A4 = rbinom(N, 1, plogis(alpha0+alpha1*L14+alpha2*L24+alpha3*L34+alpha4*A3)) 
  Q4 = I(L14==1)
  tmp = rbinom(N,1, 1-delta*EXPITT(alpha0+alpha1*L14+alpha2*L24+alpha3*L34+alpha4*A3))
  A4 = ifelse(Q4==1,tmp,A4)
Y5 = plogis((theta0+theta1*U4+theta2*A4+theta3*L14+theta4*L24+theta5*L34+theta6*4))

surv1 = mean(Y1)
surv2 = mean((Y1)*(Y2))
surv3 = mean((Y1)*(Y2)*(Y3))
surv4 = mean((Y1)*(Y2)*(Y3)*(Y4))
surv5 = mean((Y1)*(Y2)*(Y3)*(Y4)*(Y5))

c(surv1, surv2, surv3, surv4, surv5)
#0.7351876 0.5823728 0.4764706 0.3937580 0.3263147

