set.seed(12345)

##natural value
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
N <- 1000000
Ndat <- 5
delta=0.75
alpha0=-1; alpha1=-2; alpha2=-1; alpha3=1; alpha4=2; 
beta1_0=-1; beta1_1=0; beta1_2=-1; beta1_3=1; beta1_4=-1; beta1_5=1; beta1_6=0; #rectal STI indicator
beta2_0=-1; beta2_1=0; beta2_2=0; beta2_3=1; beta2_4=-1; beta2_5=1; beta2_6=0; #cd4 count
beta3_0=1; beta3_1=0; beta3_2=1; beta3_3=1; beta3_4=0; beta3_5=1; beta3_6=0;  #unprotected sexual activity (H/L)
theta0=1; theta1=0; theta2=3; theta3=-2; theta4=1; theta5=-1; theta6=0; #theta1 always 0
#cens0=-2; cens1=0; cens2=1; cens3=-0.1
sigma=1
#all time points: probability highest in ppl w/ no insurance, with rectal STI is less than 30%

#ua <- rep(TRUE, N)
ids <- as.list(1:N)
L10 <- rbinom(N, 1, plogis(-5))
A0 = rbinom(N, 1, plogis(alpha0+alpha1*L10))
Y1 <- plogis(-2+theta2*A0)
mean(Y1)

L10 <- rep(1,N)
A0 = rbinom(N, 1, plogis(alpha0+alpha1*L10))
Y1 <- plogis(-2+theta2*A0)
mean(Y1)

A0 = rbinom(N, 1, plogis(-3))
Y1 <- plogis(-2+theta2*A0)
mean(Y1)


