set.seed(100388)

rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
N <- 10000000
beta0=2; beta1=1; beta2=-1; beta3=0.5; beta4=1;
gamma0=1.5; gamma1=-1; gamma2=-0.5; gamma3=1; gamma4=0.25; gamma5=1;
alpha0=-3; alpha1=1; alpha2=-0.5; alpha3=0.25; alpha4=0.5; alpha5=0.25;
theta0=-1; theta1=2; theta2=-2; theta3=0.5; theta4=-0.25; theta5=0.5; theta6=0.75;
sigma=1

ua <- rep(TRUE, N)
U <- rbinom(N, 1, 0.7)
W1 <- W2 <- L10 <- L20 <- A0 <- Y1 <- L11 <- L21 <- A1 <- Y2 <- L12 <- L22 <- A2 <- Y3 <- L13 <- L23 <- A3 <- Y4 <- L14 <- L24 <- A4 <- Q1 <-Q2 <- Q3 <- as.numeric(rep(NA, N))
ids <- as.list(1:N)
W1 <- rbinom(N, 1, 0.5) #sex
W2 <- rnorm(N, 0, sigma)  #transformed age
L20 = rnorm(N, mean=(beta0+beta4*W1+0*U), sd=sigma) #CD4
L10 = as.numeric(rexpit(gamma0+gamma2*L20+gamma3*W1+gamma4*W2+0*U) ) #AIDS -> Lstar
A0 <- rexpit(alpha0+alpha1*L10+alpha2*L20+alpha3*L20*L10+alpha4*W1+alpha5*W2+0.5*abs(W2) )
Y1 <- plogis(theta0+theta1*A0+theta2*L10+theta3*0+theta4*L10*L20+theta5*W1+theta6*(abs(L20+W2)^1.5)+0*U)

L21 = rnorm(N, mean=(beta0+beta1*A0+beta2*L10+beta3*L20+beta4*W1+0*U), sd=sigma) #CD4
L11 = as.numeric(rexpit(gamma0+gamma1*A0+gamma2*L21+gamma3*W1+gamma4*W2+gamma5*L10+0*U) ) #AIDS
A1 = as.numeric(rexpit(alpha0+alpha1*L11+alpha2*L21+alpha3*L21*L11+alpha4*W1+alpha5*W2+0.5*abs(W2) ) | A0)
Y2 = plogis(theta0+theta1*A1+theta2*L11+theta3*0+theta4*L11*L21+theta5*W1+theta6*(abs(L21+W2)^1.5)+0*U)

L22 = rnorm(N, mean=(beta0+beta1*A1+beta2*L11+beta3*L21+beta4*W1+0*U), sd=sigma) #CD4
L12 = as.numeric(rexpit(gamma0+gamma1*A1+gamma2*L22+gamma3*W1+gamma4*W2+gamma5*L11+0*U) ) 
A2 = as.numeric(rexpit(alpha0+alpha1*L12+alpha2*L22+alpha3*L22*L12+alpha4*W1+alpha5*W2+0.5*abs(W2) ) | A1)
Y3 = plogis(theta0+theta1*A2+theta2*L12+theta3*0+theta4*L12*L22+theta5*W1+theta6*(abs(L22+W2)^1.5)+0*U)

L23 = rnorm(N, mean=(beta0+beta1*A2+beta2*L12+beta3*L22+beta4*W1+0*U), sd=sigma) #CD4
L13 = as.numeric(rexpit(gamma0+gamma1*A2+gamma2*L23+gamma3*W1+gamma4*W2+gamma5*L12+0*U) ) 
A3 = as.numeric(rexpit(alpha0+alpha1*L13+alpha2*L23+alpha3*L23*L13+alpha4*W1+alpha5*W2+0.5*abs(W2) ) | A2)
Y4 = plogis(theta0+theta1*A3+theta2*L13+theta3*0+theta4*L13*L23+theta5*W1+theta6*(abs(L23+W2)^1.5)+0*U)

L24 = rnorm(N, mean=(beta0+beta1*A3+beta2*L13+beta3*L23+beta4*W1+0*U), sd=sigma) #CD4
L14 = as.numeric(rexpit(gamma0+gamma1*A3+gamma2*L24+gamma3*W1+gamma4*W2+gamma5*L13+0*U) ) 
A4 = as.numeric(rexpit(alpha0+alpha1*L14+alpha2*L24+alpha3*L24*L14+alpha4*W1+alpha5*W2+0.5*abs(W2) ) | A3)
Y5 = plogis(theta0+theta1*A4+theta2*L14+theta3*0+theta4*L14*L24+theta5*W1+theta6*(abs(L24+W2)^1.5)+0*U)

surv1  = mean(Y1)
surv2  = mean((Y1)*(Y2))
surv3  = mean((Y1)*(Y2)*(Y3))
surv4  = mean((Y1)*(Y2)*(Y3)*(Y4))
surv5  = mean((Y1)*(Y2)*(Y3)*(Y4)*(Y5))
c(surv1, surv2, surv3, surv4, surv5)
#0.6016705 0.4951122 0.4428695 0.4092573 0.3850652

############################
delta = 0.5
set.seed(100388)

EXPITT <- function(term) {
    return( ifelse(!is.na(term),1/(1+exp(term)),NA) )
}

ua <- rep(TRUE, N)
U <- rbinom(N, 1, 0.7)
W1 <- W2 <- L10 <- L20 <- A0 <- Y1 <- L11 <- L21 <- A1 <- Y2 <- L12 <- L22 <- A2 <- Y3 <- L13 <- L23 <- A3 <- Y4 <- L14 <- L24 <- A4 <- Q1 <-Q2 <- Q3 <- as.numeric(rep(NA, N))
ids <- as.list(1:N)
W1 <- rbinom(N, 1, 0.5) #sex
W2 <- rnorm(N, 0, sigma)  #transformed age
L20 = rnorm(N, mean=(beta0+beta4*W1+0*U), sd=sigma) #CD4
L10 = as.numeric(rexpit(gamma0+gamma2*L20+gamma3*W1+gamma4*W2+0*U) ) #AIDS
A0 <- rexpit(alpha0+alpha1*L10+alpha2*L20+alpha3*L20*L10+alpha4*W1+alpha5*W2+0.5*abs(W2) )
    Q0 = I(L10==1)
    tmp = rbinom(N,1,1-delta*EXPITT(alpha0+alpha1*L10+alpha2*L20+alpha3*L20*L10+alpha4*W1+alpha5*W2+0.5*abs(W2) ))
    A0 = ifelse(Q0==1,tmp,A0)
Y1 <- plogis(theta0+theta1*A0+theta2*L10+theta3*0+theta4*L10*L20+theta5*W1+theta6*(abs(L20+W2)^1.5)+0*U)

L21 = rnorm(N, mean=(beta0+beta1*A0+beta2*L10+beta3*L20+beta4*W1+0*U), sd=sigma) #CD4
L11 = as.numeric(rexpit(gamma0+gamma1*A0+gamma2*L21+gamma3*W1+gamma4*W2+gamma5*L10+0*U) ) #AIDS
A1 = as.numeric(rexpit(alpha0+alpha1*L11+alpha2*L21+alpha3*L21*L11+alpha4*W1+alpha5*W2+0.5*abs(W2) ) | A0)
    Q1 = I(L11==1)
    tmp = rbinom(N,1, 1-delta*EXPITT(alpha0+alpha1*L11+alpha2*L21+alpha3*L21*L11+alpha4*W1+alpha5*W2+0.5*abs(W2) ))
    A1 = ifelse(Q1==1 & A0==0,tmp,A1)
Y2 = plogis(theta0+theta1*A1+theta2*L11+theta3*0+theta4*L11*L21+theta5*W1+theta6*(abs(L21+W2)^1.5)+0*U)

L22 = rnorm(N, mean=(beta0+beta1*A1+beta2*L11+beta3*L21+beta4*W1+0*U), sd=sigma) #CD4
L12 = as.numeric(rexpit(gamma0+gamma1*A1+gamma2*L22+gamma3*W1+gamma4*W2+gamma5*L11+0*U) ) 
A2 = as.numeric(rexpit(alpha0+alpha1*L12+alpha2*L22+alpha3*L22*L12+alpha4*W1+alpha5*W2+0.5*abs(W2) ) | A1)
    Q2 = I(L12==1)
    tmp = rbinom(N,1, 1-delta*EXPITT(alpha0+alpha1*L12+alpha2*L22+alpha3*L22*L12+alpha4*W1+alpha5*W2+0.5*abs(W2) ))
    A2 = ifelse(Q2==1 & A1==0,tmp,A2)
Y3 = plogis(theta0+theta1*A2+theta2*L12+theta3*0+theta4*L12*L22+theta5*W1+theta6*(abs(L22+W2)^1.5)+0*U)

L23 = rnorm(N, mean=(beta0+beta1*A2+beta2*L12+beta3*L22+beta4*W1+0*U), sd=sigma) #CD4
L13 = as.numeric(rexpit(gamma0+gamma1*A2+gamma2*L23+gamma3*W1+gamma4*W2+gamma5*L12+0*U) ) 
A3 = as.numeric(rexpit(alpha0+alpha1*L13+alpha2*L23+alpha3*L23*L13+alpha4*W1+alpha5*W2+0.5*abs(W2) ) | A2)
    Q3 = I(L13==1)
    tmp = rbinom(N,1, 1-delta*EXPITT(alpha0+alpha1*L13+alpha2*L23+alpha3*L23*L13+alpha4*W1+alpha5*W2+0.5*abs(W2) ))
    A3 = ifelse(Q3==1 & A2==0,tmp,A3)
Y4 = plogis(theta0+theta1*A3+theta2*L13+theta3*0+theta4*L13*L23+theta5*W1+theta6*(abs(L23+W2)^1.5)+0*U)

L24 = rnorm(N, mean=(beta0+beta1*A3+beta2*L13+beta3*L23+beta4*W1+0*U), sd=sigma) #CD4
L14 = as.numeric(rexpit(gamma0+gamma1*A3+gamma2*L24+gamma3*W1+gamma4*W2+gamma5*L13+0*U) ) 
A4 = as.numeric(rexpit(alpha0+alpha1*L14+alpha2*L24+alpha3*L24*L14+alpha4*W1+alpha5*W2+0.5*abs(W2) ) | A3)
    Q4 = I(L14==1)
    tmp = rbinom(N,1, 1-delta*EXPITT(alpha0+alpha1*L14+alpha2*L24+alpha3*L24*L14+alpha4*W1+alpha5*W2+0.5*abs(W2) ))
    A4 = ifelse(Q4==1 & A3==0,tmp,A4)
Y5 = plogis(theta0+theta1*A4+theta2*L14+theta3*0+theta4*L14*L24+theta5*W1+theta6*(abs(L24+W2)^1.5)+0*U)

surv1  = mean(Y1)
surv2  = mean((Y1)*(Y2))
surv3  = mean((Y1)*(Y2)*(Y3))
surv4  = mean((Y1)*(Y2)*(Y3)*(Y4))
surv5  = mean((Y1)*(Y2)*(Y3)*(Y4)*(Y5))
c(surv1, surv2, surv3, surv4, surv5)
#0.6807162 0.6117451 0.5858175 0.5725453 0.5645818





