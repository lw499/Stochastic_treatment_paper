#time varying case:
#L1: private insurance
#L2: screening of rectal STI (yes/no)
#L3: log transformed CD4 (maximal value of 10)

datagen <- function(i, K, sigma,
                    alpha0, alpha1, alpha2, alpha3, alpha4,
                    beta1_0, beta1_1, beta1_2, beta1_3, beta1_4,beta1_5, beta1_6, 
                    beta2_0, beta2_1, beta2_2, beta2_3, beta2_4, beta2_5, beta2_6,
                    beta3_0, beta3_1, beta3_2, beta3_3, beta3_4, beta3_5, beta3_6,
                    theta0, theta1, theta2, theta3, theta4, theta5, theta6,
                    cens0,cens1,cens2,cens3){
  id <- as.numeric(i) # Define individual id number
  
  # Define baseline time
  t0 <- 0
  
  # Generate baseline common cause of time-varying covariates
  U <- rbinom(1, 1, 0.7) 
  
  # Generate baseline data
  
  # Continuous covariate L1
  L1 <- rbinom(1, 1, plogis(beta1_0+beta1_1*U))
  #cavgL1 <- cumsum(L1)[1]/1 # Calculate cumavg(L1)
  # Binary covariate L2
  L2 <- rbinom(1, 1, plogis(beta2_0+beta2_1*U+beta2_2*L1))
  L3 <- rbinom(1, 1, plogis(beta3_0+beta3_1*U+beta3_3*L1))
  #min max (0,12)
  
  # Binary treatment exposure indicator A
  A <- rbinom(1, 1, plogis(alpha0+alpha1*L1+
                             alpha2*L2+alpha3*L3+alpha4*0))
  
  Cen <- rbinom(1, 1, plogis(cens0+cens1*L1+cens2*L2+cens3*L3))
  
  # Binary outcome indicator Y; write Y=2 if censored but this Y value won't be used in estimation
  
  Y <- ifelse(Cen==0,rbinom(1, 1, plogis(theta0+theta1*U+theta2*A+theta3*L1+theta4*L2+theta5*L3)),2)
  
  
  # Coerce NA to num for future data
  #enay <- 0
  #enay <- NA
  
  # Generate vectors to build data.frame for data at initial time point
  id_ <- c(id)
  t0_ <- c(t0)
  U_  <- c(U)
  L1_ <- c(L1)
  #cavgL1_ <- c(cavgL1)  
  
  L2_ = c(L2)
  L3_ = c(L3)
  A_ = c(A)
  Y_ = c(Y)
  Cen_=c(Cen)
  
  # If data for only one time point are to be generated, no more computations 
  # need to be performed
  
  # If data for multiple time points are to be generated and individual is alive
  # at end of inital measured interval, continue generating data
  if ((K > 1) && (Y==1)){
    # Generate data beyond baseline interval
    for (j in 2:K){
      # Define time interval for which data are generated
      t0 <- j-1
        Ustar=rbinom(1, 1, 1/(1+exp(-0.6+A-0.2*L1-0.4*L3))) 
        
        L1star <- rbinom(1, 1, plogis(beta1_0+beta1_1*Ustar+beta1_2*A_[j-1]+
                          beta1_3*L1_[j-1]+beta1_4*L2_[j-1]+beta1_5*L3_[j-1]+
                          beta1_6*t0))
        #temp_L1 <- c(L1_, L1star) # Store new and prior L1
        #cavgL1 <- cumsum(temp_L1)[j]/j # Calculate cumavg(L1)
        
        L2star <- rbinom(1, 1, plogis(beta2_0+beta2_1*Ustar+beta2_2*L1star+
                                        beta2_3*L2_[j-1]+beta2_4*L3_[j-1]+
                                        beta2_5*A_[j-1]+beta2_6*t0))
        L3star <- rbinom(1, 1,  plogis(beta3_0+beta3_1*Ustar+beta3_2*A_[j-1]+
                          beta3_3*L1star+beta3_4*L2star+beta3_5*L3_[j-1]+
                          beta3_6*t0))
        
        #temp_L2 <- c(L2_, L2star) # Store new nad prior L2
        #temp_L3 <- c(L3_, L3star) # Store new nad prior L3
      
      #Astar <- ifelse(A_[j-1]==0, rbinom(1, 1, plogis(alpha0+alpha1*L1+
      #                                                alpha2*L2star+alpha3*L3star+alpha4*t0)),1)
        Astar <- rbinom(1, 1, plogis(alpha0+alpha1*L1star+alpha2*L2star+alpha3*L3star+alpha4*A_[j-1]))
  
      #temp_A <- c(A_, Astar) # Store new and prior A
      
      Censtar <- rbinom(1, 1, plogis(cens0+cens1*L1star+cens2*L2star+cens3*L3star))
      
      Ystar <- ifelse(Censtar==0,rbinom(1, 1, plogis(theta0+theta1*Ustar+theta2*Astar+
                                                    theta3*L1star+theta4*L2star+theta5*L3star+
                                                    theta6*t0)),2)
      
      # Finalize new data to add to longitudinal data frame
      
      id_[j]        <- id
      t0_[j]        <- t0
      U_[j]         <- Ustar
      L1_[j]        <- L1star
      L2_[j]        <- L2star
      L3_[j]        <- L3star
      
      #cavgL1_[j]    <- cavgL1
      A_[j]         <- Astar
      
      Y_[j]         <- Ystar
      
      Cen_[j]       <- Censtar
      
      if(Ystar==0 | Ystar==2) break
      
    }
    
  }
  
  # Consolidate data in a single data frame
  temp_data <- data.frame(id = id_, t0 = t0_, U = U_, L1 = L1_, 
                          L2 = L2_, L3 = L3_, A = A_, 
                          Y = Y_, Cen=Cen_)
  return(temp_data)
}


