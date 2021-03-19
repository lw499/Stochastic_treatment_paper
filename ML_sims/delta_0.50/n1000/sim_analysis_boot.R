truer = c(0.5645818)
tmpdat = read.csv("dr_ss.csv");
tmpdat$X=NULL;
estimates = tmpdat$meany12
tmpdat = read.csv("dr_ss_var.csv");
tmpdat$X=NULL;
  standerror = tmpdat$V1; mean(standerror)
lower = estimates - 1.96*standerror
upper = estimates + 1.96*standerror

calculate95 = ifelse(truer<=upper & truer>=lower,1,0)
sum(calculate95)/1000

