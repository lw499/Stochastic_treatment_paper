truer = c(0.5645818)
tmpdat = read.csv("dr_ss.csv");
tmpdat = tmpdat[tmpdat$X!=346,]
tmpdat$X=NULL;
estimates = tmpdat$meany12
tmpdat = read.csv("dr_ss_var.csv");
tmpdat = tmpdat[tmpdat$V1<1,]
tmpdat$X=NULL;
  standerror = tmpdat$V1; mean(standerror)
lower = estimates - 1.96*standerror
upper = estimates + 1.96*standerror

calculate95 = ifelse(truer<=upper & truer>=lower,1,0)
sum(calculate95)/1000

