truer = c(0.3263147)

tmpdat = read.csv("seq_reg.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean
gformseq2 = tmpdat;
gformseqSE = apply(gformseq2, 2, var)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = 100*(gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, sqrt(MSE)*100)
mean(abs(tmpdat$meany5 - truer))*100

tmpdat = read.csv("seq_reg_non_sat.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean
gformseq2 = tmpdat;
gformseqSE = apply(gformseq2, 2, var)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = 100*(gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, sqrt(MSE)*100)
mean(abs(tmpdat$meany5 - truer))*100

tmpdat = read.csv("ipw_haz.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean[5]
gformseq2 = tmpdat$V5;
gformseqSE = var(gformseq2)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = 100*(gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, sqrt(MSE)*100)
mean(abs(tmpdat$V5 - truer))*100


tmpdat = read.csv("dr_sat.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean
gformseq2 = tmpdat;
gformseqSE = apply(gformseq2, 2, var)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = 100*(gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, sqrt(MSE)*100)
  mean(abs(tmpdat$meany5 - truer))*100


tmpdat = read.csv("dr_sat_wrongOR.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean
gformseq2 = tmpdat;
gformseqSE = apply(gformseq2, 2, var)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = 100*(gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, sqrt(MSE)*100)
mean(abs(tmpdat$meany5 - truer))*100


tmpdat = read.csv("dr_sat_wrongPS.csv")
tmpdat$X=NULL;
mean = colMeans(tmpdat)
gformseq = mean
gformseq2 = tmpdat;
gformseqSE = apply(gformseq2, 2, var)^0.5
gformseqVAR = gformseqSE^2

bpgformseq = 100*(gformseq-truer)/gformseqSE
bgformseq = (gformseq-truer)
MSE = bgformseq^2 + gformseqVAR
cbind(bgformseq*100, gformseqSE*100, sqrt(MSE)*100)
mean(abs(tmpdat$meany5 - truer))*100




