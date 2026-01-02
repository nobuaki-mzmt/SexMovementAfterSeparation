## 170829 analysis



d <- read.delim("clipboard",header=T)

step <- 180
Xrange <- c(0,step)
Yrange <- c(0,0.35)
par(mfrow=c(1,1), pin=c(4,3)*1.2)
matplot(1:(step*10)/10, t(d[,1:(step*10)+2])/100000, type="l", lty=1, ylim=Yrange, xlim=Xrange,
        xlab="time (sec)", ylab="encounter rate", las=1)
legend(100,0.15,c("b-b", "b-a", "a-b", "a-a"), col=1:4, lty=1)




d <- read.delim("clipboard",header=T)

step <- 180
Xrange <- c(0,step)
Yrange <- c(0,0.35)
par(mfrow=c(1,1), pin=c(4,3)*1.2)
matplot(1:(step*10)/10, t(d[,1:(step*10)+2])/100000, type="l", lty=1, ylim=Yrange, xlim=Xrange,
        xlab="time (sec)", ylab="encounter rate", las=1)
legend(100,0.15,c("b-b", "b-a", "a-b", "a-a"), col=1:4, lty=1)



setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Simulations")
d <- read.csv("Rs_CRWmodel_NBC_rep100000_d0-100.csv",header=T)
step <- 1800
Xrange <- c(0,step)
Yrange <- c(0,0.08)
par(mfrow=c(1,1), pin=c(4,3)*1.2)
matplot(1:(step*10)/10, t(d[,1:(step*10)+2])/100000, type="l", lty=1, ylim=Yrange, xlim=Xrange,
        xlab="time (sec)", ylab="encounter rate", las=1)
legend(100,0.15,c("b-b", "b-a", "a-b", "a-a"), col=1:4, lty=1)





