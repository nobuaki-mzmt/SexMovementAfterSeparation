## Results plot for SepSearch Sim
## 180926 N.Mizumoto
library(data.table)

## main
data.place <- "E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\SciAdv\\2ndRevision\\simulations"
time <- seq(0.2, 180, 0.2)
rep <- 1000000


## fit data
## simulations for Rs
## Random search condition
d <- read.csv(paste(data.place,"\\Rs_EnSim_rep1000000_RanSearch_L223_sec180.csv", sep=""), header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(0.2,0.4, names, col=1:4, lty=1:4)

## Reunion search 
d <- data.frame(as.matrix(fread(paste(data.place,"\\Rs_EnSim_rep1000000_SepSearch_SepDis16_sec180.csv", sep=""), header=T)))
plot_data <- matrix(as.numeric(t(d[,3:902])), ncol=4)
  
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(0.2,0.5, names, col=1:4, lty=1:4)

min(seq(0.2,200,0.2)[apply(plot_data[,1]>plot_data[,2:4], 1, sum)==3])

dad <- data.frame(as.matrix(fread(paste(
  data.place,"\\Rs_EnSim_rep10000_SepSearch_SepDis16_sec300000.csv", sep=""), header=F)))

x <- c(10:15000)
plot_data2 <- matrix(as.numeric(t(dad[,x+3])), ncol=4)
max(x[apply(plot_data2[,1]>plot_data2[,2:4], 1, sum)==3])*100*0.2
# 241840
matplot(c(time,x*100*0.2), rbind(plot_data/1000000,plot_data2/10000), type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")

matplot(x*100*0.2, plot_data2/10000, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")


## simulations for Cf
## Random search condition
d <- read.csv(paste(data.place,"\\Cf_EnSim_rep1000000_RanSearch_L223_sec180.csv", sep=""), header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(0.2,0.8, names, col=1:4, lty=1:4)

## Reunion search 
d <- read.csv(paste(data.place,"\\Cf_EnSim_rep1000000_SepSearch_SepDis22_sec200.csv", sep=""), header=T)
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
time2 <- seq(0.2,200,0.2)
matplot(time2, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(5,0.15, names, col=1:4, lty=1:4)
min(seq(0.2,200,0.2)[apply(plot_data[,1]>plot_data[,2:4], 1, sum)==3])
max(seq(0.2,200,0.2)[apply(plot_data[,1]>plot_data[,2:4], 1, sum)==3])






####################################
## for ad sim
d <- data.frame(as.matrix(fread(paste(data.place,"\\Add_EnSim_rep1000000_SepSearch_SepDis16_sec180.csv", sep=""), header=T)))
time <- seq(0.2, 180, 0.2)
names <- c("sim", "Pause", "wo pause")
plot_data <- matrix(as.numeric(t(d[,3:length(d[1,])])), ncol=3)
plot_data <- plot_data/plot_data[,1]
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", ylim=c(0.95,1.05))
legend(60,1.03, names, col=1:4, lty=1:4)

d <- data.frame(as.matrix(fread(paste(data.place,"\\Add_EnSim_rep1000000_SepSearch_SepDis22_sec180.csv", sep=""), header=T)))
time <- seq(0.2, 180, 0.2)
names <- c("sim", "Pause", "wo pause")
plot_data <- matrix(as.numeric(t(d[,3:length(d[1,])])), ncol=3)
plot_data <- plot_data/plot_data[,1]
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", ylim=c(0.95,1.05))
legend(60,1.03, names, col=1:4, lty=1:4)


