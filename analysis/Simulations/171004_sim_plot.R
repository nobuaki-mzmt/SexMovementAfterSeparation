## Results plot for SepSearch Sim
## 171004 N.Mizumoto

## main

## simulations for Rs

## Male-Female before after

## Sep search condition
data.place <- "E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Simulations\\data\\EnSim_Rs_FMBA_rep1000000_SepSearch_SepDis16.csv"
d <- read.csv(data.place, header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
time <- seq(0.2, 180, 0.2)
rep <- 1000000
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data[,c(1,4:6)]/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate")
legend(60,0.15, names[c(1,4:6)], col=1:4, lty=1:4)


## Random search 
data.place <- "E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Simulations\\data\\EnSim_Rs_FMBA_rep1000000_RanSearch_L223.csv"
d <- read.csv(data.place, header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
time <- seq(0.2, 180, 0.2)
rep <- 1000000
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data[,c(1,4)]/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate")
legend(60,0.1, names[c(1,4)], col=1:4, lty=1:4)


## Sep search condition
data.place <- "E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Simulations\\data\\EnSim_Rs_F-Surr_rep1000000_SepSearch_SepDis16.csv"
d <- read.csv(data.place, header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
time <- seq(0.2, 180, 0.2)
rep <- 1000000
plot_data <- t(d[,3:length(d[1,])])

par(pin=c(5,4)*0.6, mfrow=c(2,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate")
legend(60,0.15, names, col=1:4, lty=1:4)

plot(time[5:length(time)],plot_data[5:length(time),2]/plot_data[5:length(time),3], type="l", ylim=c(0.95,1.02), col=2, lty=2,
     xlab="time", ylab="Relative efficiency", xlim=c(1,180))
points(time[5:length(time)],plot_data[5:length(time),1]/plot_data[5:length(time),3], type="l", col=1, lty=1)
arrows(1,1,180,1, col=3, lty=3, length=0)

##
d <- read.table("clipboard",header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
time <- seq(0.2, 180, 0.2)
rep <- 100000
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4)*0.6, mfrow=c(2,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate")
legend(60,0.15, names, col=1:4, lty=1:4)

plot(time[5:length(time)],plot_data[5:length(time),1]/plot_data[5:length(time),2], type="l",  col=1, lty=2,
     xlab="time", ylab="Relative efficiency", xlim=c(1,180), ylim=c(0.99,1.01))
arrows(1,1,180,1, col=2, lty=3, length=0)
