## Results plot for SepSearch Sim
## 180926 N.Mizumoto

## main
data.place <- "E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\SciAdv\\Revision\\simulations"
time <- seq(0.2, 180, 0.2)
rep <- 1000000

## simulations for Rs
## Random search condition
d <- read.csv(paste(data.place,"\\Rs_EnSim_rep1000000_RanSearch_L223_sec180_empdata.csv", sep=""), header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(60,0.15, names, col=1:4, lty=1:4)

## Reunion search 
d <- read.csv(paste(data.place,"\\Rs_EnSim_rep1000000_SepSearch_SepDis16_sec180_empdata.csv", sep=""), header=T)
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(60,0.15, names, col=1:4, lty=1:4)


## simulations for Cf
## Random search condition
d <- read.csv(paste(data.place,"\\Cf_EnSim_rep1000000_RanSearch_L223_sec180_empdata.csv", sep=""), header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(60,0.15, names, col=1:4, lty=1:4)

## Reunion search 
d <- read.csv(paste(data.place,"\\Cf_EnSim_rep1000000_SepSearch_SepDis22_sec180_empdata.csv", sep=""), header=T)
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(60,0.15, names, col=1:4, lty=1:4)





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
legend(60,0.15, names, col=1:4, lty=1:4)

## Reunion search 
d <- read.csv(paste(data.place,"\\Rs_EnSim_rep1000000_SepSearch_SepDis16_sec180_empdata.csv", sep=""), header=T)
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(60,0.15, names, col=1:4, lty=1:4)


## simulations for Cf
## Random search condition
d <- read.csv(paste(data.place,"\\Cf_EnSim_rep1000000_RanSearch_L223_sec180_empdata.csv", sep=""), header=T)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(60,0.15, names, col=1:4, lty=1:4)

## Reunion search 
d <- read.csv(paste(data.place,"\\Cf_EnSim_rep1000000_SepSearch_SepDis22_sec180_empdata.csv", sep=""), header=T)
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(60,0.15, names, col=1:4, lty=1:4)




###################
### check Rs reverse point

d <- data.frame(as.matrix(fread(paste(data.place,"\\EnSim_rep10000_SepSearch_SepDis16_sec50000_empdata.csv", sep=""), header=T)))
d <- data.frame(as.matrix(fread(paste(data.place,"\\EnSim_rep100000_SepSearch_SepDis16_sec100000_empdata.csv", sep=""), header=T)))
x <- c( c(2,5,10)*0.1, c(2,5,10), c(2,5,10)*10, c(2,5,10)*100, c(2,5,10)*1000, c(2,5,10)*10000)*5
x <- c( 2:10*0.1, 2:10, 2:10*10, 2:10*100, 2:10*1000, 2:10*10000)
matplot(x, t(d[,x*5+2]), log="x", type="l", axes=F, las=1,
        xlab="time (sec)", ylab="Encounter rate",)
box()
axis(1, x)
axis(2)

time <- seq(0.2, 360, 0.2)
name <- d[,1:2]
names <- paste(name[,1],name[,2], sep=" and ")
plot_data <- t(d[,3:length(d[1,])])
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data/rep, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", log="x")
legend(60,0.15, names, col=1:4, lty=1:4)


d[,400000:500002]

# 84663 sec ?????????



####################################
## for ad sim

d <- read.table("clipboard",header=T)
time <- seq(0.2, 180, 0.2)
names <- c("sim", "Pause", "wo pause")
plot_data <- t(d[,3:length(d[1,])])
plot_data <- plot_data/plot_data[,1]
par(pin=c(5,4), mfrow=c(1,1))
matplot(time, plot_data, type="l", las=1,
        xlab="time (sec)", ylab="Encounter rate", ylim=c(0.90,1.10))
legend(60,0.15, names, col=1:4, lty=1:4)


