## codes for analysis on separation search
## data from R.speratus
## 170619 N. Mizumoto

## packages
library(data.table)
library(exactRankTests)

## function
se  <-  function(x)
{
  y  <-  x[!is.na(x)]  #  remove  the  missing  values
  sqrt(var(as.vector(y))/length(y))
}

## setting
anal_FPS <- 1

set1X <- 0.38522848
set1Y <- 0.384615385
set2X <- 0.379879487
set2Y <- 0.381679389
set3X <- 0.219663687
set3Y <- 0.219597153

##########
##### main
##########
## convert data (solo walk)
## 設定したFPSに調整、サイズ等スケーリングし、txtファイルで出力
setwd("F:\\tandem\\solo-tandem\\Rsperatus\\location_data\\raw_data\\")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".csv",full.names=T)
for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  
  #########################
  ##### data を解析用に修正・保存
  ########################
  d <- data.frame(as.matrix(fread(data.path, header=T)))
  
  ## 時間合わせ
  # 30FPSのdataを変更
  d[,1] <- d[,1]/30
  anal_time <- seq(1, 30*30*60+1, 30/anal_FPS)
  
  
  ## scaling
  # 満遍なく移動したと想定（max - min = 150mm）
  Xmin <- min(na.omit(d$x0)); d$x0 <- d$x0 - Xmin;
  Ymin <- min(na.omit(d$y0)); d$y0 <- d$y0 - Ymin;
  setpoint <- regexpr("set", data.path)[1]
  videoset <- as.numeric(substr(data.path, setpoint+3, setpoint+3))
  if(videoset==1){
    d$x0 <- d$x0 * set1X; d$y0 <- d$y0 * set1Y;
  } else if(videoset==2){
    d$x0 <- d$x0 * set2X; d$y0 <- d$y0 * set2Y;
  } else {
    d$x0 <- d$x0 * set3X; d$y0 <- d$y0 * set3Y;
  }
  
  
  ## データ作成, 軌跡保存
  da <- d[anal_time,]
  
  L <- length(da[,1])
  time <- da[,1]
  x <- da[,2]
  y <- da[,3]
  dis <- ( (x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2 )^0.5
  speed <- dis*5
  
  Ax <- (x[3:L-1] - x[3:L-2])
  Bx <- (x[3:L] - x[3:L-1])
  Ay <- (y[3:L-1] - y[3:L-2])
  By <- (y[3:L] - y[3:L-1])
  hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
  cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
  angle <- acos(cos)*hugo
  
  ## save
  ## sex, 0:F, 1:M
  datepoint <- regexpr("170", data.path)[1]
  id <- (substr(data.path, datepoint, datepoint+13))
  date <- (substr(data.path, datepoint, datepoint+5))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+13, datepoint+13))
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  trajectory_name <- gsub(".csv",".pdf",data.path)
  trajectory_name <- gsub("after","zafter",trajectory_name)
  
  #pdf(trajectory_name)
  
  par(pin=c(2.2, 2.2), mfrow=c(2,2))
  #par(pin=c(5, 5), mfrow=c(1,1))
  plot(y0~x0, data=da, type="l",xlim=c(0,145), ylim=c(0,145), xlab="x (mm)", ylab="y (mm)", main=paste(id,treat))
  plot(time[2:L], speed, xlim=c(0,1800), ylim=c(0,80), xlab="time (sec)", ylab="speed (mm/sec)", main=paste(id,treat))
  plot(time[3:L], abs(angle), xlim=c(0,1800), ylim=c(0,3.15),  xlab="time (sec)", ylab="angle (rad)", main=paste(id,treat))
  
  #dev.off()
  
  output_name <- gsub(".csv",paste("_", anal_FPS, "FPS.txt", sep=""),data.path)
  output <- data.frame(time = time, x = x, y = y, speed = c(NA,speed), angle = c(NA,NA,angle))
  write.table(output, output_name)
  
  #for(i in seq(1, length(anal_time), 10)){
  #  j <- max(1,i-50)
  #  plot(y0~x0, data=da[j:i,], type="l", xlim=c(0,150), ylim=c(0,150), main = paste(da[i,1], "sec"))
  #  points(y0~x0, data=da[i,], pch=19)
  #  Sys.sleep(0.05)
  #}
}


## convert data (tandem walk)
## 設定したFPSに調整、サイズ等スケーリングし、txtファイルで出力
## 注目している性のみを出力
setwd("F:\\tandem\\solo-tandem\\Rsperatus\\location_data\\rawdata_tandem")
anal_FPS <- 1
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".csv",full.names=T)
for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  
  #########################
  ##### data を解析用に修正・保存
  ########################
  d <- data.frame(as.matrix(fread(data.path, header=T)))
  
  ## 時間合わせ
  # 30FPSのdataを変更
  d[,1] <- d[,1]/30
  anal_time <- seq(1, 30*10*60+1, 30/anal_FPS)
  
  ## sex, 0:F, 1:M
  datepoint <- regexpr("170", data.path)[1]
  id <- (substr(data.path, datepoint, datepoint+13))
  date <- (substr(data.path, datepoint, datepoint+5))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+13, datepoint+13))
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  
  ## scaling
  # 満遍なく移動したと想定（max - min = 150mm）
  Xmin <- min(na.omit(d$x0)); d$x0 <- d$x0 - Xmin;
  Ymin <- min(na.omit(d$y0)); d$y0 <- d$y0 - Ymin;
  Xmin <- min(na.omit(d$x1)); d$x1 <- d$x1 - Xmin;
  Ymin <- min(na.omit(d$y1)); d$y1 <- d$y1 - Ymin;
  setpoint <- regexpr("set", data.path)[1]
  videoset <- as.numeric(substr(data.path, setpoint+3, setpoint+3))
  if(videoset==1){
    d$x0 <- d$x0 * set1X; d$y0 <- d$y0 * set1Y;
    d$x1 <- d$x1 * set1X; d$y1 <- d$y1 * set1Y;
  } else if(videoset==2){
    d$x0 <- d$x0 * set2X; d$y0 <- d$y0 * set2Y;
    d$x1 <- d$x1 * set2X; d$y1 <- d$y1 * set2Y;
  } else {
    d$x0 <- d$x0 * set3X; d$y0 <- d$y0 * set3Y;
    d$x1 <- d$x1 * set3X; d$y1 <- d$y1 * set3Y;
  }
  
  ## データ作成, 軌跡保存
  da <- d[anal_time,]
  
  L <- length(da[,1])
  time <- da[,1]
  
  if(sex == 0){
    x <- da$x0; y <- da$y0;
    xc <- da$x1; yc <- da$y1;
  } else {
    x <- da$x1; y <- da$y1; 
    xc <- da$x0; yc <- da$y0;
  }
  
  dis <- ( (x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2 )^0.5
  speed <- dis*anal_FPS
  
  dis2 <- ( (xc[2:L]-xc[2:L-1])^2 + (yc[2:L]-yc[2:L-1])^2 )^0.5
  speed2 <- dis2*anal_FPS
  
  plot(speed,speed2)
  
  Ax <- (x[3:L-1] - x[3:L-2])
  Bx <- (x[3:L] - x[3:L-1])
  Ay <- (y[3:L-1] - y[3:L-2])
  By <- (y[3:L] - y[3:L-1])
  hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
  cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
  angle <- acos(cos)*hugo
  
  ## save
  
  
  trajectory_name <- gsub(".csv",".pdf",data.path)
  
  pdf(trajectory_name)
  
  par(pin=c(2.2, 2.2), mfrow=c(2,2))
  #par(pin=c(5, 5), mfrow=c(1,1))
  plot(y0~x0, data=da, type="l",xlim=c(0,145), ylim=c(0,145), xlab="x (mm)", ylab="y (mm)", main=paste(id,treat),col=2)
  plot(y1~x1, data=da, type="l",xlim=c(0,145), ylim=c(0,145), xlab="x (mm)", ylab="y (mm)", main=paste(id,treat))
  plot(time[2:L], speed, xlab="time (sec)", ylab="speed (mm/sec)", main=paste(id,treat))
  plot(time[3:L], abs(angle), xlab="time (sec)", ylab="angle (rad)", main=paste(id,treat))
  
  dev.off()
  
  output_name <- gsub(".csv",paste("_",anal_FPS,"FPS.txt",sep=""),data.path)
  output <- data.frame(time = time, x = x, y = y, speed = c(NA,speed), angle = c(NA,NA,angle))
  write.table(output, output_name)
  
  #for(i in seq(1, length(anal_time), 10)){
  #  j <- max(1,i-50)
  #  plot(y0~x0, data=da[j:i,], type="l", xlim=c(0,150), ylim=c(0,150), main = paste(da[i,1], "sec"))
  #  points(y0~x0, data=da[i,], pch=19)
  #  Sys.sleep(0.05)
  #}
}



## 1分ごとの平均値の解析
######################
## speed, angle
Colony = Sex = ID = Minutes = Treat = Mean_speed = Mean_angle <- NULL
mincom <- seq(1,30,1)
setwd("F:\\tandem\\solo-tandem\\Rsperatus\\location_data\\5FPS")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  d <- read.table(data.path, header=T)
  
  datepoint <- regexpr("170", data.path)[1]
  id <- (substr(data.path, datepoint, datepoint+13))
  date <- (substr(data.path, datepoint, datepoint+5))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+13, datepoint+13))
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  ## 1分ごとの平均値
  L <- length(d[,1])
  minutes <- c(0, rep(seq(1,30,1), each = 5*60))
  mean_speed <- tapply(na.omit(d$speed), minutes[2:L][!is.nan(d$speed)&!is.na(d$speed)], mean)
  mean_angle <- tapply(na.omit(abs(d$angle)), minutes[3:L][!is.nan(d$angle)&!is.na(d$angle)], mean)
  if(length(mean_angle)<30){mean_angle<-c(mean_angle,rep(NA,30-length(mean_angle)))}
  if(length(mean_speed)<30){mean_speed<-c(mean_speed,rep(NA,30-length(mean_speed)))}
  Colony <- c(Colony, rep(colony, 30))
  Sex <- c(Sex, rep(sex, 30))
  ID <- c(ID, rep(paste(date,iter,sep="_"), 30))
  Treat <- c(Treat, rep(treat, 30))
  Minutes <- c(Minutes, mincom)
  Mean_speed <- c(Mean_speed, mean_speed)
  Mean_angle <- c(Mean_angle, mean_angle)
}  


res <- na.omit(data.frame(colony=Colony, sex=Sex, is=ID, treat=Treat, min=Minutes,
                          mean_speed=Mean_speed, mean_angle=Mean_angle))

## plot
res_before <- res[res$treat=="before",]
res_tandem <- res[res$treat=="tandem",]
res_after <- res[res$treat=="after",]
res_before_mean <- tapply(res_before$mean_speed, res_before[,c(5,2)], mean)
res_tandem_mean <- tapply(res_tandem$mean_speed, res_tandem[,c(5,2)], mean)
res_after_mean <- tapply(res_after$mean_speed, res_after[,c(5,2)], mean)
res_before_se <- tapply(res_before$mean_speed, res_before[,c(5,2)], se)
res_tandem_se <- tapply(res_tandem$mean_speed, res_tandem[,c(5,2)], se)
res_after_se <- tapply(res_after$mean_speed, res_after[,c(5,2)], se)

par(mfrow = c(2,3), pin=c(2.2,2.2))
plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=c(0,16))
arrows(mincom, res_before_mean-res_before_se, mincom, res_before_mean+res_before_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(mincom, res_before_mean, xlim=c(1,30), ylim=c(0,16),
        cex = 1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="speed (mm/sec)", xlab="time (min)", main="before")

plot(0, ann=F, axes=F, type="n", xlim=c(1,10), ylim=c(0,16))
arrows(1:10, res_tandem_mean-res_tandem_se, 1:10, res_tandem_mean+res_tandem_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(1:10, res_tandem_mean, xlim=c(1,10), ylim=c(0,16), 
        cex = 1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="speed (mm/sec)", xlab="time (min)", main="tandem")

plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=c(0,16))
arrows(mincom, res_after_mean-res_after_se, mincom, res_after_mean+res_after_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(res_after_mean, xlim=c(1,30), ylim=c(0,16), 
        cex=1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="speed (mm/sec)", xlab="time (min)", main="after")


res_before_mean <- tapply(res_before$mean_angle, res_before[,c(5,2)], mean)
res_tandem_mean <- tapply(res_tandem$mean_angle, res_tandem[,c(5,2)], mean)
res_after_mean <- tapply(res_after$mean_angle, res_after[,c(5,2)], mean)
res_before_se <- tapply(res_before$mean_angle, res_before[,c(5,2)], se)
res_tandem_se <- tapply(res_tandem$mean_angle, res_tandem[,c(5,2)], se)
res_after_se <- tapply(res_after$mean_angle, res_after[,c(5,2)], se)

plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=c(0,1.5))
arrows(mincom, res_before_mean-res_before_se, mincom, res_before_mean+res_before_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(mincom, res_before_mean, xlim=c(1,30), ylim=c(0,1.5), 
        cex=1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="before")

plot(0, ann=F, axes=F, type="n", xlim=c(1,10), ylim=c(0,1.5))
arrows(1:10, res_tandem_mean-res_tandem_se, 1:10, res_tandem_mean+res_tandem_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(1:10, res_tandem_mean, xlim=c(1,10), ylim=c(0,1.5), 
        cex = 1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="tandem")

plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=c(0,1.5))
arrows(mincom, res_after_mean-res_after_se, mincom, res_after_mean+res_after_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(res_after_mean, xlim=c(1,30), ylim=c(0,1.5), 
        cex=1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="after")


## statistical analysis
# speed
res_before_mean <- tapply(res_before$mean_speed, res_before[,c(5,2)], mean)
res_tandem_mean <- tapply(res_tandem$mean_speed, res_tandem[,c(5,2)], mean)
res_after_mean <- tapply(res_after$mean_speed, res_after[,c(5,2)], mean)
res_before_se <- tapply(res_before$mean_speed, res_before[,c(5,2)], se)
res_tandem_se <- tapply(res_tandem$mean_speed, res_tandem[,c(5,2)], se)
res_after_se <- tapply(res_after$mean_speed, res_after[,c(5,2)], se)

P <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_after[res_after$min == i,]
  P[i] <- wilcox.exact(mean_speed ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*90,4)

P <- rep(NULL, 30)
P2 <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_before[res_before$min == i,]
  P[i] <- wilcox.exact(mean_speed ~ sex, data=analysis, paired=F)[3]
  P2[i] <- var.test(mean_speed ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*90,4)



P <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_after[res_after$min == i,]
  P[i] <- wilcox.exact(mean_angle ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*90,4)

P <- rep(NULL, 30)
P2 <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_before[res_before$min == i,]
  P[i] <- wilcox.exact(mean_angle ~ sex, data=analysis, paired=F)[3]
  P2[i] <- var.test(mean_angle ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*90,4)

################



## タンデムbeforeの後ろ3分, afterの前3分に着目
## この部分について詳細に解析を行う

######################

## setting 
anal_FPS <- 5

## data
setwd("F:\\tandem\\solo-tandem\\Rsperatus\\location_data\\5FPS")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

## looping behaviorに着目
j <- 2
data.path <- f.namesplace[j]
d <- read.table(data.path, header=T)
datepoint <- regexpr("170", data.path)[1]
id <- (substr(data.path, datepoint, datepoint+13))
date <- (substr(data.path, datepoint, datepoint+5))
if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
colony <- (substr(data.path, datepoint+7, datepoint+9))
iter <- (substr(data.path, datepoint+13, datepoint+13))
treat <- (substr(data.path, datepoint+15, datepoint+20))
if(treat=="after-"){treat <- "after"}

# 3 minutes data
data <- d[1:(60*3*anal_FPS+1),]

# 軌跡の図示
plot(speed ~ time, data=data)
plot(y ~ x, data=data, type="l")

# 部分に注目 (20 sec 分)
par(pin=c(3,3))
plot(y ~ x, data=data[1:110,], type="l")

# 図1
# looping behavior
# そのままstep lengthをとるのは不可能
# なので、loopingを他の部分と識別したい。


# 角度を計算
d3 <- data[1:200,]
d3 <- data[400:600,]
x <- d3[,2]
y <- d3[,3]
L <- length(x)

Ax <- (x[3:L-1] - x[3:L-2])
Bx <- (x[3:L] - x[3:L-1])
Ay <- (y[3:L-1] - y[3:L-2])
By <- (y[3:L] - y[3:L-1])
hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
angle <- acos(cos)*hugo

d3$angle <- c(NA,angle,NA)

d3

diff(d3$angle)
plot(y ~ x, data=d3, type="o")
angle_diff <- c(diff(d3$angle),NA)
plot(angle)
plot(angle_diff)
loop <- 2+as.numeric(abs(angle_diff) <= pi * 40/180)
d3 <- cbind(d3,angle_diff,loop)
# 図2

plot(y ~ x, data=d3, type="o", col=(2-d3$loop))
plot(d3$y ~ d3$x, type="o", col=(2-d3$loop))
par(pin=c(4,3))
plot(angle ~ time, data=d3[3:109,], type="o", col=2-as.numeric(loop))
plot(y~x, data=na.omit(d3), type="o", col=2-loop)

for(i in 3:length(d3[,1])){
  plot(y~x, data=na.omit(d3[3:i,]), type="o", col=2-loop)
  Sys.sleep(0.1)
}

# ・角度の差がない
# ・角度の符号が同じ
# ・360度以上回転する

angle
loop
for(i in 1:length(angle)){
  if()
  
}