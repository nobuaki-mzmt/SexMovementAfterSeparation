## Sex movement after separation

## タンデム実験
## 解析ノート 2017/7/20

## pre-analysis
## data modificationの方針決め

######################

## packages
library(data.table)

## setting 
anal_FPS <- 5

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\tandem\\raw_data")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".csv",full.names=T)

ID = sep_time = sep_male_speed = sep_female_speed <- NULL

for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  d <- data.frame(as.matrix(fread(data.path, header=T)))
  datepoint <- regexpr("170", data.path)[1]
  id <- (substr(data.path, datepoint, datepoint+11))
  date <- (substr(data.path, datepoint, datepoint+5))
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+11, datepoint+11))
  
  head(d)
  
  d[,1] <- d[,1]/30
  anal_time <- seq(1, 30*30*60+1, 30/anal_FPS)
  
  #plot(y0~x0,data=d, type="l")
  #plot(y1~x1,data=d, type="l")
  
  ## scaling
  # 満遍なく移動したと想定（max - min = 150mm）#あとから修正入るかも 7/20
  Xmin <- min(c(d$x0,d$x1)); Xmax <- max(c(d$x0,d$x1)); 
  Ymin <- min(c(d$y0,d$y1)); Ymax <- max(c(d$y0,d$y1)); 
  XRange <- Xmax-Xmin;
  YRange <- Ymax-Ymin;
  d$x0 <- (d$x0-Xmin)*145/XRange
  d$x1 <- (d$x1-Xmin)*145/XRange
  d$y0 <- (d$y0-Ymin)*145/YRange
  d$y1 <- (d$y1-Ymin)*145/YRange
  
  da <- d[anal_time,]
  
  head(da)
  
  time <- da$position
  Male_X <- da$x1
  Male_Y <- da$y1
  Female_X <- da$x0
  Female_Y <- da$y0
  Male_speed <- sqrt ( diff(Male_X)^2 + diff(Male_Y)^2 ) * 5
  Female_speed <- sqrt ( diff(Female_X)^2 + diff(Female_Y)^2 ) * 5
  
  
  MF_dis <- sqrt( (Male_X-Female_X)^2 + (Male_Y-Female_Y)^2 )
  
  # plot(Female_X, Female_Y, type="l")
  #plot(MF_dis)
  
  ## ためしにオスメス間の距離が10mm以上離れるとseparateとする
  ## tandemは距離6mm以下になるとtandem再開とする
  separate <- MF_dis > 10
  tandem <- MF_dis < 6
  # plot(separate)
  mean(Male_speed[separate])
  mean(Female_speed[separate])
  
  L <- length(time)
  
  # separateイベントにラベル付け
  count <- 0
  in_sep <- 0
  sep_label <- rep(0,L)
  for(i in 1:L){
    if(in_sep > 0){
      if(tandem[i]){
        count <- count + 1
        in_sep <- 0
      } else {
        sep_label[i] <- count
      }
    } else {
      if(separate[i]){
        sep_label[i] <- count
        in_sep <- 1
      }
    }
  }
  
  if(in_sep == 1){ sep_label[sep_label == count] <- 0}
  if(count>2){
    for(i in 1:(count-1)){
      sep_time <- c(sep_time, length(MF_dis[sep_label == i]) / 5 )
      sep_male_speed <- c(sep_male_speed, mean(Male_speed[sep_label == i]) )
      sep_female_speed <- c(sep_female_speed, mean(Female_speed[sep_label == i]) )
    }
    ID <- c(ID, rep(j, count-1))
  }
  #plot(sep_female_speed, sep_time, main=id)
}

plot(sep_female_speed, sep_time,  xlim=c(0,20), ylim=c(0,400))
cbind(ID,sep_time, sep_female_speed)














plot(separate, type="l")
plot(time, tandem, type="l")
plot(time[1000:1500], tandem[1000:1500], type="l")

plot(Male_X[2700:3000], Male_Y[2700:3000], xlim=c(0,145), ylim=c(0,145))
par(new=T)
plot(Female_X[2700:3000], Female_Y[2700:3000], col=2, xlim=c(0,145), ylim=c(0,145))

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




Colony = Sex = ID = Minutes = Treat <- NULL
mincom <- seq(1,30,1)
minutes <- c(0, rep(seq(1,30,1), each = 5*60))
stoptime <- NULL
for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  d <- read.table(data.path, header=T)
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  date <- (substr(data.path, datepoint, datepoint+5))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+13, datepoint+13))
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  time <- d[,1]
  x <- d[,2]
  y <- d[,3]
  
  # 止まっている時間をみる
  L <- length(x)
  dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
  
  
  speed <- dis*anal_FPS
  
  Ax <- (x[3:L-1] - x[3:L-2])
  Bx <- (x[3:L] - x[3:L-1])
  Ay <- (y[3:L-1] - y[3:L-2])
  By <- (y[3:L] - y[3:L-1])
  hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
  cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
  angle <- c(NA, acos(cos)*hugo, NA)
  
  # 0.2秒の移動距離が0.2以下(秒速1mm以下)の場合を止まっているとみなす
  stop <- (dis < 0.2)
  stop[1] <- FALSE
  
  cbind(x,y,dis,stop,angle)
  
  
  # stopを除いたデータを使って、角度を計算する
  time2 <- time[!stop]
  x2 <- x[!stop]
  y2 <- y[!stop]
  L2 <- length(x2)
  
  Ax <- (x2[3:L2-1] - x2[3:L2-2])
  Bx <- (x2[3:L2] - x2[3:L2-1])
  Ay <- (y2[3:L2-1] - y2[3:L2-2])
  By <- (y2[3:L2] - y2[3:L2-1])
  hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
  cos <- round((Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5,14)
  kari_angle <- c(NA, acos(cos)*hugo, NA)
  
  tra_angle <- rep(NA, ds)
  for(i in 1:L){  tra_angle[time==time2[i]] <- kari_angle[i] }
  
  gsub("170623", "170720", data.path)
  
  #output <- data.frame(time = time, x = x, y = y, dis = dis, stop=stop, speed = speed, angle = angle, tra_angle=tra_angle)
  #write.table(output, output_name)
  
  
  ## MSD (time)
  
  Tau <- c(1:9, 10*(1:9), 100*(1:9))
  lt <- length(Tau)
  Msd <- rep(0,lt)
  
  for(j in 1:lt){
    point <- seq(1,L,Tau[j]*anal_FPS)  # 5FPSなのでdata構造に合わせている
    displacement <- rep(0,length(point))
    beginX = beginY = count <- 0
    for(i in point){
      moved <- (((x[i]-beginX)^2 + (y[i]-beginY)^2))^0.5
      beginX <- x[i]
      beginY <- y[i]
      count = count + 1
      displacement[count] <- moved 
    }
    Msd[j] <- mean(displacement^2)
  }
  #pdf(paste("D:",name,".pdf",sep=""))
  
  plot(log10(Tau), log10(Msd), col=2, type="l", lwd=2)
  a <- min(log10(Msd))
  abline(a=a, b=1, lty=2)
  abline(a=a, b=2, lty=2)
  r <- lm(log10(Msd) ~ log10(Tau))
  r
  
  
  
  
  i<-1
  x3 <- x[(i-1)*60<time & time<=60*i]
  y3 <- y[(i-1)*60<time & time<=60*i]
  x3 <- x[(i-1)*60<time & time<=300*i]
  y3 <- y[(i-1)*60<time & time<=300*i]
  L3 <- length(x3)
  plot(y3~x3, type="o", col=2-stop)
  lt <- length(Tau)
  Msd <- rep(0,lt)
  
  for(j in 1:lt){
    point <- seq(1,L3,Tau[j]*anal_FPS)  # 5FPSなのでdata構造に合わせている
    displacement <- rep(0,length(point))
    beginX = beginY = count <- 0
    for(i in point){
      moved <- (((x3[i]-beginX)^2 + (y3[i]-beginY)^2))^0.5
      beginX <- x3[i]
      beginY <- y3[i]
      count = count + 1
      displacement[count] <- moved 
    }
    Msd[j] <- mean(displacement^2)
  }
  #pdf(paste("D:",name,".pdf",sep=""))
  
  plot(log10(Tau), log10(Msd), col=2, type="l", lwd=2)
  a <- min(log10(Msd))
  abline(a=a, b=1, lty=2)
  abline(a=a, b=2, lty=2)
  r <- lm(log10(Msd) ~ log10(Tau))
  r
  
  
  stop
  hist(tra_angle)
  hist(angle)
  angle
  tra_angle
  
  plot(time[1:(60*5)], angle[1:(60*5)], type="l")
  angle[1:(60*5)]
  plot(na.omit(tra_angle[1:(60*5)]), type="l")
  
  da <- abs(na.omit(tra_angle[1:(60*5)]))
  acf(da[abs(da)>0.6])
  180*0.6/pi
  
  hist(da[da>0.6])
  
  plot(x[1:(60*5)], y[1:(60*5)], type="o", cex=  as.numeric(abs(angle[1:(60*5)])>0.5))
  