## Sex movement after separation

## 解析ノート 2017/7/20

## 生データからnoiseキャンセルを行い、
## speed, diffusivity, angle, stopのデータを加えて出力
## その後1分置きのデータをとる


######################

## packages
library(data.table)
library(exactRankTests)
library(CircStats)

## functions
se  <-  function(x)
{
  y  <-  x[!is.na(x)]  #  remove  the  missing  values
  sqrt(var(as.vector(y))/length(y))
}


## setting 
anal_FPS <- 5

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

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
}






## 1分ごとの平均値の解析
######################
## speed, angle

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

Colony = Sex = ID = Minutes = Treat = Mean_speed = Mean_angle <- NULL
mincom <- seq(1,30,1)

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
# 

# U-test
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
  #P[i] <- kruskal.test(mean_speed ~ colony, data=analysis)[3]
}
round(as.numeric(P)*70,4)

P <- rep(NULL, 30)
P2 <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_before[res_before$min == i,]
  P[i] <- wilcox.exact(mean_speed ~ sex, data=analysis, paired=F)[3]
  P2[i] <- var.test(mean_speed ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*70,4)



P <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_after[res_after$min == i,]
  P[i] <- wilcox.exact(mean_angle ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*70,4)

P <- rep(NULL, 30)
P2 <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_before[res_before$min == i,]
  P[i] <- wilcox.exact(mean_angle ~ sex, data=analysis, paired=F)[3]
  P2[i] <- var.test(mean_angle ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*70,4)

################



#############
## MSD (time)
#############

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

analyze_time <- 1 # min
an_sec <- analyze_time*60
par(pin=c(4,4), mfrow=c(1,1))
plot(0, type="n", axes=F, ann=F, xlim=c(0,log10(60*analyze_time)), ylim=c(-2, 4.5))
axis(1, at=0:3, labels=c(1, 10, 100, 1000))
axis(2, at=seq(-2,4,2), labels=c("10^-2", "1", "10^2", "10^4"), las=1)
box()

Tau <- seq(1, an_sec/2, 1)
#Tau <- c(1:9, 1:9*10, 1:(an_sec/100)*100)
if(round(an_sec/100) - an_sec/100 != 0){ Tau <- c(Tau, an_sec) }

lt <- length(Tau)
MSD_mean <- rep(0,lt)
for(k in 1:length(f.namesplace)){
  data.path <- f.namesplace[k]
  d <- read.table(data.path, header=T)
  
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  date <- (substr(data.path, datepoint, datepoint+5))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+13, datepoint+13))
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}

  if(treat == "before" || treat=="tandem"){next;}
  if(sex == 1){next;}
    
  x <- d[1:(analyze_time*60*anal_FPS+1),2]
  y <- d[1:(analyze_time*60*anal_FPS+1),3]
  #dis <- sqrt(diff(x)^2 + diff(y)^2)
  #plot(dis) 


  Msd <- rep(0,lt)
  
  L <- length(x)
  
  for(j in 1:lt){
    point <- seq(1,L,Tau[j]*anal_FPS)  # 5FPSなのでdata構造に合わせている
    displacement <- rep(0,length(point))
    beginX <- x[point[1]]
    beginY <- y[point[1]]
    count <- 0
    for(i in point[2:length(point)]){
      moved <- (((x[i]-beginX)^2 + (y[i]-beginY)^2))^0.5
      beginX <- x[i]
      beginY <- y[i]
      count = count + 1
      displacement[count] <- moved 
    }
    Msd[j] <- mean(displacement^2)
  }
  MSD_mean <- MSD_mean + Msd
  points(log10(Tau), log10(Msd), col="grey", type="o", lwd=2)
  a <- min(log10(Msd))
  #abline(a=a, b=1, lty=2)
  #abline(a=a, b=2, lty=2)
  r <- lm(log10(Msd) ~ log10(Tau))
  r
}

points(log10(Tau), log10(MSD_mean/37), col=2, type="o", lwd=2)

a <- min(log10(MSD_mean/37))
abline(a=a, b=1, lty=2)
abline(a=a, b=2, lty=2)

sd_list <- rep(0,lt-1)
for(i in 1:(lt-1)){
  r1 <- lm(log10(MSD_mean/37)[1:i] ~ log10(Tau)[1:i])
  r2 <- lm(log10(MSD_mean/37)[(i+1):lt] ~ log10(Tau)[(i+1):lt])
  var1 <- sum((log10(MSD_mean/37)[1:i] - (r1$coefficients[1] + log10(Tau[1:i]) * r1$coefficients[2]))^2)
  var2 <- sum((log10(MSD_mean/37)[(i+1):lt] - (r2$coefficients[1] + log10(Tau[(i+1):lt]) * r2$coefficients[2]))^2)
  sd_list[i] <- sqrt(var1+var2)/lt
}
#plot(sd_list)
sd_list
i<-10
r1 <- lm(log10(MSD_mean/37)[1:i] ~ log10(Tau)[1:i])
r2 <- lm(log10(MSD_mean/37)[(i+1):lt] ~ log10(Tau)[(i+1):lt])
abline(r1, col=3)
abline(r2, col=3)




#############
## MSD (dis)
#############

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

analyze_time <- 1 # min
an_sec <- analyze_time*60
par(pin=c(4,4), mfrow=c(1,1))
plot(0, type="n", axes=F, ann=F, xlim=c(0,log10(60*analyze_time)), ylim=c(-2, 4.5))
axis(1, at=0:3, labels=c(1, 10, 100, 1000))
axis(2, at=seq(-2,4,2), labels=c("10^-2", "1", "10^2", "10^4"), las=1)
box()


for(k in 1:length(f.namesplace)){
  data.path <- f.namesplace[k]
  d <- read.table(data.path, header=T)
  
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  date <- (substr(data.path, datepoint, datepoint+5))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+13, datepoint+13))
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  #if(treat == "before" || treat=="tandem"){next;}
  #if(sex == 1){next;}
  
  x <- d[1:(analyze_time*60*anal_FPS+1),2]
  y <- d[1:(analyze_time*60*anal_FPS+1),3]
  dis <- sqrt(diff(x)^2 + diff(y)^2)
  #plot(dis) 
  TotalDis <- sum(dis)
  Ax <- (x[3:L-1] - x[3:L-2])
  Bx <- (x[3:L] - x[3:L-1])
  Ay <- (y[3:L-1] - y[3:L-2])
  By <- (y[3:L] - y[3:L-1])
  hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
  cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
  angle <- c(NA, acos(cos)*hugo, NA)
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
  
  tra_angle <- rep(NA, L2)
  for(i in 1:L2){  tra_angle[time==time2[i]] <- kari_angle[i] }
  
  Tau <- seq(1, TotalDis, 1)
  #Tau <- c(1:9, 1:9*10, 1:(an_sec/100)*100)
  #if(round(an_sec/100) - an_sec/100 != 0){ Tau <- c(Tau, an_sec) }
  
  lt <- length(Tau)
  
  Msd <- rep(0,lt)
  
  dis <- c(0,dis)
  L <- length(x)
  
  for(j in 1:lt){
    displacement <- rep(0,sum(dis)/Tau[j])
    count <- 1
    moved <- 0
    move <- 0
    cycle <- 0
    beginX <- x[1]
    beginY <- y[1]
    for(i in 1:L){
      at.moved <- dis[i]
      moved <- moved + at.moved
      if(moved >= Tau[j]+0.00000000001){
        while(moved >= Tau[j]+0.00000000001){
          cycle = 0
          # 軌跡上の点の間の移動は等速直線運動と仮定し、tauに到達時点の点(ax,ay)を想定
          ax <- x[i-1] + (x[i]-x[i-1]) * ((Tau[j]-(moved-at.moved)+cycle)/at.moved)
          ay <- y[i-1] + (y[i]-y[i-1]) * ((Tau[j]-(moved-at.moved)+cycle)/at.moved)
          displacement[count] <- ((ax-beginX)^2 + (ay-beginY)^2)^0.5
          count <- count +1
          beginX <- ax
          beginY <- ay 
          moved <- (((x[i]-ax)^2 + (y[i]-ay)^2))^0.5
          cycle = cycle + 1
        }	
      }
    }
    Msd[j] <- mean(displacement^2)
  }
  
  plot(log10(Tau), log10(Msd), type="l", lwd=2)
  a <- min(log10(Msd))
  abline(a=a, b=1, lty=2)
  abline(a=a, b=2, lty=2)
  abline(v = log10(145))
  abline(v = log10(mean(dis)))
  abline(h = log10(145^2))
  
  hist(angle)
  
  angle <- na.omit(angle)
  angle[angle < 0]
  wrpcauchy.ml(angle, 0, 0, acc=1e-015)
  
  
  #plot(y~x, data=data, type="l", col=2-sexnum)
  points(log10(Tau), log10(Msd), col=2-sexnum, type="l", lwd=2)

  r <- lm(log10(Msd) ~ log10(Tau))
  r
  #abline(a=r[1], b=r[2], col=2)
  
}

  




#############
## stop time
#############

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

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
  
  # stopの長さと回数を計測
  count <- 0
  in_stop <- 0
  stop_label <- rep(NA,L)
  for(i in 1:L){
    if(stop[i]){
      if(in_stop==0){
        count <- count + 1
      }
      in_stop <- 1
      stop_label[i] <- count
    } else {in_stop <- 0}
  }
  
  stop_sec <- as.vector(tapply(stop, stop_label, sum))*0.2
  max(stop_sec)
  truehist(stop_sec)
  
  r <- as.numeric(fitdistr(stop_sec, densfun  = "exponential")[1])
  lines(0:120,dexp(0:120, rate=r))
  
  
  
  myu <- PL_mle(stop_sec, 0.2)
  idea.x <- seq(0,0.9999,0.0001)
  y1.est <- PL(idea.x, myu, 0.2)
  result1 <- ks.test(stop_sec, y1.est)
  
  lambda.est <- EXP_mle(stop_sec, 0.2)
  idea.x <- seq(0,0.9999,0.0001)
  y2.est <- EXP(idea.x, lambda.est, 0.2)
  result2 <- ks.test(stop_sec, y2.est)
  
  rank <- seq(length(stop_sec),1,-1)
  plot(sort(stop_sec),rank,log="xy",xlab="value",ylab="rank")

  idea.x <- seq(0,0.9999,1/length(stop_sec))
  y1.est <- PL(idea.x, myu, 0.2)
  y2.est <- EXP(idea.x, lambda.est, 0.2)
  lines(y2.est, seq(length(y2.est), 1, -1),col="blue", lwd = 3)
  lines(y1.est, seq(length(y1.est), 1, -1),col="red", lwd = 3)
  
  
  PL_mle <- function(x,xmin){
    return(1+length(x)*(sum(log(x/xmin)))^(-1))
  }
  # 最尤推定による指数の推定値 for exponential 
  EXP_mle <- function(x,xmin){
    return(length(x)*(sum(x-xmin)^(-1)))
  }
  
  
  # power law distribution のデータを生成. rが一様乱数. aが1<a<3の指数. xminが分布の最小値
  PL <- function(r,a,xmin){
    return(xmin*(1-r)^(-1/(a-1)))
  }
  # exponential distribution のデータを生成. rが一様乱数. aが指数. xminが分布の最小値
  EXP <- function(r,a,xmin){
    return(xmin-1/a*log(1-r))
  }
  
  
}