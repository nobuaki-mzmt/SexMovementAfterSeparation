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

angle_cal <- function(X, Y, Length){
  Ax <- (X[3:Length-1] - X[3:Length-2])
  Bx <- (X[3:Length] - X[3:Length-1])
  Ay <- (Y[3:Length-1] - Y[3:Length-2])
  By <- (Y[3:Length] - Y[3:Length-1])
  hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
  cos <- round((Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5,14)
  return(acos(cos)*hugo)
}


# power law distribution のデータを生成. rが一様乱数. aが1<a<3の指数. xminが分布の最小値
PL <- function(r,a,xmin){
  return(xmin*(1-r)^(-1/(a-1)))
}
# exponential distribution のデータを生成. rが一様乱数. aが指数. xminが分布の最小値
EXP <- function(r,a,xmin){
  return(xmin-1/a*log(1-r))
}
# stretched exponential distribution のデータを生成. rが一様乱数. aが指数. xminが分布の最小値
SEXP <- function(r, a, b, xmin){
  
}


# 最尤推定による指数の推定値 for power law
PL_mle <- function(x,xmin){
  return(1+length(x)*(sum(log(x/xmin)))^(-1))
}
# 最尤推定による指数の推定値 for exponential 
EXP_mle <- function(x,xmin){
  return(length(x)*(sum(x-xmin)^(-1)))
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




### 170819
#############
## beforeの後3分 と afterの前3分を使用して詳細比較
#############

# 注目するもの
# turning angle distribution
# speed
# stoptime

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

#####
## 各パラメーターを個体ごとに抽出・比較
## 各パラメーターをまとめたデータを作成

before <- f.namesplace[regexpr("_before",f.namesplace)>0]
tandem <- f.namesplace[regexpr("_tandem",f.namesplace)>0]
after <- f.namesplace[regexpr("_after",f.namesplace)>0]
file_table <- cbind(before, tandem, after)

par(mfrow=c(4,4), pin=c(2,1))
cols <- c("#4B6EFD", "#6264C8", "#795A93") #before, after

Colony = Sex = ID = Treat = Mu = Rho = Moved = Stop_time = True_speed <- NULL

for(k in 1:length(before)){
  for(j in c(1,3)){
    data.path <- file_table[k,j]
    d <- read.table(data.path, header=T)
    datepoint <- regexpr("170", data.path)[1]+7
    id <- (substr(data.path, datepoint, datepoint+13))
    date <- (substr(data.path, datepoint, datepoint+5))
    if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
    colony <- (substr(data.path, datepoint+7, datepoint+9))
    iter <- (substr(data.path, datepoint+13, datepoint+13))
    treat <- (substr(data.path, datepoint+15, datepoint+20))
    if(treat=="after-"){treat <- "after"}
    
    if(j==1){ anal_time <-  (60*27*anal_FPS):(60*30*anal_FPS)+1}
    if(j==3){ anal_time <-  (60*0):(60*3*anal_FPS)+1}
    
    time <- d[anal_time, 1]
    x <- d[anal_time, 2]
    y <- d[anal_time, 3]
    
    L <- length(x)
    dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
    speed <- dis*anal_FPS
    angle <- c(NA, angle_cal(x,y,L), NA)
    
    moved <- sum(dis, na.rm=T)
    
    # 止まっている時間をみる
    
    # 0.2秒の移動距離が0.2以下(秒速1mm以下)の場合を止まっているとみなす
    stop <- (dis < 0.2)
    stop[1] <- FALSE
    stop_time <- sum(stop[2:length(stop)])/5
    
    time2 <- time[!stop];  x2 <- x[!stop];  y2 <- y[!stop];  L2 <- length(x2)
    # stop時間を除く平均速度
    dis2 <- (c(NA,((x2[2:L2]-x2[2:L2-1])^2 + (y2[2:L2]-y2[2:L2-1])^2)^0.5))
    speed2 <- dis2*anal_FPS
    true_speed <- mean(speed2, na.rm=T)
    #truehist(speed2,col=cols[j], main=paste(id,treat), breaks=0:50)
    
    # stop時間を除くangleを計算
    kari_angle <- c(NA, angle_cal(x2, y2, L2), NA)
    
    tra_angle <- rep(NA, L2)
    for(i in 1:L2){  tra_angle[time==time2[i]] <- kari_angle[i] }
    
    d_anal <- data.frame(time=time, x=x, y=y, speed=speed, angle=angle, stop=stop, tra_angle=tra_angle)
    
    if(length(na.omit(tra_angle))>0){
      mle_wrpcauchy <- wrpcauchy.ml(na.omit(tra_angle), 0, 0, acc=1e-015)
      theta <- seq(-pi, pi, length=1000)
      probabl <- dwrpcauchy(theta, as.numeric(mle_wrpcauchy[1]), as.numeric(mle_wrpcauchy[2]))
      
      
      Range <- ceiling(max(probabl,na.rm=T))
      truehist(tra_angle, prob=T, breaks=seq(-pi, pi, length=100),
               ylim=c(0,Range), xlim=c(-pi,pi), col=cols[j],
               xlab="angle", ylab="density", main=paste(id,treat), axes=F)
      points(theta, probabl, type="l")  
      axis(1, at=seq(-pi,pi,length=3), label=c("-pi", "0", "pi"))
      axis(2, las=1)
      mu <- as.numeric(mle_wrpcauchy[1])
      rho <- as.numeric(mle_wrpcauchy[2])
    } else{ mu = rho = NA; plot(0)}
    
    Colony <- c(Colony, colony)
    Sex <- c(Sex, sex)
    ID <- c(ID, paste(date,iter,sep="_"))
    Treat <- c(Treat, treat)
    Moved <- c(Moved, moved)
    Stop_time <- c(Stop_time, stop_time)
    True_speed <- c(True_speed, true_speed)
    Mu <- c(Mu, mu)
    Rho <- c(Rho, rho)
  }
}

#Mu[Mu>pi] <- Mu[Mu>pi] - 2*pi

d <- data.frame(id=ID, colony=Colony, sex=Sex, treat=Treat, moved=Moved, stop_time=Stop_time, mean_speed=True_speed, Rho=Rho, Mu=Mu)

# plot data
par(pin=c(3,3), mfrow=c(2,2))
treat <- c(1,0,1,0)
sex <- c(0,0,1,1)

# moved distance
Mean <- (tapply(d[,5],d[,3:4],mean))[2:1,2:1]
SE <- (tapply(d[,5],d[,3:4],se))[2:1,2:1]

Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,3000), beside=T, las=1,
               ylab="moved distance (mm)")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topleft", c("male","female"), fill=c(1,0))　

wilcox.exact(moved ~ sex, data=d[d$treat=="after",], paired=F)
wilcox.exact(moved ~ sex, data=d[d$treat=="before",], paired=F)
wilcox.exact(moved ~ treat, data=d[d$sex==0,], paired=T)
wilcox.exact(moved ~ treat, data=d[d$sex==1,], paired=T)

# stop time
Mean <- (tapply(d[,6],d[,3:4],mean))[2:1,2:1]
SE <- (tapply(d[,6],d[,3:4],se))[2:1,2:1]

Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,180), beside=T, las=1,
               ylab="stop time (sec)")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topleft", c("male","female"), fill=c(1,0))　

wilcox.exact(stop_time ~ sex, data=d[d$treat=="after",], paired=F)
wilcox.exact(stop_time ~ sex, data=d[d$treat=="before",], paired=F)
wilcox.exact(stop_time ~ treat, data=d[d$sex==0,], paired=T)
wilcox.exact(stop_time ~ treat, data=d[d$sex==1,], paired=T)


# True speed
Mean <- (tapply(d[,7],d[,3:4],mean))[2:1,2:1]
SE <- (tapply(d[,7],d[,3:4],se))[2:1,2:1]

Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,20), beside=T, las=1,
               ylab="moving speed (mm/sec)")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topright", c("male","female"), fill=c(1,0))　

wilcox.exact(mean_speed ~ sex, data=d[d$treat=="after",], paired=F)
wilcox.exact(mean_speed ~ sex, data=d[d$treat=="before",], paired=F)
wilcox.exact(mean_speed ~ treat, data=d[d$sex==0,], paired=T)
wilcox.exact(mean_speed ~ treat, data=d[d$sex==1,], paired=T)


# Rho (wrapped caucy distribution)
Mean <- (tapply(d[,8],d[,3:4],mean, na.rm=T))[2:1,2:1]
SE <- (tapply(d[,8],d[,3:4],se))[2:1,2:1]

Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,1), beside=T, las=1,
               ylab="ρ")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topright", c("male","female"), fill=c(1,0))　

wilcox.exact(Rho ~ sex, data=d[d$treat=="after",], paired=F)
wilcox.exact(Rho ~ sex, data=d[d$treat=="before",], paired=F)
wilcox.exact(Rho ~ treat, data=d[d$sex==0 & d$id!="170601_1",], paired=T)
wilcox.exact(Rho ~ treat, data=d[d$sex==1,], paired=T)

# write.table(d, "Individual_data_summary_Rs.txt")



###################
## 170823
## 各個体30分間のデータについて、
## 止まっている時間と動いている時間について
## 指数分布・べき分布でfitしたい。

## 止まった回数とmaxの時間を求めて、
## ほとんど休止が起こらなかった個体は仕方ないので解析から外す

## 1. NA含むデータの処理
for(k in 1:length(before)){
  for(j in c(1,3)){
    data.path <- file_table[k,j]
    d <- read.table(data.path, header=T)
    anal_time <-  (60*0*anal_FPS):(60*30*anal_FPS)+1
    x <- d[anal_time, 2]
    if(sum(is.na(x))>0){print(data.path)}
  }
}

## とりあえずNA込みのデータを抜いて解析してみる


## 2. 止まった回数と最長時間を示す

par(mfrow=c(4,4), pin=c(2,1))
for(k in 1:length(before)){
  for(j in c(1,3)){
    if(k==8||k==11||k==23){next;}
    data.path <- file_table[k,j]
    d <- read.table(data.path, header=T)
    datepoint <- regexpr("170", data.path)[1]+7
    id <- (substr(data.path, datepoint, datepoint+13))
    date <- (substr(data.path, datepoint, datepoint+5))
    if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
    colony <- (substr(data.path, datepoint+7, datepoint+9))
    iter <- (substr(data.path, datepoint+13, datepoint+13))
    treat <- (substr(data.path, datepoint+15, datepoint+20))
    if(treat=="after-"){treat <- "after"}
    
    anal_time <-  (60*0*anal_FPS):(60*30*anal_FPS)+1
    time <- d[anal_time, 1]
    x <- d[anal_time, 2]
    y <- d[anal_time, 3]
    
    L <- length(x)
    dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
    stop <- (dis < 0.2)
    stop[1] <- FALSE
    
    # stopの長さと回数を計測
    count <- 0
    in_stop <- 0

    L <- length(stop)
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
    #print(max(stop_label, na.rm=T))
    
    stop_sec <- as.vector(tapply(stop, stop_label, sum))*0.2
    #print(max(stop_sec, na.rm=T))
    
    if(max(stop_sec)<6.5){next;}
    
    min_sec <- 0.2
    max_sec <- 3600

    stop_sec <- as.vector(tapply(stop, stop_label, sum))*0.2
    stop_sec <- stop_sec[stop_sec >= min_sec & stop_sec < max_sec]

    (table(stop_sec))
    stop_lab <- sort(unique(stop_sec))
    stop_count <- as.vector(table(stop_sec))
    stop_ratio <- stop_count / sum(stop_count)
    stop_cum <- rep(0, length(stop_ratio))
    for(i in 1:length(stop_ratio)){
      stop_cum[i] <- sum(stop_ratio[i:length(stop_ratio)])
    }
    
    myu <- PL_mle(stop_sec, min_sec)
    idea.x <- seq(0,0.9999,0.0001)
    y1.est <- PL(idea.x, myu, min_sec)
    result1 <- ks.test(stop_sec, y1.est)
    
    lambda.est <- EXP_mle(stop_sec, min_sec)
    idea.x <- seq(0,0.9999,0.0001)
    y2.est <- EXP(idea.x, lambda.est, min_sec)
    result2 <- ks.test(stop_sec, y2.est)
    
    # 1 culm plot
    idea.x <- seq(0,0.9999,1/length(stop_sec))
    y1.est <- PL(idea.x, myu, min_sec)
    y2.est <- EXP(idea.x, lambda.est, min_sec)
    
    #plot(stop_lab,stop_cum, log="xy", main=paste(id,treat))
    #lines(y1.est, rev(idea.x), col="red", lwd = 1)
    #lines(y2.est, rev(idea.x), col="blue", lwd = 1)
    
    # 2 rank plot
    rank <- seq(length(stop_sec),1,-1)
    plot(sort(stop_sec),rank,log="xy",xlab="value",ylab="rank", main=paste(id,treat))
    
    idea.x <- seq(0,0.9999,1/length(stop_sec))
    y1.est <- PL(idea.x, myu, min_sec)
    y2.est <- EXP(idea.x, lambda.est, min_sec)
    lines(y1.est, seq(length(y1.est), 1, -1),col="red", lwd = 1)
    lines(y2.est, seq(length(y2.est), 1, -1),col="blue", lwd = 1)

    
  }
}




d2 <- data.frame(stop=stop, stop_label=stop_label, id=id)

tapply(d2[,1], d2[,2:3], sum)*0.2

stop_sec <- as.vector(tapply(stop, stop_label, sum))*0.2
stop_sec <- stop_sec[stop_sec >= min_sec]











# stopの長さと回数を計測
count <- 0
in_move <- 0
move <- !f_bef[,"stop"]
id <- f_bef[,"id"]

L <- length(move)
id_label <- id[1]
move_label <- rep(NA,L)
for(i in 1:L){
  if(id_label==id[i]){
    if(move[i]){
      if(in_move==0){
        count <- count + 1
      }
      in_move <- 1
      move_label[i] <- count
    } else {in_move <- 0}
  } else {
    if(in_move <- 1){ move_label[move_label==count] <- NA }
    id_label <- id[i]
    count <- count + 1
    if(move[i]){ in_move = 1; move_label[i] <- count; } else { in_move = 0 }
  }
}

min_sec <- 0.2
max_sec <- 180

move_sec <- as.vector(tapply(move, move_label, sum))*0.2
max(move_sec)
truehist(move_sec)

r <- as.numeric(fitdistr(move_sec, densfun  = "exponential")[1])
lines(0:120,dexp(0:120, rate=r))


move_sec <- move_sec[move_sec > min_sec]
truehist(move_sec)
r <- as.numeric(fitdistr(move_sec, densfun  = "exponential")[1])
lines(0:120,dexp(0:120, rate=r))

myu <- PL_mle(move_sec, min_sec)
idea.x <- seq(0,0.9999,0.0001)
y1.est <- PL(idea.x, myu, min_sec)
result1 <- ks.test(move_sec, y1.est)

lambda.est <- EXP_mle(move_sec, min_sec)
idea.x <- seq(0,0.9999,0.0001)
y2.est <- EXP(idea.x, lambda.est, min_sec)
result2 <- ks.test(move_sec, y2.est)

myu.est <- TP_mle(move_sec, min_sec, max_sec)
idea.x <- seq(0,0.9999,0.0001)
y3.est <- TP(idea.x, myu.est, min_sec, max_sec)
result2 <- ks.test(move_sec, y2.est)

rank <- seq(length(move_sec),1,-1)
plot(sort(move_sec),rank,log="xy",xlab="value",ylab="rank")

idea.x <- seq(0,0.9999,1/length(move_sec))
y1.est <- PL(idea.x, myu, min_sec)
y2.est <- EXP(idea.x, lambda.est, min_sec)
y3.est <- sort(TP(idea.x, myu, min_sec, max_sec))

lines(y1.est, seq(length(y1.est), 1, -1),col="red", lwd = 3)
lines(y2.est, seq(length(y2.est), 1, -1),col="blue", lwd = 3)
lines(y3.est, seq(length(y3.est), 1, -1),col="green", lwd = 3)




###########170826

log_likelihood_for_exp <- function(x,a){
  return(function(par){
    alpha <- par[1]
    C <- par[2]
    sum( log10( alpha*exp(-alpha*(x-a)) ) )
  })
}


log_likelihood_for_Sexp <- function(x,a){
  return(function(par){
    alpha <- par[1]
    beta <- par[2]
    A <- par[3]
    sum( log10( A*exp( -(alpha*(x-a))^beta ) ) )
  })
}



x <- rexp(1000, 0.5)
hist(x)

r <- optim(par = c(0.0001,0.0001), fn = log_likelihood_for_exp(x,0), control = list(fnscale = -1))
# control = list(fnscale = -1) は、optim関数はディフォルトでは最小化をすることになっているので、そうではなく最大化してくれということを伝えているだけ
t <- seq(0,20,length=1001)
plot(t, (1/r$par[1])*exp(-(1/r$par[1])*(t-0)))
plot(t, 0.5*exp(-0.5*(t-0)))

optim(par = c(0.1, 0.1, 0.1), fn = log_likelihood_for_Sexp(x,0), control = list(fnscale = -1))
# control = list(fnscale = -1) は、optim関数はディフォルトでは最小化をすることになっているので、そうではなく最大化してくれということを伝えているだけ

