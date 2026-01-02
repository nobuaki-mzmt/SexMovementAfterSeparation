## Sex movement after separation

## 解析ノート for fitting 2017/9/23

## pause / move の時間をそれぞれfittingして、指数を出す

######################

## packages
library(data.table)
library(exactRankTests)
library(CircStats)

## functions

# power law distribution のデータを生成. rが一様乱数. aが1<a<3の指数. xminが分布の最小値
PL <- function(r,a,xmin){
  return(xmin*(1-r)^(-1/(a-1)))
}
# exponential distribution のデータを生成. rが一様乱数. aが指数. xminが分布の最小値
EXP <- function(r,a,xmin){
  return(xmin-1/a*log(1-r))
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



###################

## 0. 止まっている判定の閾値決め

par(mfrow=c(1,1), pin=c(4,3))
Dis <- NULL
for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  d <- read.table(data.path, header=T)
  datepoint <- regexpr("170", data.path)[1]+7
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  if(treat!="after" || sex==1){ next; }
  
  anal_time <-  (60*0*anal_FPS):(60*3*anal_FPS)+1
  time <- d[anal_time, 1]
  x <- d[anal_time, 2]
  y <- d[anal_time, 3]
  
  L <- length(x)
  dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
  Dis <- c(Dis, dis)
}
max(Dis, na.rm=T)
par(mfrow=c(1,2), pin=c(3,3))
truehist(Dis, breaks=seq(0,3.4,length=341), ylim=c(0,4))
truehist(Dis, breaks=seq(0,3.4,length=341), xlim=c(0,1))
truehist(Dis, breaks=seq(0,3.4,length=341), xlim=c(0,0.5))

stop_thresh <- 0.5

## 1. Experimental data

setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

Stop_sec = Move_sec <- NULL
Treat <- "after"
Sex <- 0
for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  d <- read.table(data.path, header=T)
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  if(treat!=Treat || sex!=Sex){ next; }
  
  time <- d[,1]; x <- d[,2]; y <- d[,3]
  L <- length(x)
  dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
  stop <- (dis < stop_thresh)
  stop[1] <- FALSE

  # 3min時のmove/pauseが終了するまでのデータを見る
  min3 <- 60*3*anal_FPS
  mem <- stop[min3]
  an_mem <- mem
  count <- 0
  while(mem == an_mem){
    count <- count +1
    an_mem <- stop[min3+count]
  }
  L <- min3+count
  
  count <- 0
  in_stop <- 0
  stop_label <- rep(NA,L)
  move_label <- rep(NA,L)
  for(i in 2:L){
    if(stop[i]){
      if(in_stop==0){
        count <- count + 1
      }
      in_stop <- 1
      stop_label[i] <- count
    } else {
      in_stop <- 0
      move_label[i] <- count
    }
  }
  
  ## stopで3分に達した場合の処理が必要。
  stop_sec <- as.vector(tapply(stop[1:L], stop_label[1:L], sum))*0.2
  Stop_sec <- c(Stop_sec, stop_sec)
  
  move_sec <- as.vector(tapply(!stop[1:L], move_label[1:L], sum))*0.2
  Move_sec <- c(Move_sec, move_sec)
}

par(mfrow=c(1,2))

## pause
min_sec <- 0.2
Stop_lab <- sort(unique(Stop_sec))
Stop_count <- as.vector(table(Stop_sec))
Stop_ratio <- Stop_count / sum(Stop_count)
Stop_cum <- rep(0, length(Stop_ratio))
for(i in 1:length(Stop_ratio)){
  Stop_cum[i] <- sum(Stop_ratio[i:length(Stop_ratio)])
}

myu <- PL_mle(Stop_sec, min_sec)
# 1.816256
idea.x <- seq(0,0.99999,0.00001)
y1.est <- PL(idea.x, myu, min_sec)
result1 <- ks.test(Stop_sec, y1.est)

lambda.est <- EXP_mle(Stop_sec, min_sec)
idea.x <- seq(0,0.99999,0.00001)
y2.est <- EXP(idea.x, lambda.est, min_sec)
result2 <- ks.test(Stop_sec, y2.est)

idea.x <- seq(0,0.99999,0.00001)
y1.est <- PL(idea.x, myu, min_sec)
y2.est <- EXP(idea.x, lambda.est, min_sec)

plot(Stop_lab,Stop_cum, log="xy", main=paste("pause", Treat,Sex))
lines(y1.est, rev(idea.x), col="red", lwd = 1.5)
lines(y2.est, rev(idea.x), col="blue", lwd = 1.5)


## move
min_sec <- 0.2
Move_lab <- sort(unique(Move_sec))
Move_count <- as.vector(table(Move_sec))
Move_ratio <- Move_count / sum(Move_count)
Move_cum <- rep(0, length(Move_ratio))
for(i in 1:length(Move_ratio)){
  Move_cum[i] <- sum(Move_ratio[i:length(Move_ratio)])
}

myu <- PL_mle(Move_sec, min_sec)
# 2.258725
idea.x <- seq(0,0.99999,0.00001)
y1.est <- PL(idea.x, myu, min_sec)
result1 <- ks.test(Move_sec, y1.est)

lambda.est <- EXP_mle(Move_sec, min_sec)
idea.x <- seq(0,0.99999,0.00001)
y2.est <- EXP(idea.x, lambda.est, min_sec)
result2 <- ks.test(Move_sec, y2.est)

idea.x <- seq(0,0.99999,0.00001)
y1.est <- PL(idea.x, myu, min_sec)
y2.est <- EXP(idea.x, lambda.est, min_sec)

plot(Move_lab,Move_cum, log="xy", main=paste("move",Treat,Sex))
lines(y1.est, rev(idea.x), col="red", lwd = 1.5)
lines(y2.est, rev(idea.x), col="blue", lwd = 1.5)



## 2. surrogate data

setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

randomize_replace <- TRUE
sur_size <- 1000000

stop_thresh <- 0.5
Stop_sec <- NULL
Treat <- "after"
Sex <- 0
Stop <- NULL
for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  d <- read.table(data.path, header=T)
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  if(treat!=Treat || sex!=Sex){ next; }
  
  time <- d[,1]
  x <- d[,2]
  y <- d[,3]
  
  L <- length(x)
  dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
  stop <- (dis < stop_thresh)
  stop[1] <- FALSE
  
  min3 <- 60*3*anal_FPS
  mem <- stop[min3]
  an_mem <- mem
  count <- 0
  while(mem == an_mem){
    count <- count +1
    an_mem <- stop[min3+count]
  }
  L <- min3+count
  Stop <- c(Stop,stop[1:L])
}

sur_stop <- sample(Stop, sur_size, replace=randomize_replace)

# stopの長さと回数を計測
count <- 0
in_stop <- 0
stop_label <- rep(NA,length(sur_stop))
move_label <- rep(NA,length(sur_stop))
for(i in 1:length(Stop)){
  if(sur_stop[i]){
    if(in_stop==0){
      count <- count + 1
    }
    in_stop <- 1
    stop_label[i] <- count
  }else {
    in_stop <- 0
    move_label[i] <- count
  }
}

Stop_sec <- as.vector(tapply(sur_stop, stop_label, sum))*0.2
Move_sec <- as.vector(tapply(!sur_stop, move_label, sum))*0.2

min_sec <- 0.2
Stop_lab <- sort(unique(Stop_sec))
Stop_count <- as.vector(table(Stop_sec))
Stop_ratio <- Stop_count / sum(Stop_count)
Stop_cum <- rep(0, length(Stop_ratio))
for(i in 1:length(Stop_ratio)){
  Stop_cum[i] <- sum(Stop_ratio[i:length(Stop_ratio)])
}

myu <- PL_mle(Stop_sec, min_sec)
idea.x <- seq(0,0.99999,0.00001)
y1.est <- PL(idea.x, myu, min_sec)
result1 <- ks.test(Stop_sec, y1.est)

lambda.est <- EXP_mle(Stop_sec, min_sec)
# 1.088313
idea.x <- seq(0,0.99999,0.00001)
y2.est <- EXP(idea.x, lambda.est, min_sec)
result2 <- ks.test(Stop_sec, y2.est)

idea.x <- seq(0,0.99999,0.00001)
y1.est <- PL(idea.x, myu, min_sec)
y2.est <- EXP(idea.x, lambda.est, min_sec)

plot(Stop_lab,Stop_cum, log="xy", main=paste("pause", Treat,Sex))
lines(y1.est, rev(idea.x), col="red", lwd = 1.5)
lines(y2.est, rev(idea.x), col="blue", lwd = 1.5)

## move
min_sec <- 0.2
Move_lab <- sort(unique(Move_sec))
Move_count <- as.vector(table(Move_sec))
Move_ratio <- Move_count / sum(Move_count)
Move_cum <- rep(0, length(Move_ratio))
for(i in 1:length(Move_ratio)){
  Move_cum[i] <- sum(Move_ratio[i:length(Move_ratio)])
}

myu <- PL_mle(Move_sec, min_sec)
idea.x <- seq(0,0.99999,0.00001)
y1.est <- PL(idea.x, myu, min_sec)
result1 <- ks.test(Move_sec, y1.est)

lambda.est <- EXP_mle(Move_sec, min_sec)
idea.x <- seq(0,0.99999,0.00001)
y2.est <- EXP(idea.x, lambda.est, min_sec)
result2 <- ks.test(Move_sec, y2.est)

idea.x <- seq(0,0.99999,0.00001)
y1.est <- PL(idea.x, myu, min_sec)
y2.est <- EXP(idea.x, lambda.est, min_sec)

plot(Move_lab,Move_cum, log="xy", main=paste("move",Treat,Sex))
lines(y1.est, rev(idea.x), col="red", lwd = 1.5)
lines(y2.est, rev(idea.x), col="blue", lwd = 1.5)

# 0.2, 0.4, 0.6, 0.8, 1.0
# 0.820866142, 0.145669291, 0.025918635, 0.005905512, 0.001640420



## MSD calcuration

pause_myu <- 1.816256
move_myu <- 2.258725

N <- 10000

x <- rep(0, N)
y <- rep(0, N)
time <- rep(0,N)

move <- runif(1,0,1) <= 0.5
for(i in 2:N){
  if(move){
    angle <- 2*pi*runif(1,0,1)
    L <- PL(runif(1,0,1),move_myu,0.2)
    time[i] <- time[i-1] + L
    x[i] <- x[i-1] + cos(angle) * L
    y[i] <- y[i-1] + sin(angle) * L
    move <- !move
  } else {
    L <- PL(runif(1,0,1),pause_myu,0.2)
    time[i] <- time[i-1] + L
    x[i] <- x[i-1]
    y[i] <- y[i-1]
    move <- !move
  }
}

plot(x,y, type="l")
plot(log10(time), log10((x^2+y^2)), log="xy")
lm(log10((x[2:N]^2+y[2:N]^2)) ~ log10(time[2:N]))
abline(0,1)
