## 170911 MSD analysis

#############
## MSD (time)
#############

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

analyze_time <- 3 # min
an_sec <- analyze_time*60
par(pin=c(4,4), mfrow=c(1,1))
plot(0, type="n", axes=F, ann=F, ylim=c(-3, 4.5), xlim=c(log10(0.2),log10(an_sec)))
axis(1, at=0:3, labels=c(1, 10, 100, 1000))
axis(2, at=seq(-2,4,2), labels=c("10^-2", "1", "10^2", "10^4"), las=1)
box()

Tau <- seq(0.2, an_sec, 0.2)
#Tau <- c(1:9, 1:9*10, 1:(an_sec/100)*100)
if(round(an_sec/100) - an_sec/100 != 0){ Tau <- c(Tau, an_sec) }

lt <- length(Tau)
MSD_mean <- rep(0,lt)
for(k in 1:length(f.namesplace)){
  data.path <- f.namesplace[k]
  d <- read.table(data.path, header=T)
  
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  if(treat == "after" || treat=="tandem"){next;}
  if(sex == 1){next;}
  
  x <- d[1:(analyze_time*60*anal_FPS+1),2]
  y <- d[1:(analyze_time*60*anal_FPS+1),3]

  beginx <- x[1]
  beginy <- y[1]
  
  Msd <- (x-beginx)*(x-beginx) + (y-beginy)*(y-beginy)

  MSD_mean <- MSD_mean + Msd
  points(log10(Tau), log10(Msd), col="grey", type="o", lwd=2)
  a <- min(log10(Msd))
  #abline(a=a, b=1, lty=2)
  #abline(a=a, b=2, lty=2)
  r <- lm(log10(Msd)[Msd>0] ~ log10(Tau)[Msd>0])
  r
}

points(log10(Tau), log10(MSD_mean/37), col=2, type="o", lwd=2)

a <- (log10(MSD_mean/37))[Tau==1]
abline(a=a, b=1, lty=2)
abline(a=a, b=2, lty=2)

lm(log10(MSD_mean/37)[Msd>0]~log10(Tau)[Msd>0])



#############
## MSD (dis)
#############

## data
# 170623に5FPSに変換したもの（scaling 済み）を使用
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

analyze_time <- 3 # min
an_sec <- analyze_time*60
par(pin=c(4,4), mfrow=c(1,1))
plot(0, type="n", axes=F, ann=F, xlim=c(0,3), ylim=c(-3, 4.5))
axis(1, at=0:3, labels=c(1, 10, 100, 1000))
axis(2, at=seq(-2,4,2), labels=c("10^-2", "1", "10^2", "10^4"), las=1)
box()

Tau <- seq(1, 1000, 1)
lt <- length(Tau)
MSD_mean <- rep(0,lt)
for(k in 1:length(f.namesplace)){
  data.path <- f.namesplace[k]
  d <- read.table(data.path, header=T)
  
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  if(treat == "after" || treat=="tandem"){next;}
  if(sex == 1){next;}
  
  x <- d[1:(analyze_time*60*anal_FPS+1),2]
  y <- d[1:(analyze_time*60*anal_FPS+1),3]
  beginx <- x[1]
  beginy <- y[1]
  dis <- c(0,sqrt(diff(x)^2 + diff(y)^2))
  Msd <- rep(0,lt)
  x2 = y2 <- rep(0,lt)
  count <- 1
  for(i in 1:length(x)){
    while(sum(dis[1:i]) > Tau[count]){
      ratio <- ( (Tau[count]-sum(dis[1:(i-1)])) / (dis[i]) )
      ax <- x[i-1] + (x[i]-x[i-1]) * ratio
      ay <- y[i-1] + (y[i]-y[i-1]) * ratio
      Msd[count] <- (ax-beginx)^2 + (ay-beginy)^2
      count <- count +1
      if(sum(dis[1:i])>max(Tau)){break}
    }
    if(sum(dis[1:i])>max(Tau)){break}
  }
  points(log10(Tau), log10(Msd), type="o", col="grey")
  MSD_mean <- MSD_mean + Msd
}

points(log10(Tau), log10(MSD_mean/37), col=2, type="o", lwd=2)

a <- (log10(MSD_mean/37))[Tau==1]
abline(a=a, b=1, lty=2)
abline(a=a, b=2, lty=2)

lm(log10(MSD_mean/37)[Msd>0]~log10(Tau)[Msd>0])

  
  
