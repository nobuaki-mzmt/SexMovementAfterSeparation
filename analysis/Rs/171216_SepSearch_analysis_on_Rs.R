## Sex movement after separation
## Formal analysis
## 2017/10/01 N. Mizumoto
## Data of Reticulitermes speratus

#### packages
library(data.table)
library(exactRankTests)
library(CircStats)
library(survival)

#### functions
# standard error
se  <-  function(x)
{
  y  <-  x[!is.na(x)]  #  remove  the  missing  values
  sqrt(var(as.vector(y))/length(y))
}

# angle calcuration
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

# 最尤推定による指数の推定値 for power law
PL_mle <- function(x,xmin){
  return(1+length(x)*(sum(log(x/xmin)))^(-1))
}
# 最尤推定による指数の推定値 for exponential 
EXP_mle <- function(x,xmin){
  return(length(x)*(sum(x-xmin)^(-1)))
}

# 対数尤度関数 for power law
PL_llh <- function(x,xmin,a){
  n <- length(x)
  return(n*log(a-1)+n*(a-1)*log(xmin)-a*sum(log(x)))
}

# 対数尤度関数 for power law
EXP_llh <- function(x,xmin,a){
  n <- length(x)
  return(n*(log(a)+a*xmin)-a*sum(x))
}


#### data
## tandem data (60min observation)
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\tandem\\raw_data")
Folder.path <- getwd()
f.namesplace_tandem <- list.files(Folder.path, pattern=".csv",full.names=T)

## solo tandem (30min before - 10min tandem - 30min after)
setwd("E:\\Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623")
Folder.path <- getwd()
f.namesplace_solo <- list.files(Folder.path, pattern=".txt",full.names=T)

#### setting 
anal_FPS <- 5
cut_time <- 3
##############



#### main

###########################################################
###########################################################
###########################################################
#### 1. observation of tandem running (60min observation)
## porpose: extracting several data on tandem running

ID = ID2 = ID3 = ID4 = sep_time = tan_time = cens = cens2 = sep_male_speed = sep_female_speed = sep_dis <- NULL
fem_speed_dynamics = MF_dis_dynamics <- NULL

for(j in 1:length(f.namesplace_tandem)){
  data.path <- f.namesplace_tandem[j]
  d <- data.frame(as.matrix(fread(data.path, header=T)))
  datepoint <- regexpr("170", data.path)[1]
  id <- (substr(data.path, datepoint, datepoint+11))
  
  ## analyze 60 minutes data (if less than 60 minutes, use as long as possible)
  d[,1] <- d[,1]/30
  if(length(d[,1]) >= 30*60*60+1){
    anal_time <- seq(1, 30*60*60+1, 30/anal_FPS)
  } else {
    anal_time <- seq(1, length(d[,1]), 30/anal_FPS)
  }
  
  ## scaling
  # 満遍なく移動したと想定（max - min = 145mm )
  Xmin <- min(c(d$x0,d$x1)); Xmax <- max(c(d$x0,d$x1)); 
  Ymin <- min(c(d$y0,d$y1)); Ymax <- max(c(d$y0,d$y1)); 
  XRange <- Xmax-Xmin;
  YRange <- Ymax-Ymin;
  d$x0 <- (d$x0-Xmin)*145/XRange
  d$x1 <- (d$x1-Xmin)*145/XRange
  d$y0 <- (d$y0-Ymin)*145/YRange
  d$y1 <- (d$y1-Ymin)*145/YRange
  
  da <- d[anal_time,]
  
  # calcurate each component
  time <- da$position
  Male_X <- da$x1
  Male_Y <- da$y1
  Female_X <- da$x0
  Female_Y <- da$y0
  Male_speed <- c(NA, sqrt ( diff(Male_X)^2 + diff(Male_Y)^2 ) * 5)
  Female_speed <- c(NA, sqrt ( diff(Female_X)^2 + diff(Female_Y)^2 ) * 5)
  
  L <- length(time)
  
  # move direction (0-360 degree)
  Man <- atan( (Male_Y[2:L] - Male_Y[2:L-1])/(Male_X[2:L] - Male_X[2:L-1]) )
  Man[ Male_X[2:L] - Male_X[2:L-1] > 0 & !is.na(Man)] <- Man[ Male_X[2:L] - Male_X[2:L-1] > 0 & !is.na(Man)] + pi
  Man[Man<0 & !is.na(Man)] <- Man[Man<0 & !is.na(Man)] + 2*pi
  Man <- (Man / pi * 180)
  
  Feman <- atan( (Female_Y[2:L] - Female_Y[2:L-1])/(Female_X[2:L] - Female_X[2:L-1]) )
  Feman[ Female_X[2:L] - Female_X[2:L-1] > 0 & !is.na(Feman)] <- Feman[ Female_X[2:L] - Female_X[2:L-1] > 0 & !is.na(Feman)] + pi
  Feman[Feman<0 & !is.na(Feman)] <- Feman[Feman<0 & !is.na(Feman)] + 2*pi
  Feman <- (Feman / pi * 180)
  
  MF_dis <- sqrt( (Male_X-Female_X)^2 + (Male_Y-Female_Y)^2 )
  
  
  ## definition of tandem running
  
  ## def1: 距離7mm以下(approximate the largest size of Rs dealates, icluding the length of antena)が
  ## 3sec以上続いて、初めてタンデムとみなす
  ## (similarlly a pair is considered to get separated when the distance >7mm for >3sec)
  tandem <- MF_dis <= 7
  sep_beg <- 1
  for(i in 2:L){
    if(!tandem[i] & tandem[i-1]){
      sep_beg <- i
    }
    if(!tandem[i-1] & tandem[i]){
      sep_en <- i-1
      if(time[sep_en] - time[sep_beg] <= 3){
        tandem[sep_beg:sep_en] <- TRUE
      }
    }
  }
  
  tan_beg <- 1
  tan_en <- 1
  for(i in 2:L){
    if(tandem[i] & !tandem[i-1]){
      tan_beg <- i
    }
    if(tandem[i-1] & !tandem[i]){
      tan_en <- i-1
      if(time[tan_en] - time[tan_beg] <= 3){
        tandem[tan_en:tan_beg] <- FALSE
      }
    }
  }
  
  ## def2: 近くにいる期間に、進んだ距離が30mm以下ならタンデムとみなさない
  ## (if they don't move we cannot consider that to be tandem "running")
  
  # labeling tandem events
  tandem_label <- rep(0,L)
  count <- 1
  for(i in 1:L){
    if(tandem[i]){
      tandem_label[i] <- count
    } else if(i>1&&tandem[i-1]) {
      count <- count + 1
    }
  }
  for(i in 1:max(tandem_label)){
    if( sum(Female_speed[tandem_label==i]/5, na.rm=T) <= 30){
      tandem[tandem_label==i] <- FALSE
    }
  }
  
  # separate
  separate <- !tandem
  for(i in 1:L){
    if(!separate[i]){
      separate[1:i] <- FALSE
      break;
    }
  }
  # separateイベントにラベル付け
  separate_label <- rep(0,L)
  count <- 1
  for(i in 1:L){
    if(separate[i]){
      separate_label[i] <- count
    } else if(i>1&&separate[i-1]) {
      count <- count + 1
    }
  }
  
  ## def3: sep中2個体の方向の差の最大が45度以下ならsepとみなさない
  ## Although the distance between individuals sometimes become >7mm,
  ## actually they sometimes looks to perform tandem running.
  for(i in 1:max(separate_label)){
    an_dif <- abs(Man[separate_label==i] - Feman[separate_label==i])
    rev_an_dif <- abs(360-an_dif)
    an_dif[ an_dif > rev_an_dif & !is.na(an_dif)] <- rev_an_dif[ an_dif > rev_an_dif & !is.na(an_dif)]
    if( max(an_dif, na.rm=T) < 45 ){
      separate[separate_label==i] <- FALSE
    }
  }
  tandem <- !separate
  
  ## 最ラベルつけ
  separate_label <- rep(0,L)
  count <- 1
  for(i in 1:L){
    if(separate[i]){
      separate_label[i] <- count
    } else if(i>1&&separate[i-1]) {
      count <- count + 1
    }
  }
  tandem_label <- rep(0,L)
  count <- 1
  for(i in 1:L){
    if(tandem[i]){
      tandem_label[i] <- count
    } else if(i>1&&tandem[i-1]) {
      count <- count + 1
    }
  }
  
  ## The first not separation but before encounter
  if(separate[L]){
    separate[separate_label==max(separate_label)] <- FALSE
  }
  ## output
  if(max(separate_label)==0){ # no separation
    sep_time <- c(sep_time, 0)
    tan_time <- c(tan_time, sum(tandem>0) / anal_FPS)
    cens <- c(cens, 0)
    cens2 <- c(cens2, 0)
    #sep_male_speed <- c(sep_male_speed, NA)
    #sep_female_speed <- c(sep_female_speed, NA)
    sep_dis <- c(sep_dis, NA)
    ID <- c(ID, rep(j, 1))
    ID2 <- c(ID2, rep(id, 1))
    ID3 <- c(ID3, rep(j, 1))
    ID4 <- c(ID4, rep(id, 1))
  } else {
    for(i in 1:max(separate_label)){
      sep_time <- c(sep_time, length(MF_dis[separate_label == i]) / anal_FPS )
      #sep_male_speed <- c(sep_male_speed, mean(Male_speed[separate_label == i]) )
      #sep_female_speed <- c(sep_female_speed, mean(Female_speed[separate_label == i]) )
      sep_dis <- c(sep_dis, mean(MF_dis[separate_label == i]) )
      fem_speed_dynamics <- c(fem_speed_dynamics, list( Female_speed[separate_label == i] ) )
      MF_dis_dynamics <- c(MF_dis_dynamics, list( MF_dis[separate_label == i] ) )
    }
    for(i in 1:max(tandem_label)){
      tan_time <- c(tan_time, length(MF_dis[tandem_label == i]) / anal_FPS )
      if(tandem_label[length(tandem_label)] == i){
        cens <- c(cens, 0)
      }else {
        cens <- c(cens, 1)
      }
    }
    for(i in 1:max(separate_label)){
      if(separate_label[length(separate_label)] == i){
        cens2 <- c(cens2, 0)
      }else {
        cens2 <- c(cens2, 1)
      }
    }
    
    ID <- c(ID, rep(j, max(separate_label)))
    ID2 <- c(ID2, rep(id, max(separate_label)))
    ID3 <- c(ID3, rep(j, max(tandem_label)))
    ID4 <- c(ID4, rep(id, max(tandem_label)))
    
  }
  print(j)
}

# analyzed pair: 29 pairs
# separation was observed in 26/29


res1 <- cbind(ID, ID2, sep_time, cens2)
res2 <- cbind(ID3, ID4, tan_time, cens)

# write.table(res1, "E:\\res1.txt")

# observed separation events
sum(sep_time>0)
# 283

sum(sep_time==0)
# 3

# survival curve for tandem continueing time
ge.sf<-survfit(Surv(tan_time,cens)~0)
par(mfrow=c(1,1), pin=c(4,3))
plot(ge.sf, xlab="time (sec)", ylab="probability of tandem running", main="R. speratus", las=1)
Rs_tantime <- tan_time
Rs_cens <- cens

# survival curve for separation time
ge.sf<-survfit(Surv(sep_time[sep_time>0],cens2[sep_time>0])~0)
par(mfrow=c(1,1), pin=c(4,3))
plot(ge.sf, xlab="time (sec)", ylab="probability of tandem running", main="R. speratus", las=1)
Rs_septime <- sep_time
Rs_cens2 <- cens2


# expected distance 
# calcurate the mode of separated distance (during whole of events)
sep_def_dis <- NULL
for(i in 1:length(MF_dis_dynamics)){
  dd <- MF_dis_dynamics[[i]]
  sep_def_dis <- c(sep_def_dis, dd)
}
truehist(sep_def_dis)
y <- density(sep_def_dis)
plot(y)
y$x[y$y == max(y$y)]
# 16.09702


###########################################################
###########################################################
###########################################################
#### 2. observation of searching behavior (before encounter / get separated)
## porpose: comparing the searching patterns

###########################################################
## 2-1 Comparison of movement patterns
## movement patterns = move speed, turning angle, pause time
Colony = Sex = ID = Minutes = Treat = Mean_speed = Mean_angle = Mean_pause <- NULL
mincom <- 1:30
for(j in 1:length(f.namesplace_solo)){
  data.path <- f.namesplace_solo[j]
  d <- read.table(data.path, header=T)
  
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  date <- (substr(data.path, datepoint, datepoint+5))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+13, datepoint+13))
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  L <- length(d[,1])
  stop <- (d$speed <= 0.7*5)
  time2 <- d$time[!stop];  x2 <- d$x[!stop];  y2 <- d$y[!stop];  L2 <- length(x2)
  
  # moving turning angle (stop時間を除くangleを計算)
  # fitting with wrapped Cauchy 
  if(length(x2)==1){ tra_angle <- NA } else{
    kari_angle <- c(NA, angle_cal(x2, y2, L2), NA)
    tra_angle <- rep(NA, L)
    for(i in 2:L2){  tra_angle[d$time==time2[i]] <- kari_angle[i] }
  }
  speed2 <- d$speed[!stop]
  
  ## 1分ごとの平均値
  minutes <- c(0, rep(seq(1,30,1), each = 5*60))
  
  obj <- unique(na.omit(minutes[1:L][!is.na(d$speed)&!stop]))
  mean_speed <- cbind(obj,  tapply(speed2, minutes[1:L][!stop], mean, na.rm=T))
  lack <- cbind(setdiff(1:30, obj))
  mean_speed <- rbind(mean_speed, cbind(lack, rep(NA,length(lack))))
  mean_speed <- mean_speed[order(mean_speed[,1]),]
  
  obj <- unique(na.omit(minutes[!is.na(tra_angle)]))
  mean_angle <- cbind(obj,  tapply(na.omit(abs(tra_angle)), minutes[1:L][!is.na(tra_angle)], mean))
  lack <- cbind(setdiff(1:30, obj))
  mean_angle <- rbind(mean_angle, cbind(lack, rep(NA,length(lack))))
  mean_angle <- mean_angle[order(mean_angle[,1]),]
  
  mean_pause <- tapply(stop[2:L], minutes[2:L], sum, na.rm=T)*0.2
  if(length(mean_pause)<30){mean_pause<-c(mean_pause,rep(NA,30-length(mean_pause)))}
  
  
  Colony <- c(Colony, rep(colony, 30))
  Sex <- c(Sex, rep(sex, 30))
  ID <- c(ID, rep(paste(date,iter,sep="_"), 30))
  Treat <- c(Treat, rep(treat, 30))
  Minutes <- c(Minutes, mincom)
  Mean_speed <- c(Mean_speed, as.vector(mean_speed[,2]))
  Mean_angle <- c(Mean_angle, as.vector(mean_angle[,2]))
  Mean_pause <- c(Mean_pause, as.vector(mean_pause))
  
}  

res <- na.omit(data.frame(colony=Colony, sex=Sex, is=ID, treat=Treat, min=Minutes, mean_speed=Mean_speed, mean_angle=Mean_angle, mean_pause=Mean_pause))
res_before <- res[res$treat=="before",]
res_tandem <- res[res$treat=="tandem",]
res_after <- res[res$treat=="after",]

par(mfrow = c(2,3), pin=c(2.2,2.2))

## plot (for pause time)
res_before_mean <- tapply(res_before$mean_pause, res_before[,c(5,2)], mean, na.rm=T)
res_tandem_mean <- tapply(res_tandem$mean_pause, res_tandem[,c(5,2)], mean, na.rm=T)[1:10,]
res_after_mean <- tapply(res_after$mean_pause, res_after[,c(5,2)], mean)
res_before_se <- tapply(res_before$mean_pause, res_before[,c(5,2)], sd)
res_tandem_se <- tapply(res_tandem$mean_pause, res_tandem[,c(5,2)], sd)[1:10,]
res_after_se <- tapply(res_after$mean_pause, res_after[,c(5,2)], sd)
range <- c(-5,70)
plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=range)
arrows(mincom, res_before_mean-res_before_se, mincom, res_before_mean+res_before_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(mincom, res_before_mean, xlim=c(1,30), ylim=range, cex=1.5, type="o", col=c(2,1),
        lty=c(2,1), pch=20, ylab="pause time (sec)", xlab="time (min)", main="before")

plot(0, ann=F, axes=F, type="n", xlim=c(1,10), ylim=range)
arrows(1:10, res_tandem_mean-res_tandem_se, 1:10, res_tandem_mean+res_tandem_se, code=3,
       angle=90, length=0.05, col=c(rep(2,10),rep(1,10)))
par(new=T)
matplot(1:10, res_tandem_mean, xlim=c(1,10), ylim=range, 
        cex = 1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="pause time (sec)", xlab="time (min)", main="tandem")

plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=range)
arrows(mincom, res_after_mean-res_after_se, mincom, res_after_mean+res_after_se, code=3,
       angle=90, length=0.05, col=c(rep(2,30),rep(1,30)))
par(new=T)
matplot(res_after_mean, xlim=c(1,30), ylim=range, 
        cex=1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="pause time (sec)", xlab="time (min)", main="after")


# ver2
matplot(mincom, res_before_mean, xlim=c(1,30), ylim=range, cex=1.5, type="o", col=c(2,1),
        lty=c(2,1), pch=20, ylab="pause time (sec)", xlab="time (min)", main="before")
par(new=T)
matplot(mincom, res_before_mean-res_before_se, xlim=c(1,30), ylim=range,type="l", col=c(2,1),
        lty=c(2,1), ann=F, axes=F)
par(new=T)
matplot(mincom, res_before_mean+res_before_se, xlim=c(1,30), ylim=range,type="l", col=c(2,1),
        lty=c(2,1), ann=F, axes=F)

matplot(1:10, res_tandem_mean, xlim=c(1,10), ylim=range, cex=1.5, type="o", col=c(2,1),
        lty=c(2,1), pch=20, ylab="pause time (sec)", xlab="time (min)", main="before")
par(new=T)
matplot(1:10, res_tandem_mean-res_tandem_se, xlim=c(1,10), ylim=range,type="l", col=c(2,1),
        lty=c(2,1), ann=F, axes=F)
par(new=T)
matplot(1:10, res_tandem_mean+res_tandem_se, xlim=c(1,10), ylim=range,type="l", col=c(2,1),
        lty=c(2,1), ann=F, axes=F)

matplot(mincom, res_after_mean, xlim=c(1,30), ylim=range, cex=1.5, type="o", col=c(2,1),
        lty=c(2,1), pch=20, ylab="pause time (sec)", xlab="time (min)", main="before")
par(new=T)
matplot(mincom, res_after_mean-res_after_se, xlim=c(1,30), ylim=range,type="l", col=c(2,1),
        lty=c(2,1), ann=F, axes=F)
par(new=T)
matplot(mincom, res_after_mean+res_after_se, xlim=c(1,30), ylim=range,type="l", col=c(2,1),
        lty=c(2,1), ann=F, axes=F)




## plot (for move speed)
res_before_mean <- tapply(res_before$mean_speed, res_before[,c(5,2)], mean, na.rm=T)
res_tandem_mean <- tapply(res_tandem$mean_speed, res_tandem[,c(5,2)], mean, na.rm=T)[1:10,]
res_after_mean <- tapply(res_after$mean_speed, res_after[,c(5,2)], mean)
res_before_se <- tapply(res_before$mean_speed, res_before[,c(5,2)], sd)
res_tandem_se <- tapply(res_tandem$mean_speed, res_tandem[,c(5,2)], sd)[1:10,]
res_after_se <- tapply(res_after$mean_speed, res_after[,c(5,2)], sd)
range <- c(0,25)
plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=range)
arrows(mincom, res_before_mean-res_before_se, mincom, res_before_mean+res_before_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(mincom, res_before_mean, xlim=c(1,30), ylim=range, 
        cex=1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="before")

plot(0, ann=F, axes=F, type="n", xlim=c(1,10), ylim=range)
arrows(1:10, res_tandem_mean-res_tandem_se, 1:10, res_tandem_mean+res_tandem_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(1:10, res_tandem_mean, xlim=c(1,10), ylim=range, 
        cex = 1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="tandem")

plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=range)
arrows(mincom, res_after_mean-res_after_se, mincom, res_after_mean+res_after_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(res_after_mean, xlim=c(1,30), ylim=range, 
        cex=1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="after")


## plot (for turning angle)
res_before_mean <- tapply(res_before$mean_angle, res_before[,c(5,2)], mean, na.rm=T)
res_tandem_mean <- tapply(res_tandem$mean_angle, res_tandem[,c(5,2)], mean, na.rm=T)[1:10,]
res_after_mean <- tapply(res_after$mean_angle, res_after[,c(5,2)], mean)
res_before_se <- tapply(res_before$mean_angle, res_before[,c(5,2)], sd)
res_tandem_se <- tapply(res_tandem$mean_angle, res_tandem[,c(5,2)], sd)[1:10,]
res_after_se <- tapply(res_after$mean_angle, res_after[,c(5,2)], sd)
range <- c(0,1.5)
plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=range)
arrows(mincom, res_before_mean-res_before_se, mincom, res_before_mean+res_before_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(mincom, res_before_mean, xlim=c(1,30), ylim=range, 
        cex=1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="before")

plot(0, ann=F, axes=F, type="n", xlim=c(1,10), ylim=range)
arrows(1:10, res_tandem_mean-res_tandem_se, 1:10, res_tandem_mean+res_tandem_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(1:10, res_tandem_mean, xlim=c(1,10), ylim=range, 
        cex = 1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="tandem")

plot(0, ann=F, axes=F, type="n", xlim=c(1,30), ylim=range)
arrows(mincom, res_after_mean-res_after_se, mincom, res_after_mean+res_after_se, code=3,
       angle=90, length=0.05)
par(new=T)
matplot(res_after_mean, xlim=c(1,30), ylim=range, 
        cex=1.5, type="o", col=c(2,1), lty=c(2,1), pch=20 ,
        ylab="angle (rad)", xlab="time (min)", main="after")



###########################################################
## 2-2 Detailed comparison
## Porpose: modeling their movement (mimicking)
## There are significant differences between sexes both in speed and turning angle in 3 minutes after separation.
## analysis on the first 3 minutes after separation and the last 3 minutes before encounters (control)

## intermittency (move/pause)
## move: CRW (speed, turning angle every 0.2 sec)

##### 2-2-1 parameter for CRW
before <- f.namesplace_solo[regexpr("_before",f.namesplace_solo)>0]
tandem <- f.namesplace_solo[regexpr("_tandem",f.namesplace_solo)>0]
after <- f.namesplace_solo[regexpr("_after",f.namesplace_solo)>0]
file_table <- cbind(before, tandem, after)

#par(mfrow=c(4,4), pin=c(2,1))
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
    
    if(j==1){ anal_time <-  (60*(30-cut_time)*anal_FPS):(60*30*anal_FPS)+1}
    if(j==3){ anal_time <-  (60*0):(60*cut_time*anal_FPS)+1}
    
    time <- d[anal_time, 1]; x <- d[anal_time, 2]; y <- d[anal_time, 3]; L <- length(x)
    dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
    moved <- sum(dis, na.rm=T)
    # pause time
    # 0.2秒の移動距離が0.7以下の場合を止まっているとみなす (= 20% of the mode of the 2nd peak of speed)
    stop <- (dis < 0.7)
    stop[1] <- FALSE
    stop_time <- sum(stop[2:length(stop)])/anal_FPS
    
    time2 <- time[!stop];  x2 <- x[!stop];  y2 <- y[!stop];  L2 <- length(x2)
    
    # moving speed (stop時間を除く平均速度)
    #if(length(x2)==1){ true_speed = NA; } else{
    #  dis2 <- (c(NA,((x2[2:L2]-x2[2:L2-1])^2 + (y2[2:L2]-y2[2:L2-1])^2)^0.5))
    #  speed2 <- dis2*anal_FPS
    #  true_speed <- mean(speed2, na.rm=T)
    #}
    speed2 <- d$speed[!stop]
    true_speed <- mean(speed2, na.rm=T)
    
    # moving turning angle (stop時間を除くangleを計算)
    # fitting with wrapped Cauchy 
    if(length(x2)==1){ tra_angle <- NA } else{
      kari_angle <- c(NA, angle_cal(x2, y2, L2), NA)
      tra_angle <- rep(NA, L2)
      for(i in 1:L2){  tra_angle[time==time2[i]] <- kari_angle[i] }
    }
    if(length(na.omit(tra_angle))>0){
      mle_wrpcauchy <- wrpcauchy.ml(na.omit(tra_angle), 0, 0, acc=1e-015)
      #theta <- seq(-pi, pi, length=1000)
      #probabl <- dwrpcauchy(theta, as.numeric(mle_wrpcauchy[1]), as.numeric(mle_wrpcauchy[2]))
      #Range <- ceiling(max(probabl,na.rm=T))
      #truehist(tra_angle, prob=T, breaks=seq(-pi, pi, length=100),
      #         ylim=c(0,Range), xlim=c(-pi,pi), col=cols[j],
      #         xlab="angle", ylab="density", main=paste(id,treat), axes=F)
      #points(theta, probabl, type="l")  
      #axis(1, at=seq(-pi,pi,length=3), label=c("-pi", "0", "pi"))
      #axis(2, las=1)
      mu <- as.numeric(mle_wrpcauchy[1])
      rho <- as.numeric(mle_wrpcauchy[2])
    } else{ mu = rho = NA; }
    
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

d <- data.frame(id=ID, colony=Colony, sex=Sex, treat=Treat, moved=Moved, stop_time=Stop_time, mean_speed=True_speed, Rho=Rho, Mu=Mu)

### plot data and statistical analysis
par(pin=c(2,2), mfrow=c(2,2))
treat <- c(1,0,1,0)
sex <- c(0,0,1,1)


### True speed
Mean <- (tapply(d[,7],d[,3:4],mean,na.rm=T))[2:1,2:1]
#      treat
#sex   before    after
#  1 11.70606 11.570073
#  0 12.52358  5.622711
SE <- (tapply(d[,7],d[,3:4],se))[2:1,2:1]
#       treat
#sex    before     after
#  1  0.8081606 0.6861828
#  0  1.2375496 0.6757735
Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,20), beside=T, las=1,
               ylab="moving speed (mm/sec)")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topright", c("male","female"), fill=c(1,0))　
wilcox.exact(mean_speed ~ sex, data=d[d$treat=="after",], paired=F)
#W = 31, p-value = 3.844e-06
wilcox.exact(mean_speed ~ sex, data=d[d$treat=="before",], paired=F)
#W = 198, p-value = 0.4252
wilcox.exact(mean_speed ~ treat, data=d[d$sex==0,], paired=T)
#V = 7, p-value = 0.000145
wilcox.exact(mean_speed ~ treat, data=d[d$sex==1,], paired=T)
#V = 93, p-value = 0.953

### Rho (wrapped caucy distribution)
Mean <- (tapply(d[,8],d[,3:4],mean, na.rm=T))[2:1,2:1]
#       treat
#sex    before     after
#  1 0.8467515 0.7939391
#  0 0.8659342 0.7576250
SE <- (tapply(d[,8],d[,3:4],se))[2:1,2:1]
#         treat
#sex      before      after
#  1 0.009184568 0.01192549
#  0 0.007299856 0.04768086
Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,1), beside=T, las=1,
               ylab="ρ")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topright", c("male","female"), fill=c(1,0))　
wilcox.exact(Rho ~ sex, data=d[d$treat=="after",], paired=F)
#W = 122, p-value = 0.4929
wilcox.exact(Rho ~ sex, data=d[d$treat=="before",], paired=F)
#W = 189, p-value = 0.23
wilcox.exact(Rho ~ treat, data=d[d$sex==0 & d$id!="170526_3" & d$id!=" 170601_1" & d$id!="170607_2" & d$id!="170608_2",], paired=T)
#V = 17, p-value = 0.02454 (n -= 4)
wilcox.exact(Rho ~ treat, data=d[d$sex==1,], paired=T)
#V = 17, p-value = 0.0007896



###### 2-2-2 parameter for Intermittency

### 2-2-2-0. 止まっている判定の閾値決め
Dis <- NULL
for(j in 1:length(f.namesplace_solo)){
  data.path <- f.namesplace_solo[j]
  d <- read.table(data.path, header=T)
  datepoint <- regexpr("170", data.path)[1]+7
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  #if(treat!="before" || sex==1){ next; }
  time <- d[, 1];  x <- d[, 2];  y <- d[, 3];  L <- length(x)
  dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
  Dis <- c(Dis, dis)
}
MaxBin <- ceiling(max(Dis, na.rm=T)*10)/10
par(mfrow=c(1,1), pin=c(3,3))
h <- truehist(Dis)
h <- truehist(Dis, xlim=c(0,2), breaks=seq(0,MaxBin,0.01))
y <- density(Dis[Dis>2], na.rm=T)
y$x[y$y == max(y$y)]

stop_thresh <- 0.7

### 2-2-2-1 Fitting to Experimental data
## generating step data
Stop_sec = Move_sec <- NULL
Treat <- "after"
Sex <- 0
for(j in 1:length(f.namesplace_solo)){
  data.path <- f.namesplace_solo[j]
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
  stop <- (dis <= stop_thresh)
  stop[1] <- FALSE
  # 3min時のmove/pauseが終了するまでのデータを見る
  min3 <- 60*cut_time*anal_FPS
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
  stop_sec <- as.vector(tapply(stop[1:L], stop_label[1:L], sum))*(1/anal_FPS)
  Stop_sec <- c(Stop_sec, stop_sec)
  move_sec <- as.vector(tapply(!stop[1:L], move_label[1:L], sum))*(1/anal_FPS)
  Move_sec <- c(Move_sec, move_sec)
}

par(mfrow=c(2,1), pin=c(4,3)*0.8)

## fitting pause
min_sec <- 0.2
Stop_lab <- sort(unique(Stop_sec[Stop_sec>=min_sec]))
Stop_count <- as.vector(table(Stop_sec[Stop_sec>=min_sec]))
Stop_ratio <- Stop_count / sum(Stop_count)
Stop_cum <- rep(0, length(Stop_ratio))
for(i in 1:length(Stop_ratio)){
  Stop_cum[i] <- sum(Stop_ratio[i:length(Stop_ratio)])
}

idea.x <- seq(0,0.999,0.001)
myu_P <- PL_mle(Stop_sec[Stop_sec>=min_sec], min_sec)
myu_P # 1.739645
y1.est_P <- PL(idea.x, myu_P, min_sec)

lambda.est_P <- EXP_mle(Stop_sec[Stop_sec>=min_sec], min_sec)
lambda.est_P # 0.1217179
y2.est_P <- EXP(idea.x, lambda.est_P, min_sec)

plot(Stop_lab,Stop_cum, log="xy", main=paste("pause", Treat,Sex), axes=F,
     ylab="inverse cumulative frequency", xlab="pause time", ylim=c(0.001,1))
axis(1, c(0.2,0.5,1,2,5,10,20,50,100,200,500))
axis(2, c(0.001,0.01,0.1,1))
points(y1.est_P, rev(idea.x), col="red", lwd = 1.5, type="l")
points(y2.est_P, rev(idea.x), col="blue", lwd = 1.5, type="l")


## fitting
K1 <- 1
K2 <- 1
AIC1.exp <- -2*PL_llh(Stop_sec[Stop_sec>=min_sec], min_sec, myu_P)+2*K1
AIC2.exp <- -2*EXP_llh(Stop_sec[Stop_sec>=min_sec], min_sec, lambda.est_P)+2*K2
AIC <- c(AIC1.exp, AIC2.exp)
delta <- AIC - min(AIC)
w1.exp <- exp(-delta[1]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w2.exp <- exp(-delta[2]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w1.exp # 1
w2.exp # 0
## fit to PL

## fitting move
Move_lab <- sort(unique(Move_sec[Move_sec>=min_sec]))
Move_count <- as.vector(table(Move_sec[Move_sec>=min_sec]))
Move_ratio <- Move_count / sum(Move_count)
Move_cum <- rep(0, length(Move_ratio))
for(i in 1:length(Move_ratio)){
  Move_cum[i] <- sum(Move_ratio[i:length(Move_ratio)])
}

idea.x <- seq(0,0.99999,0.00001)
myu <- PL_mle(Move_sec[Move_sec>=min_sec], min_sec)
myu # 2.19908
y1.est <- PL(idea.x, myu, min_sec)

lambda.est <- EXP_mle(Move_sec[Move_sec>=min_sec], min_sec)
lambda.est # 1.261352
y2.est <- EXP(idea.x, lambda.est, min_sec)

plot(Move_lab,Move_cum, log="xy", main=paste("move", Treat,Sex), axes=F,
     ylab="inverse cumulative frequency", xlab="move time", ylim=c(0.001,1))
axis(1, c(0.2,0.5,1,2,5,10,20,50,100,200,500))
axis(2, c(0.001,0.01,0.1,1))
points(y1.est, rev(idea.x), col="red", lwd = 1.5, type="l")
points(y2.est, rev(idea.x), col="blue", lwd = 1.5, type="l")

K1 <- 1
K2 <- 1
AIC1.exp <- -2*PL_llh(Move_sec[Move_sec>=min_sec], min_sec, myu)+2*K1
AIC2.exp <- -2*EXP_llh(Move_sec[Move_sec>=min_sec], min_sec, lambda.est)+2*K2
AIC <- c(AIC1.exp, AIC2.exp)
delta <- AIC - min(AIC)
w1.exp <- exp(-delta[1]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w2.exp <- exp(-delta[2]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w1.exp # 1
w2.exp # 4.017707e-158
## fit to PL


### 2-2-2-1 Making randomized data
## コントロールとして、3分間分のpauseとmove
Stop <- NULL
Treat <- "after"
Sex <- 0
for(j in 1:length(f.namesplace_solo)){
  data.path <- f.namesplace_solo[j]
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
  stop <- (dis <= stop_thresh)
  stop[1] <- FALSE
  
  min3 <- 60*3*anal_FPS
  Stop <- c(Stop,stop[1:min3])
}






Stop <- na.omit(Stop)
sum(Stop) / length(Stop)
# pause: 0.852037
1 - sum(Stop) / length(Stop)
# move: 0.147963

P <- sum( PL(runif(100,0,1), 2.0, 0.2) )
M <- sum( PL(runif(100,0,1), 3.0, 0.2) )
P/(M+P)

###########################################################
## 2-3 MSD calcuration
## Porpose: pausing of separated female contribute low diffusivity?

#############
## 2-3-1 MSD (time)
analyze_time <- 3 # min
an_sec <- analyze_time*60
par(pin=c(4,4), mfrow=c(1,1))
plot(0, type="n", axes=F,  ylim=c(-3, 4.5), xlim=c(log10(0.2),log10(an_sec)),
     xlab="tau (sec)", ylab="MSD")
axis(1, at=0:2, labels=c(1, 10, 100))
axis(2, at=seq(-2,4,2), labels=c("10^-2", "1", "10^2", "10^4"), las=1)

Tau <- seq(0.2, an_sec, 0.2)
#Tau <- c(1:9, 1:9*10, 1:(an_sec/100)*100)
if(round(an_sec/100) - an_sec/100 != 0){ Tau <- c(Tau, an_sec) }

lt <- length(Tau)
MSD_mean <- rep(0,lt)
for(k in 1:length(f.namesplace_solo)){
  data.path <- f.namesplace_solo[k]
  d <- read.table(data.path, header=T)
  
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  
  if(treat == "before" || treat=="tandem"){next;}
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
## 2-3-2 MSD (dis)

analyze_time <- 3 # min
an_sec <- analyze_time*60
par(pin=c(4,4), mfrow=c(1,1))
plot(0, type="n", axes=F,   xlim=c(0,3), ylim=c(-1.5, 4.5),
     xlab="tau (mm)", ylab="MSD")
axis(1, at=0:3, labels=c(1, 10, 100, 1000))
axis(2, at=seq(0,4,2), labels=c("1", "10^2", "10^4"), las=1)

Tau <- seq(1, 1000, 1)
lt <- length(Tau)
MSD_mean <- rep(0,lt)
for(k in 1:length(f.namesplace_solo)){
  data.path <- f.namesplace_solo[k]
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



