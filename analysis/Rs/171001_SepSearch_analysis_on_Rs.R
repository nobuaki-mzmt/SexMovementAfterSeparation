## Sex movement after separation
## Formal analysis
## 2017/10/01 N. Mizumoto
## Data of Reticulitermes speratus

#### packages
library(data.table)
library(exactRankTests)
library(CircStats)

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

ID = sep_time = tan_time = cens = sep_male_speed = sep_female_speed = sep_dis <- NULL
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
    #sep_male_speed <- c(sep_male_speed, NA)
    #sep_female_speed <- c(sep_female_speed, NA)
    sep_dis <- c(sep_dis, NA)
    ID <- c(ID, rep(j, 1))
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
    ID <- c(ID, rep(j, max(separate_label)))
  }
  print(j)
}

# analyzed pair: 28 pairs
# separation was observed in 27/28

# observed separation events
sum(sep_time>0)
# 204

# survival curve for tandem continueing time
ge.sf<-survfit(Surv(tan_time,cens)~0)
par(mfrow=c(1,1), pin=c(4,3))
plot(ge.sf, xlab="time (sec)", ylab="probability of tandem running", main="R. speratus", las=1)

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
# 16.76747


###########################################################
###########################################################
###########################################################
#### 2. observation of searching behavior (before encounter / get separated)
## porpose: comparing the searching patterns

###########################################################
## 2-1 Comparison of basic characteristics (speed, turning angle)
Colony = Sex = ID = Minutes = Treat = Mean_speed = Mean_angle <- NULL
mincom <- seq(1,30,1)

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
  Mean_speed <- c(Mean_speed, as.vector(mean_speed))
  Mean_angle <- c(Mean_angle, as.vector(mean_angle))
}  

res <- na.omit(data.frame(colony=Colony, sex=Sex, is=ID, treat=Treat, min=Minutes, mean_speed=Mean_speed, mean_angle=Mean_angle))

## plot (for speed)
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

## plot (for turning angle)
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
# U-test with Bonferroni correction
# speed-after
P <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_after[res_after$min == i,]
  P[i] <- wilcox.exact(mean_speed ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*30,4)
# [1]  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0001  0.0000  0.0025  0.0006  0.0001  0.0161  0.0782
# [15]  0.0053  0.0001  0.0053  0.0004  0.0034  0.0034  0.0123  0.2466  0.7273  0.2716  0.0874  0.0181  0.6940 20.5205
# [29]  0.1445  1.4596

# speed-before
P <- rep(NULL, 30)
P2 <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_before[res_before$min == i,]
  P[i] <- wilcox.exact(mean_speed ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*30,4)
# [1] 13.2896 10.2804 23.2644 21.2215  4.0128  8.5324 28.2057 27.4901 14.9608  5.3536 12.2359  2.4157  2.4157  1.2696
# [15]  5.9721 17.9706  7.3583 12.7565 22.5769 13.8348  3.7780 18.6031  7.7367 25.3592 17.3477 16.7349 19.2449 17.9706
# [29] 25.3592 26.7768

# angle-after
P <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_after[res_after$min == i,]
  P[i] <- wilcox.exact(mean_angle ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*30,4)
# [1]  0.0000  0.0000  0.0025  0.0976  0.1210  0.0976  0.4729  1.2696  0.3603  6.6400  2.4157  1.7065 13.2896 17.3477
# [15]  5.3536  6.2998  4.5149  7.3583 11.7278  5.6568  3.7780 28.9228 18.6031 23.2644 28.9228 11.3847 21.9185  0.5809
# [29] 25.7337 19.7953

# angle-before
P <- rep(NULL, 30)
P2 <- rep(NULL, 30)
for(i in 1:30){
  analysis <- res_before[res_before$min == i,]
  P[i] <- wilcox.exact(mean_angle ~ sex, data=analysis, paired=F)[3]
  P2[i] <- var.test(mean_angle ~ sex, data=analysis, paired=F)[3]
}
round(as.numeric(P)*30,4)
# [1]  4.2583  6.2998 13.8348  4.7829  0.9293  0.5166  5.9721  4.7829  7.7367  3.5539  2.1082  0.0699  1.0885  0.0782
# [15]  0.8573  6.9928  1.8328  2.1082 15.5412  5.3536  2.1082 10.2804  2.1082  3.7780  3.1362  5.0623  9.8238  6.2998
# [29]  8.5324 19.8956


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
    if(length(x2)==1){ true_speed = NA; } else{
    dis2 <- (c(NA,((x2[2:L2]-x2[2:L2-1])^2 + (y2[2:L2]-y2[2:L2-1])^2)^0.5))
    speed2 <- dis2*anal_FPS
    true_speed <- mean(speed2, na.rm=T)
    }
    
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

### moved distance
Mean <- (tapply(d[,5],d[,3:4],mean))[2:1,2:1]
#      treat
#sex   before     after
#  1 2103.157 2596.4913
#  0 2088.458  297.1482
SE <- (tapply(d[,5],d[,3:4],se))[2:1,2:1]
Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,3000*7), beside=T, las=1,
               ylab="moved distance (mm)")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topleft", c("male","female"), fill=c(1,0))　
wilcox.exact(moved ~ sex, data=d[d$treat=="after",], paired=F)
# W = 0, p-value = 1.132e-10
wilcox.exact(moved ~ sex, data=d[d$treat=="before",], paired=F)
# W = 175, p-value = 0.9163
wilcox.exact(moved ~ treat, data=d[d$sex==0,], paired=T)
# V = 3, p-value = 3.815e-05
wilcox.exact(moved ~ treat, data=d[d$sex==1,], paired=T)
# V = 153, p-value = 0.01808

### stop time
Mean <- (tapply(d[,6],d[,3:4],mean))[2:1,2:1]
#      treat
#sex   before     after
#  1 44.44211  14.76842
#  0 46.43333 153.47778
SE <- (tapply(d[,6],d[,3:4],se))[2:1,2:1]
Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,180*7), beside=T, las=1,
               ylab="stop time (sec)")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topleft", c("male","female"), fill=c(1,0))　
wilcox.exact(stop_time ~ sex, data=d[d$treat=="after",], paired=F)
#W = 340, p-value = 4.527e-10
wilcox.exact(stop_time ~ sex, data=d[d$treat=="before",], paired=F)
#W = 129, p-value = 0.2073
wilcox.exact(stop_time ~ treat, data=d[d$sex==0,], paired=T)
#V = 152, p-value = 3.052e-05
wilcox.exact(stop_time ~ treat, data=d[d$sex==1,], paired=T)
#V = 7, p-value = 7.248e-05

### True speed
Mean <- (tapply(d[,7],d[,3:4],mean,na.rm=T))[2:1,2:1]
#      treat
#sex   before    after
#  1 15.01067 15.48564
#  0 15.28351 10.52885
SE <- (tapply(d[,7],d[,3:4],se))[2:1,2:1]
Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,20), beside=T, las=1,
               ylab="moving speed (mm/sec)")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topright", c("male","female"), fill=c(1,0))　
wilcox.exact(mean_speed ~ sex, data=d[d$treat=="after",], paired=F)
#W = 61, p-value = 0.002004
wilcox.exact(mean_speed ~ sex, data=d[d$treat=="before",], paired=F)
#W = 145, p-value = 0.8318
wilcox.exact(mean_speed ~ treat, data=d[d$sex==0,], paired=T)
#V = 18, p-value = 0.007629
wilcox.exact(mean_speed ~ treat, data=d[d$sex==1,], paired=T)
#V = 112, p-value = 0.5153

### Rho (wrapped caucy distribution)
Mean <- (tapply(d[,8],d[,3:4],mean, na.rm=T))[2:1,2:1]
#       treat
#sex    before     after
#  1 0.8467515 0.7939391
#  0 0.8659342 0.7576250
SE <- (tapply(d[,8],d[,3:4],se))[2:1,2:1]
Fig <- barplot(Mean, col=c(1,0), pch=19, ylim=c(0,1), beside=T, las=1,
               ylab="ρ")
arrows(Fig , Mean, Fig , Mean + SE, angle=90,length=0.1)
legend("topright", c("male","female"), fill=c(1,0))　
wilcox.exact(Rho ~ sex, data=d[d$treat=="after",], paired=F)
#W = 122, p-value = 0.4929
wilcox.exact(Rho ~ sex, data=d[d$treat=="before",], paired=F)
#W = 189, p-value = 0.23
wilcox.exact(Rho ~ treat, data=d[d$sex==0 & d$id!="170526_3" & d$id!=" 170601_1" & d$id!="170607_2" & d$id!="170608_2",], paired=T)
#V = 17, p-value = 0.02454
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
Stop_lab <- sort(unique(Stop_sec))
Stop_count <- as.vector(table(Stop_sec))
Stop_ratio <- Stop_count / sum(Stop_count)
Stop_cum <- rep(0, length(Stop_ratio))
for(i in 1:length(Stop_ratio)){
  Stop_cum[i] <- sum(Stop_ratio[i:length(Stop_ratio)])
}

idea.x <- seq(0,0.999,0.001)
myu_P <- PL_mle(Stop_sec, min_sec)
myu_P # 1.739645
y1.est_P <- PL(idea.x, myu_P, min_sec)

lambda.est_P <- EXP_mle(Stop_sec, min_sec)
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
AIC1.exp <- -2*PL_llh(Stop_sec, min_sec, myu_P)+2*K1
AIC2.exp <- -2*EXP_llh(Stop_sec, min_sec, lambda.est_P)+2*K2
AIC <- c(AIC1.exp, AIC2.exp)
delta <- AIC - min(AIC)
w1.exp <- exp(-delta[1]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w2.exp <- exp(-delta[2]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w1.exp # 1
w2.exp # 0
## fit to PL
 
## fitting move
Move_lab <- sort(unique(Move_sec))
Move_count <- as.vector(table(Move_sec))
Move_ratio <- Move_count / sum(Move_count)
Move_cum <- rep(0, length(Move_ratio))
for(i in 1:length(Move_ratio)){
  Move_cum[i] <- sum(Move_ratio[i:length(Move_ratio)])
}

idea.x <- seq(0,0.999,0.001)
myu <- PL_mle(Move_sec, min_sec)
myu # 2.19908
y1.est <- PL(idea.x, myu, min_sec)

lambda.est <- EXP_mle(Move_sec, min_sec)
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
AIC1.exp <- -2*PL_llh(Move_sec, min_sec, myu)+2*K1
AIC2.exp <- -2*EXP_llh(Move_sec, min_sec, lambda.est)+2*K2
AIC <- c(AIC1.exp, AIC2.exp)
delta <- AIC - min(AIC)
w1.exp <- exp(-delta[1]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w2.exp <- exp(-delta[2]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w1.exp # 1
w2.exp # 4.017707e-158
## fit to PL

############################
### 2-2-2-1 Fitting to surrogate data
## generating step data
randomize_replace <- TRUE
sur_size <- 100000000
Treat <- "after"
Sex <- 0
Stop <- NULL
for(j in 1:length(f.namesplace_solo)){
  data.path <- f.namesplace_solo[j]
  d <- read.table(data.path, header=T)
  datepoint <- regexpr("170", data.path)[1]+7
  id <- (substr(data.path, datepoint, datepoint+13))
  if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
  treat <- (substr(data.path, datepoint+15, datepoint+20))
  if(treat=="after-"){treat <- "after"}
  if(treat!=Treat || sex!=Sex){ next; }
  
  time <- d[,1]; x <- d[,2]; y <- d[,3];  L <- length(x)
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

par(mfrow=c(2,1), pin=c(4,3)*0.8)


## fitting to pause
min_sec <- 0.2
Stop_lab <- sort(unique(Stop_sec))
Stop_count <- as.vector(table(Stop_sec))
Stop_ratio <- Stop_count / sum(Stop_count)
Stop_cum <- rep(0, length(Stop_ratio))
for(i in 1:length(Stop_ratio)){
  Stop_cum[i] <- sum(Stop_ratio[i:length(Stop_ratio)])
}

idea.x <- seq(0,0.999,0.0001)
myu <- PL_mle(Stop_sec, min_sec)
myu # 1.552449
y1.est <- PL(idea.x, myu, min_sec)

lambda.est <- EXP_mle(Stop_sec, min_sec)
lambda.est # 0.5998498
y2.est <- EXP(idea.x, lambda.est, min_sec)

plot(Stop_lab,Stop_cum, log="xy", main=paste("pause", Treat,Sex), axes=F,
     ylab="inverse cumulative frequency", xlab="pause time", xlim=c(0.2,20))
axis(1, c(0.2,0.5,1,2,5,10,20))
axis(2, c(0.001,0.01,0.1,1))
points(y1.est, rev(idea.x), col="red", lwd = 1.5, type="l")
points(y2.est, rev(idea.x), col="blue", lwd = 1.5, type="l")

K1 <- 1
K2 <- 1
AIC1.exp <- -2*PL_llh(Stop_sec, min_sec, myu)+2*K1
AIC2.exp <- -2*EXP_llh(Stop_sec, min_sec, lambda.est)+2*K2
AIC <- c(AIC1.exp, AIC2.exp)
delta <- AIC - min(AIC)
w1.exp <- exp(-delta[1]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w2.exp <- exp(-delta[2]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w1.exp # 1.518509e-275
w2.exp # 1
## fit to Exp

## move
Move_lab <- sort(unique(Move_sec))
Move_count <- as.vector(table(Move_sec))
Move_ratio <- Move_count / sum(Move_count)
Move_cum <- rep(0, length(Move_ratio))
for(i in 1:length(Move_ratio)){
  Move_cum[i] <- sum(Move_ratio[i:length(Move_ratio)])
}

idea.x <- seq(0,0.999,0.001)
myu <- PL_mle(Move_sec, min_sec)
myu # 13.77231
y1.est <- PL(idea.x, myu, min_sec)

lambda.est <- EXP_mle(Move_sec, min_sec)
lambda.est # 42.65267
y2.est <- EXP(idea.x, lambda.est, min_sec)

plot(Move_lab,Move_cum, log="xy", main=paste("move", Treat,Sex), axes=F,
     ylab="inverse cumulative frequency", xlab="move time", xlim=c(0.2,1))
axis(1, c(0.2,0.5,1))
axis(2, c(0.01,0.1,1))
points(y1.est, rev(idea.x), col="red", lwd = 1.5, type="l")
points(y2.est, rev(idea.x), col="blue", lwd = 1.5, type="l")

K1 <- 1
K2 <- 1
AIC1.exp <- -2*PL_llh(Stop_sec, min_sec, myu)+2*K1
AIC2.exp <- -2*EXP_llh(Stop_sec, min_sec, lambda.est)+2*K2
AIC <- c(AIC1.exp, AIC2.exp)
delta <- AIC - min(AIC)
w1.exp <- exp(-delta[1]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w2.exp <- exp(-delta[2]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
w1.exp # 1
w2.exp # 0
## not fit both
Move_lab
Move_ratio
#0.2 0.4 0.6 0.8
#0.9315589354 0.0621039290 0.0057034221 0.0006337136

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



