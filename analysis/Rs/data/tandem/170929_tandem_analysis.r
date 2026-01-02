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

ID = sep_time = sep_male_speed = sep_female_speed = sep_dis <- NULL
fem_speed_dynamics = MF_dis_dynamics <- NULL

for(j in 1:length(f.namesplace)){
  data.path <- f.namesplace[j]
  d <- data.frame(as.matrix(fread(data.path, header=T)))
  datepoint <- regexpr("170", data.path)[1]
  id <- (substr(data.path, datepoint, datepoint+11))
  date <- (substr(data.path, datepoint, datepoint+5))
  colony <- (substr(data.path, datepoint+7, datepoint+9))
  iter <- (substr(data.path, datepoint+11, datepoint+11))
  
  d[,1] <- d[,1]/30
  if(length(d[,1]) >= 30*60*60+1){
    anal_time <- seq(1, 30*60*60+1, 30/anal_FPS)
  } else {
    anal_time <- seq(1, length(d[,1]), 30/anal_FPS)
  }
  
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
  
  time <- da$position
  Male_X <- da$x1
  Male_Y <- da$y1
  Female_X <- da$x0
  Female_Y <- da$y0
  Male_speed <- c(NA, sqrt ( diff(Male_X)^2 + diff(Male_Y)^2 ) * 5)
  Female_speed <- c(NA, sqrt ( diff(Female_X)^2 + diff(Female_Y)^2 ) * 5)
  
  L <- length(time)
  Man <- atan( (Male_Y[2:L] - Male_Y[2:L-1])/(Male_X[2:L] - Male_X[2:L-1]) )
  Man[ Male_X[2:L] - Male_X[2:L-1] > 0 & !is.na(Man)] <- Man[ Male_X[2:L] - Male_X[2:L-1] > 0 & !is.na(Man)] + pi
  Man[Man<0 & !is.na(Man)] <- Man[Man<0 & !is.na(Man)] + 2*pi
  Man <- (Man / pi * 180)
  
  Feman <- atan( (Female_Y[2:L] - Female_Y[2:L-1])/(Female_X[2:L] - Female_X[2:L-1]) )
  Feman[ Female_X[2:L] - Female_X[2:L-1] > 0 & !is.na(Feman)] <- Feman[ Female_X[2:L] - Female_X[2:L-1] > 0 & !is.na(Feman)] + pi
  Feman[Feman<0 & !is.na(Feman)] <- Feman[Feman<0 & !is.na(Feman)] + 2*pi
  Feman <- (Feman / pi * 180)
  
  MF_dis <- sqrt( (Male_X-Female_X)^2 + (Male_Y-Female_Y)^2 )
  
  # plot(Female_X, Female_Y, type="l")
  #plot(MF_dis)
  #plot(time, MF_dis)
  #plot(time[time > 2500 & time < 3000], MF_dis[time > 2500 & time < 3000])
  #plot(Female_speed)
  
  
  ## 距離7mm以下が3sec以上続いて、初めてタンデムとみなす
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
  
  ## 近くにいる期間に、進んだ距離が30mm以下ならタンデムとみなさない
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
  
  
  ## sep中2個体の方向の差の最大が30度以下ならsepとみなさない
  for(i in 1:max(separate_label)){
    an_dif <- abs(Man[separate_label==i] - Feman[separate_label==i])
    rev_an_dif <- abs(360-an_dif)
    an_dif[ an_dif > rev_an_dif & !is.na(an_dif)] <- rev_an_dif[ an_dif > rev_an_dif & !is.na(an_dif)]
    if( max(an_dif, na.rm=T) < 45 ){
      separate[separate_label==i] <- FALSE
    }
  }
  
  
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
  
  if(separate[L]){
    separate[separate_label==max(separate_label)] <- FALSE
  }
  
  if(max(separate_label)==0){next;}
  
  for(i in 1:max(separate_label)){
    sep_time <- c(sep_time, length(MF_dis[separate_label == i]) / 5 )
    sep_male_speed <- c(sep_male_speed, mean(Male_speed[separate_label == i]) )
    sep_female_speed <- c(sep_female_speed, mean(Female_speed[separate_label == i]) )
    sep_dis <- c(sep_dis, mean(MF_dis[separate_label == i]) )
  }
  ID <- c(ID, rep(j, max(separate_label)))
  
  for(i in 1:max(separate_label)){
    fem_speed_dynamics <- c(fem_speed_dynamics, list( Female_speed[separate_label == i] ) )
    MF_dis_dynamics <- c(MF_dis_dynamics, list( MF_dis[separate_label == i] ) )
  }
  
  #plot(time,separate, type="o")
  #plot(time,separate, type="o", xlim=c(0,300))
  print(j)
}

res <- cbind(ID,sep_time, sep_female_speed, sep_male_speed)
plot(sep_female_speed, sep_time)
plot(sep_male_speed, sep_time)
plot(sep_female_speed-sep_male_speed, sep_time)
plot(abs(sep_male_speed-sep_female_speed), sep_time)
cor.test(res[,3], res[,2])

plot(1, type="n", xlim=c(0,150), ylim=c(0,50), xlab="time (sec)", ylab="speed (mm/sec)")
for(i in 1:length(fem_speed_dynamics)){
  df <- fem_speed_dynamics[[i]]
  points(seq(0.2, length(df)*0.2, 0.2), df, type="l")
}

dr = drn <- rep(0, 150*5)
for(i in 1:length(fem_speed_dynamics)){
  df <- fem_speed_dynamics[[i]]
  dr[1:length(df)] <- dr[1:length(df)]+df
  drn[1:length(df)] <- drn[1:length(df)]+1
}
plot( dr[drn>0]/drn[drn>0] )

sep_def_dis <- NULL
for(i in 1:length(fem_speed_dynamics)){
  df <- fem_speed_dynamics[[i]]
  dd <- MF_dis_dynamics[[i]]
  sep_def_dis <- c(sep_def_dis, dd[df/5 <= 0.5][1])
}
mean(sep_def_dis, na.rm=T)

sep_def_dis <- NULL
for(i in 1:length(fem_speed_dynamics)){
  dd <- MF_dis_dynamics[[i]]
  sep_def_dis <- c(sep_def_dis, dd)
}
truehist(sep_def_dis)
y <- density(sep_def_dis)
plot(y)
y$x[y$y == max(y$y)]

