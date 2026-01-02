## Sex movement after separation
## analysis for 2nd revise
## 2019/02/15 N. Mizumoto
## Data of Reticulitermes speratus

## Data
## solo tandem (30min before - 10min tandem - 30min after)
MMM <- "E:\\"
Laptop <- "C:\\Users\\nobuaki02\\"
place <- "Dropbox\\research\\with_dropbox\\papers and projects\\2017\\Sex movement after separation\\Rs\\solo-tandem\\location_data\\5FPS\\170623"
setwd(paste(MMM, place, sep=""))
Folder.path <- getwd()
f.namesplace_solo <- list.files(Folder.path, pattern=".txt",full.names=T)

#### setting 
anal_FPS <- 5
stop_thresh <- 0.7
cut_time <- 1
###############

## Truncated power-law
# r: 0-1
TP <- function(r,myu,xmin,xmax){
  return( ( xmax^(1-myu) - (1-r)*(xmax^(1-myu)-xmin^(1-myu) ) )^(1/(1-myu)) )
}


## Stretched exponention
# r (0-1)
SE <- function(r, lambda, beta, xmin){
  return( (xmin^beta - 1/lambda*log(1-r)  )^(1/beta))
}


TP_bin_LLF <- function(param, data, xmin, xmax){
  j = 1:(max(data)/0.2)
  dj <- j
  for(i in j){
    dj[i] = sum(round(data*5) == i)
  }
  return(-length(data)*log(xmin^(1-param)-xmax^(1-param)) +sum(dj*log( (xmin+(j-1)*0.2)^(1-param) - (xmin+j*0.2)^(1-param) )))
} 

SE_bin_LLF <- function(param, data, xmin){
  j = 1:(max(data)/0.2)
  dj <- j
  for(i in 1:length(j)){
    dj[i] = sum(round(data*5) == i)
  }
  return( length(data)*param[2]*xmin^param[1] + sum( dj*log( exp(-param[2]*(xmin+(j-1)*0.2)^param[1]) - exp(-param[2]*(xmin+j*0.2)^param[1]) ) ) )
} 


## Truncated power-law
TP_CDF <- function(myu, x, xmin, xmax){
  (xmin^(1-myu)-x^(1-myu)) / (xmin^(1-myu)-xmax^(1-myu))
}

## Stretched Exponential
SE_CDF <- function(lambda, beta, x, xmin){
  1-exp(-lambda*(x^beta-xmin^beta))
}

###############
#### Relationship between move duration and 
SexList <- c(0,1)
TreatList <- c("after","before")

###############
## generating step data
StepObtain <- function(Sex, Treat) {
  Stop_sec = Move_sec = Move_dis <- NULL
  for(j in 1:length(f.namesplace_solo)){
    data.path <- f.namesplace_solo[j]
    d <- read.table(data.path, header=T)
    datepoint <- regexpr("170", data.path)[1]+7
    id <- (substr(data.path, datepoint, datepoint+13))
    if(regexpr("_F_",data.path)>0) { sex <- 0;} else { sex <- 1; }
    treat <- (substr(data.path, datepoint+15, datepoint+20))
    if(treat=="after-"){treat <- "after"}
    if(treat!=Treat || sex!=Sex){ next; }
    
    d <- na.omit(d[,1:3])
    time <- d[,1]; x <- d[,2]; y <- d[,3]
    ini <- 1
    L <- length(x)
    dis <- (c(NA,((x[2:L]-x[2:L-1])^2 + (y[2:L]-y[2:L-1])^2)^0.5))
    stop <- (dis <= stop_thresh)
    stop[1] <- FALSE
    
    min3 <- 60*cut_time*anal_FPS+1
    if(min3>L){min3 <- L}
    
    if(Treat == "before"){
      # 29-30 min
      # look from the begining of behavior (move/pause) at 29 min
      # cut the last behavior
      mem <- stop[L]
      an_mem <- mem
      count <- 0
      while(mem == an_mem){
        count <- count -1
        an_mem <- stop[L+count]
      }
      L <- L+count
      
      mem <- stop[L-min3]
      an_mem <- mem
      count <- 0
      while(mem == an_mem){
        count <- count -1
        an_mem <- stop[L-min3+count]
      }
      ini <- L - min3 + count
      
      # 0-30 min
      # discard last behavior
      #mem <- stop[min3]
      #an_mem <- mem
      #count <- 0
      #while(mem == an_mem){
      #  count <- count -1
      #  an_mem <- stop[min3+count]
      #}
      #L <- min3+count
      
    } else {
      # 0-1 min
      # look until the behavior (move/pause) at 1 min finishes
      mem <- stop[min3]
      an_mem <- mem
      count <- 0
      while(mem == an_mem){
        count <- count +1
        an_mem <- stop[min3+count]
      }
      L <- min3+count
    }
    
    
    
    count <- 0
    in_stop <- 0
    stop_label <- rep(NA,L)
    move_label <- rep(NA,L)
    for(i in (ini+1):L){
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
    move_dis <- as.vector(tapply(dis[(ini+1):L], move_label[(ini+1):L], sum))
    Move_dis <- c(Move_dis, move_dis)
  }
  return(list(Stop_sec, Move_sec, Move_dis))
}

par(mfrow=c(1,4), pin=c(2,2))
k2 = 2
k3 = 2
Data <- StepObtain(SexList[k2], TreatList[k3])
#plot(Data[[2]], Data[[3]], xlab="Move duration (sec)", ylab="Moved distance (mm)")
plot(Data[[2]], Data[[3]]/Data[[2]], xlab="Move duration (sec)", ylab="Average displacements (mm)")
s <- Data[[3]]/Data[[2]]
Anova( lm(s~Data[[2]]))

#################
### fitting
Fitting <- function(StepData, k2, k3, k1){
  if(k1 == 1){ scheme <- "Pause"} else {scheme <- "Move"}
  if(k2 == 1){ SexLab <- "Female"} else {SexLab <- "Male"}
  
  idea.x <- seq(0,0.999,0.001)
  
  min_sec <- min(StepData)
  max_sec <- max(StepData) # Max x set as Max value
  Step_lab <- sort(unique(StepData))
  Step_count <- as.vector(table(StepData))
  Step_ratio <- Step_count / sum(Step_count)
  Step_cum <- rep(0, length(Step_ratio))
  for(i in 1:length(Step_ratio)){
    Step_cum[i] <- sum(Step_ratio[i:length(Step_ratio)])
  }
  
  # fit with truncated Power-law
  #param_TP <- optimize(TP_LLF, interval=c(1,10), data=StepData, xmin=min_sec, xmax=max_sec, maximum=T)
  param_TP <- optimize(TP_bin_LLF, interval=c(0,10), data=StepData, xmin=min_sec, xmax=max_sec, maximum=T)
  y1.est_TP <- TP(idea.x, param_TP$maximum[1], min_sec, max_sec)
  AIC1.TP <- -2*TP_bin_LLF(param_TP$maximum[1], StepData, min_sec, max_sec)+2*1
  
  (TP_bin_LLF(1.1, data=StepData, xmin=min_sec, xmax=max_sec))
  
  # fit with stretched exponential
  param_SE <- optim(c(0.1,0.1), SE_bin_LLF, data=StepData, xmin=min_sec, control = list(fnscale = -1), method="Nelder-Mead")
  #param_SE <- optim(c(0.2,0.2), SE_LLF, data=StepData, xmin=min_sec, control = list(fnscale = -1), method="Nelder-Mead")
  y2.est_SE <- SE(idea.x, param_SE$par[2], param_SE$par[1], min_sec)
  AIC2.SE <- -2*SE_bin_LLF(param_SE$par, StepData, min_sec)+2*2
  
  
  
  ## Selection
  AIC <- c(AIC1.TP, AIC2.SE)
  delta <- AIC - min(AIC)
  w1.TP_exp <- exp(-delta[1]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
  w2.SE_exp <- exp(-delta[2]/2)/(exp(-delta[1]/2)+exp(-delta[2]/2)) # Akaike weight
  
  if(w1.TP_exp > w2.SE_exp){
    judge <- "TP";}   else { judge <- "SE";
  }
  
  ## result plot
  if(PlotCombine){
    points(Step_lab,Step_cum, pch=19, col=(k2-1)*2+k3)
    lines(y1.est_TP, rev(idea.x), col="red", lwd = 1.5)
    lines(y2.est_SE, rev(idea.x), col="blue", lwd = 1.5)
  } else {
    plot(Step_lab,Step_cum, log="xy", main=paste(scheme, SexLab, TreatList[k3]), axes=F,
         ylab="inverse cumulative frequency", xlab="time", ylim=c(0.001,1), xlim=c(0.1, 1000))
    axis(1, c(0.1,1,10,100,1000))
    axis(2, c(0.001,0.01,0.1,1))
    lines(y1.est_TP, rev(idea.x), col="red", lwd = 1.5)
    lines(y2.est_SE, rev(idea.x), col="blue", lwd = 1.5)
  }
  
  ## goodness of fit
  GoF <- NULL
  if(judge == "TP"){
    Dvalue <- ks.test(StepData, TP_CDF, xmin=min_sec, xmax=max_sec, myu=param_TP$maximum[1])$statistic
    for(i in 1:2500){
      saro <- TP(runif(length(StepData),0,1), myu=param_TP$maximum[1], xmin=min_sec, xmax=max_sec)
      saro <- floor(saro*5)/5
      GoF <- c(GoF, as.numeric(ks.test(saro, TP_CDF, xmin=min_sec, xmax=max_sec, myu=param_TP$maximum[1])$statistic))
    }
  } else {
    Dvalue <- ks.test(StepData, SE_CDF, xmin=min_sec, beta=param_SE$par[1], lambda=param_SE$par[2])$statistic
    for(i in 1:2500){
      saro <- SE(runif(length(StepData),0,1), beta=param_SE$par[1], lambda=param_SE$par[2], min_sec)
      saro <- floor(saro*5)/5
      GoF <- c(GoF, as.numeric(ks.test(saro, SE_CDF, xmin=min_sec, beta=param_SE$par[1], lambda=param_SE$par[2])$statistic))
    }
  }
  
  GoF <- sum(GoF>Dvalue)/2500
  GoF
  return(c(SexLab, TreatList[k3], scheme, length(StepData), min_sec, max_sec, param_TP$maximum[1], param_SE$par[2], param_SE$par[1], w1.TP_exp, w2.SE_exp, judge, GoF) )
}

#################
### Calcuration
PlotCombine <- FALSE #TRUE #
par(mfrow=c(2,4), pin=c(4,3)*0.5)
#par(mfrow=c(1,2), pin=c(4,3)*0.8)
res <- NULL
for(k1 in 1:2){
  if(PlotCombine){
    if(k1 == 1){
      plot(0.2,1, log="xy", main=paste("Pause"), axes=F,
           ylab="inverse cumulative frequency", xlab="time",
           ylim=c(0.001,1), xlim=c(0.2, 500), type="n")
    } else {
      plot(0.2,1, log="xy", main=paste("Move"), axes=F,
           ylab="inverse cumulative frequency", xlab="time",
           ylim=c(0.001,1), xlim=c(0.2, 200), type="n")
    }
    axis(1, c(0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000))
    axis(2, c(0.001,0.01,0.1,1))
  }
  for(k2 in 1:2){
    for(k3 in 1:2){
      Data <- StepObtain(SexList[k2], TreatList[k3])
      res <- rbind(res, Fitting(Data[[k1]], k2, k3, k1))
    }
  }
}

res
#write.table(res, "E:\\Rsres.csv")

