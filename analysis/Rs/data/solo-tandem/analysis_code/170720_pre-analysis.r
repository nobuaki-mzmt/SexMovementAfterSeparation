## Sex movement after separation

## 解析ノート 2017/7/17

## まずどのような処理をするかを、詳細な観察を通して決める。
## 動きのパターンを調べる上で重要な要素としては、
## ・角度
## ・速度
## ・止まっている時間
## ・step length
## の4つがあると考えられる。

## これらについてextractする方法を決定し、モデル上で動きを再現することを目的とする。

## これらについて方法が確定したら、生データの軌跡を解析用軌跡に変更し、
## 性・状況の効果を見るとともに、MSDの解析にうつる


######################

## setting 
anal_FPS <- 5

## data
setwd("F:\\tandem\\solo-tandem\\Rsperatus\\location_data\\5FPS")
setwd("C:\\Users\\惟暁\\Desktop\\location_data\\5FPS")
Folder.path <- getwd()
f.namesplace <- list.files(Folder.path, pattern=".txt",full.names=T)

## ある特定の個体に注目する
## "170526_163_F_1" (after)
j <- 1
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
if(treat == "after"){
  data <- d[1:(60*3*anal_FPS+1),]
}

# 軌跡の図示
plot(speed ~ time, data=data)
plot(y ~ x, data=data, type="l")



# 細かい部分に注目 (20 sec 分)
par(pin=c(3,3))
plot(y ~ x, data=data[1:101,], type="l" ,xlim=c(50,55), ylim=c(98,103))
# 図1
# 動画解析上のノイズがあったり、動かない方向転換があるため、
# めちゃくちゃ小さい距離の移動があるように思われる。
# 個体サイズが5mm程度あることを考慮すると、この動きは全て個体の身体の範囲内である。
# 動きのパターンを考える上では問題はない？

# 止まっている時間をみる
ds <- length(data[,1])
dis <- (c(((data[2:ds,2]-data[2:ds-1,2])^2 + (data[2:ds,3]-data[2:ds-1,3])^2)^0.5))
plot(dis)

# 0.2秒の移動距離が0.2以下(秒速1mm以下)の場合を止まっているとみなす
stop <- (dis < 0.2)

# stopは動いていないので、ここのデータを除いで軌跡をみると
plot(y ~ x, data=data[1:101,][c(TRUE,!stop),], type="l" ,xlim=c(50,55), ylim=c(98,103))



# 図2
# ノイズがキャンセルされてスムーズになっていることが分かる

data[na.omit((1:101)[c(TRUE,!stop)]),]

# このデータを使って、角度を再計算する
# この処理がない場合には、止まっている部分の角度が求められないため、
# 移動→止まる→方向転換
# のような状況に対応することが出来ない。

x <- data[na.omit((1:101)[c(TRUE,!stop)]),2]
y <- data[na.omit((1:101)[c(TRUE,!stop)]),3]
L <- length(x)

Ax <- (x[3:L-1] - x[3:L-2])
Bx <- (x[3:L] - x[3:L-1])
Ay <- (y[3:L-1] - y[3:L-2])
By <- (y[3:L] - y[3:L-1])
hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
angle <- acos(cos)*hugo

d3 <- data[na.omit((1:101)[c(TRUE,!stop)]),]
plot(y ~ x, data=d3, type="o" ,xlim=c(50,55), ylim=c(98,103))


# 図3


# Turchin’s ‘angle method’ (Turchin 1998) を用いて、角度でstep lengthを求めてみる
threshold <- pi/3 # e.g. 60度を閾値にしてみる
plot(y ~ x, data=d3[c(TRUE,abs(angle) >= threshold),], type="o" ,xlim=c(50,55), ylim=c(98,103))



# thresholdの決め方は、de jager 2011によると、
# 少しずつ閾値角度を大きくしていき、角度の自己相関(autocorrlation)が消えるとこを閾値としている


## ここまでは、ほとんど動かない期間に着目した。 "170526_163_F_1" (before)
## もう一つの動きのパターンとしてよく見られたloopingに今度は着目してみる。"170526_163_F_1" (after)

## "170526_163_F_1" (before)
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

# 部分に注目 (40 sec 分)
par(pin=c(3,3))
d3 <- data[1:201,]
plot(y ~ x, data=d3, type="l" ,xlim=c(20,110), ylim=c(50,140))
# 図1

# loopの部分を他の部分と区別したい

# 止まっている時間を消去し、angleの再設定
ds <- length(data[,1])
dis <- (c(((data[2:ds,2]-data[2:ds-1,2])^2 + (data[2:ds,3]-data[2:ds-1,3])^2)^0.5))
stop <- (dis < 0.2)

x <- data[,2]
y <- data[,3]
L <- length(x)
stop <- c(FALSE, stop)

Ax <- (x[3:L-1] - x[3:L-2])
Bx <- (x[3:L] - x[3:L-1])
Ay <- (y[3:L-1] - y[3:L-2])
By <- (y[3:L] - y[3:L-1])
hugo <- (Ax * By - Ay * Bx + 0.000001)/abs(Ax * By - Ay * Bx + 0.000001)
cos <- (Ax * Bx + Ay * By) / ((Ax^2 + Ay^2)*(Bx^2 + By^2))^0.5
angle <- acos(cos)*hugo

# loopの定義
# ・ 一定の角度で方向転換し続ける (threshold以下の角度変化を続ける: loop1)
# ・ 動いている (stop)
# ・ 途中で角度の符号が変わらない
# ・ 360度以上の回転をする回転をする

d3 <- data.frame(x=x, y=y, stop=stop, angle=c(NA,angle,NA))

thresh_loop <- pi*45/360
loop1 <- c(NA, diff(d3$angle) < thresh_loop)

# loop1の塊ないのdisplaceの性質に着目する
# numbering
loop_num <- 1
in_loop <- 0
for(i in 1:length(loop1)){
  if(loop1[i]){

}



plot(d3$angle)

plot(y ~ x, data=d3, type="o")


# 図3


# Turchin’s ‘angle method’ (Turchin 1998) を用いて、角度でstep lengthを求めてみる
threshold <- pi/3 # e.g. 60度を閾値にしてみる
plot(y ~ x, data=d3[c(TRUE,abs(angle) >= threshold),], type="o" ,xlim=c(50,55), ylim=c(98,103))




## 次は以上の方法を使って、
## 各生データに対して
## 止まっている時間の確定によるノイズキャンセル
## ノイズキャンセル後の速度・角度の計算
## を行う。

