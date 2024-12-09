## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----rock---------------------------------------------------------------------
library(TSA)
plot(rock)

## -----------------------------------------------------------------------------
r = lm(perm~.,rock)
summary(r)

## -----------------------------------------------------------------------------
qqnorm(r$res)

## -----------------------------------------------------------------------------

p2=c(1,-1.2,0.6,-0.2,0.05)
r2 = polyroot(p2)
r2

## -----------------------------------------------------------------------------
a1=ARMAacf(ar=c(1.2,-0.6,0.2,-0.05),lag=15)
plot(a1,type='o')

## -----------------------------------------------------------------------------
data(arma11.s)
plot(arma11.s,type='o',ylab=expression(Yt))

## -----------------------------------------------------------------------------
acf(arma11.s)

## -----------------------------------------------------------------------------
pacf(arma11.s)

## -----------------------------------------------------------------------------
sigma = c(1,2,3,4,10,15) #选取6个sigma
for(i in 1:length(sigma)){
  x1 = rnorm(1e4,0,sigma[i]) # 1000个服从正态分布的随机数
  x2 = rnorm(1e4,0,sigma[i]) # 两个独立同正态分布的随机向量
  r = (x1^2+x2^2)^0.5 # 生成服从Rayleigh分布的随机数
  hist(r,prob=TRUE, main = paste('sigma=',sigma[i]))   # 样本直方图
  
  x = seq(0,4*sigma[i],length.out=1e4)
  f = (x/sigma[i]^2)*exp(-(x^2)/(2*sigma[i]^2))
  lines(x,f) # 密度函数
}


## -----------------------------------------------------------------------------
pp = seq(0,1,0.05) 
p = 0.75
x1 = rnorm(1000,0,1)
x2 = rnorm(1000,3,1)
m = p*x1+(1-p)*x2   #N(0,1)和N(3,1)分布的混合
hist(m,prob=TRUE,breaks = seq(-5,5,length.out=100),main = 'p1=0.75')         #当p=0.75时直方图
for(i in 1:length(pp)){
  mm = pp[i]*x1+(1-pp[i])*x2 # 不同p值下的直方图
  hist(mm,prob = TRUE, breaks = seq(-5,8,length.out = 100),main = paste('p1=',pp[i]))
}

## -----------------------------------------------------------------------------
lambda = 5
shape = 2
scale = 4
size = 10000
eps = 0.1
xss = numeric(1e4)  # 初始化
ns = numeric(1e4)
t = 10

n = qpois(1-eps, lambda = lambda * t) #泊松过程：次数

# 生成服从泊松过程的Nt：
for(i in 1:size){
  
  Tn = rexp(1000, lambda) #  到达时间间隔为参数lambda的指数分布
  Sn = cumsum(Tn) # 累计求和
  ns[i] = min(which(Sn > t)) - 1
  
}

# 泊松-伽马过程：
xss = sapply(ns, function (n) {
  ys = c(rgamma(n = n, shape = shape, scale = scale))
  sum(ys[1:n])
})

means = mean(xss) # 样本均值和样本方差
vars = var(xss)

meant = lambda * t * shape * scale #理论均值和方差
vart = (shape + 1) * shape * scale^2

means;meant;vars;vart

## -----------------------------------------------------------------------------
m = 1e5 #times of simulations
t = runif(m) 
b = mean(t^2*(1-t)^2) #calculate the constant B(3,3)
for(x in seq(0.1,0.9,0.1)){
  u = runif(m,0,x)
  theta.hat = mean(x*u^2*(1-u)^2)/b
  theta.t = pbeta(x,3,3)
  print(c(theta.hat,theta.t, theta.hat/theta.t))
  
}

## -----------------------------------------------------------------------------
sigma2 = 4
MC.Phi <- function(x, R = 10000, antithetic = FALSE) {
  u <- runif(R/2)
  if (antithetic) v <- 1 - u else v <- runif(R/2)
  u <- c(u, v)
  g <- x * exp(-(u * x)^2 / 2*sigma2) # x*u ~ U(0,x)
  cdf <- mean(g) / sigma2
  cdf
}

m <- 1000
MC1 <- MC2 <- numeric(m)
for (x in seq(1,10)){
  for (i in 1:m){
    MC1[i] <- MC.Phi(x, R = 1000, antithetic = FALSE)
    MC2[i] <- MC.Phi(x, R = 1000, antithetic = TRUE)
  }
  print(c(sd(MC1),sd(MC2),sd(MC2)/sd(MC1)))
}


## -----------------------------------------------------------------------------
x = seq(1,5,0.1)
g = x^2*exp(-x^2/2)/sqrt(2*pi)
f1 = exp(-x^2/2)/sqrt(2*pi) # x~N(0,1) in f1
f2 = x*exp(-x^2/2)  # x~Rayleigh(1) in f2, x>=0
gs <- c(expression(g(x)==x^2*e^(-x^2/2)/sqrt(2*pi)),
                   expression(f[1](x)==e^(-x^2/2)/sqrt(2*pi)),
                  expression(f[2](x)==x*e^(-x^2/2)))
plot(x,g,type = 'l',lwd = 2)
lines(x,f1,lwd = 2, col=2,lty =2)
lines(x,f2,lwd = 2, col = 3,lty = 3 )
legend("topright",legend= gs, lwd = 2, lty = 1:3)


## -----------------------------------------------------------------------------
m3 = 1e4
est = sd = numeric(2)
g <- function(x) {  #原函数
  x^2*exp(-x^2/2)/sqrt(2*pi) * (x > 1)
}

x = rnorm(m,0,1) # using f1:N(0,1)
fg = g(x)/exp(-x^2/2)/sqrt(2*pi)
est[1] = mean(fg)
sd[1] = sd(fg)

x1 = rnorm(1e4,0,1) # 1000个服从正态分布的随机数
x2 = rnorm(1e4,0,1) # 两个独立同正态分布的随机向量
r = (x1^2+x2^2)^0.5 # 生成服从Rayleigh(1)分布的随机数
fg = g(r)/r*exp(-r^2/2)
est[2] = mean(fg)
sd[2] = sd(fg)
print("f1:")
print(sd[1])
print("f2:")
print(sd[2])


## ----echo=TRUE----------------------------------------------------------------
count = 0
c = numeric(5)
quick_sort = function(arr,count) {
  
  if (length(arr) <= 1) {
    return(arr)
  }
  a = arr[1]
  left = arr[arr < a]
  right = arr[arr > a]
  count =  count+length(left)+length(right)
  count <<- count
  #print(count)
  return(c(quick_sort(left,count), a, quick_sort(right,count)))
}

## ----echo=TRUE----------------------------------------------------------------
count = 0
for (i in seq(1,5)){
  n = c(1e4,2e4,4e4,6e4,8e4)
  quick_sort(sample(1:n[i],size=n[i],replace = FALSE),count = count)
  c[i] = count
  print(c(c[i],n[i]*log(n[i])))
}
plot(c,n*log(n),xlab='an')
an = c
tn = n*log(n)
modelreg = lm(tn~an)

abline(modelreg)


## -----------------------------------------------------------------------------
# 设置参数
n <- 1000  # 样本大小
m <- 10000  # Monte Carlo模拟次数
quantile_probs <- c(0.025, 0.05, 0.95, 0.975)

# 定义偏度函数，避免负数的平方根
calc_skewness <- function(x) {
  mean_x <- mean(x)
  sd_x <- sd(x)
  n <- length(x)
  skewness <- sum((x - mean_x)^3) / (n * sd_x^3)
  if (skewness >= 0) {
    return(sqrt(skewness))  # 取平方根
  } else {
    return(NA)  # 返回NA避免负值
  }
}

# Monte Carlo模拟
set.seed(123)
skew_vals <- numeric(m)
for (i in 1:m) {
  sample_data <- rnorm(n)
  skew_vals[i] <- calc_skewness(sample_data)
}

# 去除产生NA的情况
skew_vals <- skew_vals[!is.na(skew_vals)]

# 估计分位数
estimated_quantiles <- quantile(skew_vals, probs = quantile_probs)

# 计算标准误差，根据公式(2.14)
calc_standard_error <- function(q, n, f_xq) {
  sqrt(q * (1 - q) / (n * f_xq^2))
}

# 获取分位数处的密度值
f_xq_vals <- dnorm(qnorm(quantile_probs))

# 标准误差
se_vals <- mapply(calc_standard_error, quantile_probs, MoreArgs = list(n = n, f_xq = f_xq_vals))

# 计算大样本近似的分位数
large_sample_quantiles <- qnorm(quantile_probs, mean = 0, sd = sqrt(6 / n))

# 输出结果
cat("Monte Carlo估计的分位数:\n", estimated_quantiles, "\n")
cat("标准误差:\n", se_vals, "\n")
cat("大样本近似分位数:\n", large_sample_quantiles, "\n")



## -----------------------------------------------------------------------------
asstest <- function(seed){
seed <- set.seed(123)
x <- rnorm(20,2,10)
sigma <- rnorm(20,5,50)  #生成正态分布的随机数
y <- 3*x+sigma
cor(x,y)
pearson <- cor.test(x,y) #变量相关性显著性检验
kendall <- cor.test(x,y,method = 'kendall')
spearman <- cor.test(x,y,method = 'spearman')
data.frame(x,y)
return(list(pearson=pearson,kendall=kendall,spearman=spearman))
}
asstest()

## -----------------------------------------------------------------------------
# Given powers for the two tests
p1 <- 0.651  # Power for Z-test
p2 <- 0.676  # Power for paired-t test
n1 <- 10000  # Number of experiments for Z-test
n2 <- 10000  # Number of experiments for paired-t test

# Pooled proportion
p_hat <- (p1 * n1 + p2 * n2) / (n1 + n2)

# Standard error
se <- sqrt(p_hat * (1 - p_hat) * (1/n1 + 1/n2))

# Test statistic (z-score)
z <- (p1 - p2) / se

# Two-tailed p-value
p_value <- 2 * pnorm(-abs(z))

# Output the z-score and p-value
cat("Z-score:", z, "\n")
cat("P-value:", p_value, "\n")

# Decision based on p-value at 0.05 significance level
if (p_value < 0.05) {
  cat("Reject the null hypothesis: the powers are significantly different.\n")
} else {
  cat("Fail to reject the null hypothesis: no significant difference in powers.\n")
}


## -----------------------------------------------------------------------------
N = 1000; N0 = 950; N1 = 50; a = 0.1; M = 1e4
temp1 = temp2 = numeric(N)
fwer_bon = fwer_bh = fdr_bon = fdr_bh = tpr_bon = tpr_bh = numeric(M) #初始化

for (m in 1:M){  # 10000次遍历
  
  pnull = runif(950); palt = rbeta(50,0.1,1); pall = c(pnull,palt) #生成原假设和备择假设下p值的随机数
  padj_bon = p.adjust(pall, method = 'bonferroni') #用两种方法计算调整的p值
  padj_bh = p.adjust(pall, method = 'fdr')
  fwer_bon[m] = mean(1-(1-padj_bon)**N)  # 计算两种修正p值下的fwer
  fwer_bh[m] = mean(1-(1-padj_bh)**N)

  p_bon_sort = sort(padj_bon) # 计算bonferroni方法的fdr值
  for(i in seq(N)){
    if(p_bon_sort[i] <= i*a/N){
      temp1[i] = p_bon_sort[i]
    }
  }
  fdr_bon[m] = N0*max(temp1)/which(temp1 == max(temp1))

  p_bh_sort = sort(padj_bh) # 计算B-H方法的fdr值
  for(i in seq(N)){
    if(p_bh_sort[i] <= i*a/N){
      temp2[i] = p_bh_sort[i]
    }
  }
  fdr_bh = N0*max(temp2)/which(temp2 == max(temp2))

  paltadj_bon = p.adjust(palt, method = 'bonferroni')
  paltadj_bh = p.adjust(palt, method = 'fdr')
  tpr_bon[m] = length(which(paltadj_bon<a))/N1  # 分别计算tpr
  tpr_bh[m] = length(which(paltadj_bh<a))/N1
  
}

t = rbind(
  FWER = c(mean(fwer_bon), mean(fwer_bh)),  # 输出表格
  FDR = c(mean(fdr_bon),mean(fdr_bh)),
  TPR = c(mean(tpr_bon),mean(tpr_bh))
)
colnames(t) = c('Bonferroni correction','B-H correction')
knitr::kable(t,format = 'markdown')

## -----------------------------------------------------------------------------
obs = c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487) # aircondit数据
B = 1e4; thetastar = numeric(B)
theta = mean(obs)  # 样本均值

for(b in 1:B){
  obsstar = sample(obs, replace = TRUE) # 每一次抽样，共抽10000次
  thetastar[b] = mean(obsstar)   # 计算抽取的样本均值
}

round(c(
  mle = 1/mean(obs), bias = mean(thetastar)-theta, se.boot = sd(thetastar), se.samp = sd(obs)/sqrt(length(obs))
),3)

## -----------------------------------------------------------------------------
m<-1e3;library(boot);set.seed(12345) 
obs = c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)  # aircondit数据
lambda = 1/mean(obs)  # 指数分布lambda的估计
boot.mean <- function(x,i) 1/mean(x[i])  # 定义均值函数
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2) 
for(i in 1:m){
  R = sample(obs, replace = TRUE)  #抽取样本
  de <- boot(data=R,statistic=boot.mean, R = 999) # 自举法
  ci <- boot.ci(de,type=c("norm","basic","perc","bca"))  #自举法置信区间
  ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5];ci.bca[i,]<-ci$bca[4:5]
}

ciinterval = round(rbind(norm = c(mean(ci.norm[,1]),mean(ci.norm[,2])),  # 输出4种方法的置信区间
                         basic = c(mean(ci.basic[,1]),mean(ci.basic[,2])),
                         perc = c(mean(ci.perc[,1]),mean(ci.perc[,2])),
                         BCa = c(mean(ci.bca[,1]), mean(ci.bca[,2]))
                         ),3)
colnames(ciinterval) = c('lowbound', 'upbound')
ciinterval

cat('norm =',mean(ci.norm[,1]<=lambda & ci.norm[,2]>=lambda),   # 计算每种方法的正确率
    'basic =',mean(ci.basic[,1]<=lambda & ci.basic[,2]>=lambda),
    'perc =',mean(ci.perc[,1]<=lambda & ci.perc[,2]>=lambda),
    'BCa =',mean(ci.bca[,1]<=lambda & ci.bca[,2]>=lambda))


## -----------------------------------------------------------------------------
library(MASS); library(bootstrap); data(scor)
# 定义生成函数
datagene = function(data){
  cov_matrix = cov(data)
  eigenvalues = eigen(cov_matrix)$values
  theta = eigenvalues[1]/sum(eigenvalues)
  return(theta)
}


# 定义推断函数
statisinference = function(data){
  n = nrow(data) # scor数据的样本量
  theta_jack = numeric(n)
  for(i in 1:n){
    scor_jack = data[-i,]   # 将scor数据去除第i个观测值
    theta_jack[i] = datagene(scor_jack) # 计算新的theta
  }
  return(theta_jack)
}


# 定义结果函数，输出偏差和标准差
resultout = function(data){
  n = nrow(data)
  # 计算偏差估计
  bias = (n-1)*(mean(statisinference(data))-datagene(data))
  # 计算标准差估计
  se = sqrt((n-1)/n*sum((statisinference(data)-mean(statisinference(data)))^2))
  round(c(bias = bias, se = se),5)
}

# 调用结果函数
resultout(scor)


## -----------------------------------------------------------------------------
library(DAAG, quietly = TRUE); data(ironslag)
## 定义生成函数
datagene2 = function(){
  magnetic = ironslag$magnetic
  
  n = length(magnetic)
  return(n)
}

## 定义推断函数
modelreg = function(){
  n = datagene2(); magnetic = ironslag$magnetic; chemical = ironslag$chemical
  e1 = numeric(n)
  for(k in 1:n){
    y = magnetic[-k] # 交叉验证，剔除第k个值
    x = chemical[-k]
    
    # 四个模型拟合
    J1 = lm(y ~ x)
    yhat1 = J1$coef[1] + J1$coef[2] * chemical[k]
    e1[k] = magnetic[k] -yhat1
    
    round(c(linear_bias = mean(e1^2)),4)
    
    return(summary(J1))
  }
}

## 定义输出函数
result2out = function(){
  n = datagene2(); magnetic = ironslag$magnetic; chemical = ironslag$chemical; e2 = numeric(n)
  for(k in 1:n){
    y = magnetic[-k] # 交叉验证，剔除第k个值
    x = chemical[-k]
    
    J2 = lm(y ~ x + I(x^2))
    yhat2 = J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
    e2[k] = magnetic[k] - yhat2
  }
  round(c(quadratic = mean(e2^2)), 4)
  return(summary(J2))
}

result3out = function(){
  n = datagene2(); magnetic = ironslag$magnetic; chemical = ironslag$chemical; e3 = numeric(n)
  for(k in 1:n){
    y = magnetic[-k] # 交叉验证，剔除第k个值
    x = chemical[-k]
    
    J3 = lm(log(y) ~ x)
    yhat3 = J3$coef[1] + J3$coef[2] * chemical[k]
    e3[k] = magnetic[k] - yhat3
  }
  round(c(expo = mean(e3^2)), 4)
  return(summary(J3))
}

result4out = function(){
  n = datagene2(); magnetic = ironslag$magnetic; chemical = ironslag$chemical; e4 = numeric(n)
  for(k in 1:n){
    y = magnetic[-k] # 交叉验证，剔除第k个值
    x = chemical[-k]
    
    J4 = lm(y ~ x + I(x^2) + I(x^3))
    yhat4 = J4$coef[1] +J4$coef[2]*chemical[k] + J4$coef[3]*chemical[k]^2 + J4$coef[4]*chemical[k]^3
    e4[k] = magnetic[k] - yhat4
  }
  round(c(cubic = mean(e4^2)),4)
  return(summary(J4))
  #return(mean(e4^2))
}

## 调用输出函数
modelreg()
result2out()
result3out()
result4out()
 

## -----------------------------------------------------------------------------
## 定义生成函数
rep = 1000
genexs = function(){
  feed = chickwts$feed
  feed1 = "soybean"
  attach(chickwts)
  xs = sort(as.vector(weight[feed == "soybean"])) # 生成soybean的重排
  detach(chickwts)
  return(xs)
}

geneys = function(){
  feed = chickwts$feed
  attach(chickwts)
  ys = sort(as.vector(weight[feed == 'linseed'])) # 生成linseed的重排
  detach(chickwts)
  return(ys)
}


## 定义推断函数
permu_pro = function(){
  library(cramer)
  xs = genexs()
  ys = geneys()
  zs = c(xs, ys)
  n1 = length(xs)
  n2 = length(ys)
  n = n1 + n2

  Ts = numeric(rep)
  for (i in 1:rep) {
    ks = sample(1:n, n1, replace = FALSE)
    zs1 = zs[ks]
    zs2 = zs[-ks]
    Ts[i] = cramer.test(zs1, zs2)$statistic
  }

  cra = cramer.test(xs, ys)
  T.hat = cra$statistic
  round(c(p.hat = mean(abs(T.hat) < abs(c(T.hat, Ts)))),4)
  
  
}

permu_sum = function(){
  library(cramer)
  xs = genexs()
  ys = geneys()
  zs = c(xs, ys)
  n1 = length(xs)
  n2 = length(ys)
  n = n1 + n2

  Ts = numeric(rep)
  for (i in 1:rep) {
    ks = sample(1:n, n1, replace = FALSE)
    zs1 = zs[ks]
    zs2 = zs[-ks]
    Ts[i] = cramer.test(zs1, zs2)$statistic
  }

  cra = cramer.test(xs, ys)
  T.hat = cra$statistic
  #round(c(p.hat = mean(abs(T.hat) < abs(c(T.hat, Ts)))),4)
  return(cra)
  #hist(Ts)
  
}

permu_hist = function(){
  library(cramer)
  xs = genexs()
  ys = geneys()
  zs = c(xs, ys)
  n1 = length(xs)
  n2 = length(ys)
  n = n1 + n2

  Ts = numeric(rep)
  for (i in 1:rep) {
    ks = sample(1:n, n1, replace = FALSE)
    zs1 = zs[ks]
    zs2 = zs[-ks]
    Ts[i] = cramer.test(zs1, zs2)$statistic
  }

  cra = cramer.test(xs, ys)
  T.hat = cra$statistic
  
  hist(Ts)
  
}


# 调用函数
permu_pro()
permu_sum()
permu_hist()

## -----------------------------------------------------------------------------
## 定义生成函数
data(chickwts)
datagene4= function(){
  attach(chickwts)
  soybean = chickwts$weight[chickwts$feed=="soybean"]
  linseed = chickwts$weight[chickwts$feed=="linseed"]
  detach(chickwts)
  n = length(soybean)
  m = length(linseed)

  tmp = min(n, m)
  soybean = sort(soybean[1:tmp])
  linseed = sort(linseed[1:tmp])
  
  return(list(soybean = soybean, linseed = linseed))
}

## 定义推断函数
spearman_fun = function(B){
  soybean = datagene4()$soybean
  linseed = datagene4()$linseed
  zs = c(soybean, linseed)
  spearman_cortest = cor.test(x = soybean, y = linseed, method = "spearman")

  k = length(zs)

  rhos = numeric(B)

  for (b in 1:B) {
    i = sample(1:k, k/2, replace = FALSE)
    xs = zs[i]
    ys = zs[-i]
    rhos[b] = cor(x = xs, y = ys, method = "spearman")
  }

  hist(rhos, breaks = 100)
  theta.hat = spearman_cortest$estimate
  p.hat = mean(abs(rhos) > abs(theta.hat))
  return(spearman_cortest)
}

spearman_fun2 = function(B){
  soybean = datagene4()$soybean
  linseed = datagene4()$linseed
  zs = c(soybean, linseed)
  spearman_cortest = cor.test(x = soybean, y = linseed, method = "spearman")

  k = length(zs)

  rhos = numeric(B)

  for (b in 1:B) {
    i = sample(1:k, k/2, replace = FALSE)
    xs = zs[i]
    ys = zs[-i]
    rhos[b] = cor(x = xs, y = ys, method = "spearman")
  }

  theta.hat = spearman_cortest$estimate
  p.hat = mean(abs(rhos) > abs(theta.hat))
  return(p.hat)
  }

## 调用推断函数，画直方图
spearman_fun(1000)

## 定义结果函数
resultout4 = function(){
  spearman_cortest = spearman_fun(1000)
  theta.hat = spearman_cortest$estimate
  p.value = spearman_cortest$p.value
  
  
  round(c( theta.hat = theta.hat, p.value = p.value),3)
}

resultout4() # 调用函数
spearman_fun2(1000) # 输出p.hat


## -----------------------------------------------------------------------------
## 定义生成函数：柯西分布
 f <- function(x) {
        
        return(1/(pi*(1+x^2)))
    }

# 定义统计推断函数：
statinf = function(){
    m <- 10000
    x <- numeric(m)
    x[1] <- rnorm(1, mean=1)
    k <- 0
    u <- runif(m)

    for (i in 2:m) {
        xt <- x[i-1]
        y <- rnorm(1, mean = xt)
        num <- f(y) * dnorm(xt, mean = y)
        den <- f(xt) * dnorm(y, mean = xt)
        if (u[i] <= num/den){
          x[i] <- y
        } else {
          x[i] <- xt
          k <- k+1     #y is rejected
        }
    }
    print(round(c('p.reject' = k/m),4))
    return(x)
}
    
# 定义结果输出函数
out = function(){
    m = 1e4
    b <- 1001      #discard the burn-in sample
    y <- statinf()[b:m]
    plot(1001:m, y, type = 'l', main='', ylab='x')
    # a <- ppoints(100)
    # Q <- quantile(x, a)
    hist(statinf(),probability = TRUE, breaks = 100)
    
    round(c('deciles of generation observation'=quantile(y,0.1), 
            'deciles of standard qcauchy'=qcauchy(0.1)),4)
}

out()
    

## -----------------------------------------------------------------------------

# 定义生成函数，生成标准柯西分布的马尔科夫链

    cauchy.chain <- function(sigma, N, X1) {
        #generates a Metropolis chain for Cauchy(0,1)
        #with Normal(X[t], sigma) proposal distribution
        #and starting value X1
        x <- rep(0, N)
        x[1] <- X1
        u <- runif(N)

        for (i in 2:N) {
            xt <- x[i-1]
            y <- rnorm(1, xt, sigma)     #candidate point
            r1 <- dcauchy(y, 0, 1) * dnorm(xt, y, sigma)
            r2 <- dcauchy(xt, 0, 1) * dnorm(y, xt, sigma)
            r <- r1 / r2
            if (u[i] <= r) x[i] <- y else
                 x[i] <- xt
            }
        return(x)
    }
  

## 定义推断函数，即计算G_R的统计量
Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
        }

## 定义结果输出函数，输出四条链轨迹图及R-hat统计量序列
out2 = function(){
  
    sigma <- 2     #parameter of proposal distribution
    k <- 4          #number of chains to generate
    n <- 15000      #length of chains
    b <- 1000       #burn-in length

    #choose overdispersed initial values
    x0 <- c(-10,-5,5,10)

    #generate the chains
    set.seed(12345)
    X <- matrix(0, nrow=k, ncol=n)
    for (i in 1:k)
        X[i, ] <- cauchy.chain(sigma, n, x0[i])
    
    
    #compute diagnostic statistics
    psi <- t(apply(X, 1, cumsum))
    for (i in 1:nrow(psi))
        psi[i,] <- psi[i,] / (1:ncol(psi))

    #plot psi for the four chains

    for (i in 1:k)
      if(i==1){
        plot((b+1):n,psi[i, (b+1):n],ylim=c(-5,5), type="l",
            xlab='Index', ylab=bquote(phi))
      }else{
        lines(psi[i, (b+1):n], col=i)
    }
    par(mfrow=c(1,1)) #restore default
    
    #plot the sequence of R-hat statistics
    rhat <- rep(0, n)
    for (j in (b+1):n){
      rhat[j] <- Gelman.Rubin(psi[,1:j])
    }
        
    plot(rhat[(b+1):n], type="l", xlab="", ylab="R", ylim = c(1,4))
    abline(h=1.2, lty=2)
}
out2()
    

## -----------------------------------------------------------------------------
# 设置参数
n <- 10      # 样本大小
a <- 2       # Beta分布的参数
b <- 2       # Beta分布的参数
N <- 5000  # Gibbs采样的迭代次数
burn = 1000

# 1. 生成Gibbs采样的数据
generate_samples <- function(n, a, b, iterations) {
  x_samples <- numeric(iterations)
  y_samples <- numeric(iterations)
  y <- 0.5  # 初始化 y 的初始值
  
  for (i in 1:iterations) {
    x <- rbinom(1, n, y)              # 从 x | y 的条件分布采样
    y <- rbeta(1, x + a, n - x + b)   # 从 y | x 的条件分布采样
    
    x_samples[i] <- x
    y_samples[i] <- y
    
    xx_samples = x_samples[(burn+1):N] # 取1001-5000样本
    yy_samples = y_samples[(burn+1):N]
  }
  
  return(list(xx_samples = xx_samples, yy_samples = yy_samples))
}

# 2. 进行统计推断
statistical_analysis <- function(samples) {
  xx = samples$xx_samples
  yy = samples$yy_samples
  x_mean <- mean(xx)
  y_mean <- mean(yy)
  x_var <- var(xx)
  y_var <- var(yy)
  
  return(list(x_mean = x_mean, y_mean = y_mean, x_var = x_var, y_var = y_var))
}

# 3. 输出结果，包括轨迹图和直方图
plot_results <- function(samples) {
  par(mar = c(5, 5, 2, 2))
  
  # x 和 y 的轨迹图
  plot(samples$xx_samples, type = "l", main = "Trace of x", xlab = "Iteration", ylab = "x")
  plot(samples$yy_samples, type = "l", main = "Trace of y", xlab = "Iteration", ylab = "y")
  
  # x 和 y 的直方图
  hist(samples$xx_samples, breaks = seq(-0.5, n + 0.5, by = 1), main = "Histogram of x", xlab = "x", probability = TRUE)
  hist(samples$yy_samples, breaks = 30, main = "Histogram of y", xlab = "y", probability = TRUE)
}

# 主程序
samples <- generate_samples(n, a, b, N)        # 生成样本
analysis_results <- statistical_analysis(samples)        # 统计推断
plot_results(samples)                                    # 输出结果

# 打印统计推断的结果
print(analysis_results)


## ----fig.width=7, fig.height=5------------------------------------------------
# 加载coda包
if (!require(coda)) install.packages("coda")
library(coda)

# 设置参数
n <- 10
a <- 2
b <- 2
iterations <- 10000  # 每条链的最大迭代次数
chains <- 4          # 使用的链数量

# 1. 生成多条独立的Gibbs采样链
generate_samples <- function(n, a, b, iterations) {
  x_samples <- numeric(iterations)
  y_samples <- numeric(iterations)
  y <- 0.5  # 初始化 y 的初始值
  
  for (i in 1:iterations) {
    x <- rbinom(1, n, y)              # 从 x | y 的条件分布采样
    y <- rbeta(1, x + a, n - x + b)   # 从 y | x 的条件分布采样
    
    x_samples[i] <- x
    y_samples[i] <- y
  }
  
  return(list(x_samples = x_samples, y_samples = y_samples))
}

# 生成多条链
chains_x <- lapply(1:chains, function(i) generate_samples(n, a, b, iterations)$x_samples)
chains_y <- lapply(1:chains, function(i) generate_samples(n, a, b, iterations)$y_samples)

# 将链转换为mcmc.list对象
mcmc_chains_x <- mcmc.list(lapply(chains_x, mcmc))
mcmc_chains_y <- mcmc.list(lapply(chains_y, mcmc))

# 2. 使用Gelman-Rubin方法监测收敛性
check_convergence <- function(mcmc_chains, variable_name) {
  gelman_diag <- gelman.diag(mcmc_chains, autoburnin = FALSE)
  cat(paste("Gelman-Rubin statistic for", variable_name, ":", gelman_diag$psrf[1,1], "\n"))
  return(gelman_diag$psrf[1,1] < 1.2)  # 返回是否收敛
}

# 3. 检查链的收敛性
converged_x <- check_convergence(mcmc_chains_x, "x")
converged_y <- check_convergence(mcmc_chains_y, "y")

# 输出收敛结果
if (converged_x && converged_y) {
  cat("Both chains have converged (Gelman-Rubin < 1.2)\n")
} else {
  cat("Chains have not yet converged. Consider running more iterations.\n")
}

# 4. 绘制链的轨迹图
plot_trace <- function(chains, variable_name) {
  par(mar = c(5, 5, 2, 2))  # 设置为2x2图形布局
  for (i in 1:length(chains)) {
    plot(chains[[i]], type = "l", main = paste("Trace of", variable_name, "Chain", i),
         xlab = "Iteration", ylab = variable_name)
  }
  par(mfrow = c(1, 1))  # 重置图形布局
}

# 绘制 x 和 y 的轨迹图
plot_trace(chains_x, "x")
plot_trace(chains_y, "y")



## ----echo=FALSE---------------------------------------------------------------
# 加载需要的库
library(jpeg)
library(grid) 

# 读取图片文件
img <- readJPEG("Exercise-Proof.jpg")

grid.newpage()

# 插入图片
grid.raster(img)


## -----------------------------------------------------------------------------
g = function(a,k){  # 计算每一项
  d = length(a)
  gg = (-1)^k * exp((2*k+2)*log(norm(a, type = '2')) + lgamma((d+1)/2) + lgamma(k+1.5) - lgamma(k+1) - k*log(2) - log(2*k+1) - log(2*k+2) - lgamma(k+d/2+1))
  return(gg)
}

## -----------------------------------------------------------------------------
gsum = function(n){  # k从0到n项求和
  sg = 0
  for(i in 0:n) sg = sg + g(a,i)
  return(sg)
}

## -----------------------------------------------------------------------------
N=1e3
a = matrix(c(1,2))
gsum(N)

## -----------------------------------------------------------------------------
## 定义若干生成函数
integral = function(u, k){ #积分函数
  (1+u^2/(k-1))^(-k/2)
}

getc = function(a,k){ # 计算a
  
  sqrt(a^2*k/(k+1-a^2))
   
}

## 定义表达式函数
expr = function(a,k){
  
  iintegral = function(u){
    integral(u,k)
  }
  c = getc(a,k-1)
  
  res = 2/sqrt(pi*(k-1)) * exp(lgamma(k/2) - lgamma((k-1)/2)) * integrate(iintegral,lower = 0, upper = c)$value
  return(res)
}

g = function(a){
    expr(a,k) - expr(a,k+1)
  }
## 求解方程并输出结果

resl = function(k){
  e = 0.01
 if (g(e) < 0 && g(sqrt(k) - e) > 0 || g(e) > 0 && g(sqrt(k) - e) < 0) {
    r = uniroot(g, c(e,sqrt(k)-e))$root
    } else {
      r = NA
    }
  return(r)
  
}

for (k in c(4:25,100,500,1000)){
  print(resl(k))
}

## -----------------------------------------------------------------------------
## 定义生成函数
expr2 = function(a,k){
  1-pt(sqrt(a^2*k/(k+1-a^2)), df = k)
}

g2 = function(a){  # 等式函数
  expr2(a,k) - expr2(a,k-1)
}

reslt = function(k){ # 解方程
  e = 0.01
  r = uniroot(g2, interval = c(e,sqrt(k) - e))$root
  return(r)
}
 
 # 输出结果对比
  for(k in c(4:25,100,500,1000)){
    
    p = round(c('integral solution' = resl(k),'probability solution' = reslt(k)), 3)
    print(p)
}


## -----------------------------------------------------------------------------
yo = c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
n = length(yo)

## 定义生成函数
lambda = function(){ # 原始观测数据的lambda估计
  1/mean(yo)
}

## 统计推断函数: E-M算法
em = function(tau){
  # E算法
  lambda0 = lambda() 
  yad = yo
  for(i in 1:1000){ # 迭代1000次
    # 对于修改的数据，其期望E[T|T>tau] = tao + 1/lambda
    expectedgv = tau + 1/lambda0
    # 将值为1的数据改为期望值
    yad[yo == tau] = expectedgv
    # 更新lambda
    lambda0 = 1/mean(yad)
  }
  return(lambda0)
}
## 结果输出函数
reslt3 = function(){
  round(rbind(c('E-M Algorithm' = em(1), 'Observed Data MLE' = lambda())),3)
}
reslt3()

## -----------------------------------------------------------------------------
#install.packages('nloptr')
library(nloptr) 
eval_f = function(x){
  return('objective'=4*x[1]+2*x[2]+9*x[3])
}
eval_g = function(x){
  return(list('constraints'=c(2*x[1]+x[2]+x[3]-2,x[1]-x[2]+3*x[3]-3),
              'jacobian'=rbind(c(2,1,1),c(1,-1,3))))
}
res2=nloptr(x0 = c(1,1,1),
            eval_f = eval_f,
            lb = c(0,0,0),
            ub = c(Inf,Inf,Inf),
            eval_g_ineq = eval_g,
             opts = list('algorithm'='NLOPT_LN_COBYLA',"xtol_rel"=1.0e-6))
res2$solution
res2$objective


## -----------------------------------------------------------------------------
data(mtcars)

formulas <- list(
mtcars$mpg ~ mtcars$disp,
mtcars$mpg ~ I(1 / mtcars$disp),
mtcars$mpg ~ mtcars$disp + mtcars$wt,
mtcars$mpg ~ I(1 / mtcars$disp) + mtcars$wt
)


## 使用for循环

loop_res = list()
for(i in 1:4){
  loop_res[[i]] = lm(formulas[[i]])
  
}
cat("循环输出结果为：",'\n')
print(loop_res)

## 使用lapply函数
cat("lapply函数输出结果为：",'\n')
lapply(formulas, function(f) lm(f, data = mtcars))


## -----------------------------------------------------------------------------
data(mtcars)
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})

## 使用for循环
loop_res2 = list()

for(i in 1:10){
  loop_res2[[i]] = lm(mtcars$mpg~mtcars$disp,data=bootstraps[[i]])
  
}
cat("循环输出结果为：",'\n')
print(loop_res2)

## 使用lapply函数，并且不使用匿名函数

a = function(data){
  lm(mtcars$mpg~mtcars$disp, data = data)
}

cat("lapply函数输出结果为：",'\n')
lapply(bootstraps, a)


## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
res_rsquare = function(){
for(i in 1:4) {
  r1 = rsq(lm(formulas[[i]]))
  p1 = formulas[[i]]
  print(c(formulas[[i]],r1))
}
for(j in 1:10) {
  r2 = rsq(lm(bootstraps[[j]]))
  p2 = paste('bootstraps', j,':')
  print(c(p2, r2))
}
}
res_rsquare()

## -----------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)

## 使用slapply函数：
cat("使用sapply函数：",'\n')
sapply(trials, function(x) x$p.value)

## extra challenge：
p = numeric(100)
cat("直接调用p值：",'\n')
for (i in 1:100){
  p[i] = trials[[i]]$p.value
}
print(p)


## -----------------------------------------------------------------------------
x1 = 1:4
x2 = runif(4)

f = function(x,y){
  x/y
}

res_it = simplify2array(Map(f, x1, x2))

vapply(res_it, quantile, c(1,2,5,6,8))


## -----------------------------------------------------------------------------

## 定义统计推断函数（计算chi-square）
infer_chi = function(inputs){
  # 获取inputs矩阵的行和与列和
  nrow = nrow(inputs) # 行和
  ncol = ncol(inputs) # 列和
  
  #获取总样本量
  N = sum(inputs)
  
  #获取每行与每列之和
  sumrow = rowSums(inputs)
  sumcol = colSums(inputs)
  
  #计算期望函数
  exp = outer(sumrow,sumcol)/N
  
  #计算卡方统计量
  chi_square = (inputs-exp)^2/exp
  return(sum(chi_square))
  
}

## 输出结果
inputs = matrix(c(2,3,15,10,20,10), nrow = 2, byrow = TRUE)
infer_chi(inputs)

## -----------------------------------------------------------------------------
## 定义table函数
tablef = function(x,y){
  # 获取输入向量的最大值和范围
  max_x <- max(x)
  max_y <- max(y)
  
  # 初始化一个计数矩阵
  count_matrix <- matrix(0, nrow = max_x, ncol = max_y)
  
  # 使用矩阵索引计数
  for (i in max(seq_along(x),seq_along(y)) ){
    count_matrix[x[i], y[i]] <- count_matrix[x[i], y[i]] + 1
  }
  
  return(count_matrix)

}

## 结果输出函数
res_table = function(x,y){
  result <- tablef(x, y)
  return(result)
}

x = c(1,3,6,9,12)
y = c(2,4,5,7,15)
cat("重新定义的table函数：",'\n')
res_table(x,y)


## 加速卡方检验

fast_chi = function(observed){
  inputs = tablef(observed[1], observed[2])
  
  # 获取inputs矩阵的行和与列和
  nrow = nrow(inputs) # 行和
  ncol = ncol(inputs) # 列和
  
  #获取总样本量
  N = sum(inputs)
  
  #获取每行与每列之和
  sumrow = rowSums(inputs)
  sumcol = colSums(inputs)
  
  #计算期望函数
  exp = outer(sumrow,sumcol)/N
  
  #计算卡方统计量
  chi_square = (inputs-exp)^2/exp
  return(sum(chi_square))
}

observed = matrix(c(1,4,3,2,5,5,3,2,4,1,5,3), nrow = 2, byrow = TRUE)

cat("加速卡方检验的统计量：",'\n')
fast_chi(observed)


## -----------------------------------------------------------------------------
tablef <- function(numbers) {
  # 创建一个空的列表来存储数字及其对应的频数
freq_list <- list()
 
# 遍历输入的数字
for (num in numbers) {
  # 如果数字已经在列表中，增加其频数
  if (num %in% names(freq_list)) {
    freq_list[[as.character(num)]] <- freq_list[[as.character(num)]] + 1
  } else {
    # 如果数字不在列表中，添加它并设置频数为1
    freq_list[[as.character(num)]] <- 1
  }
}
 
# 将列表转换为两个向量：唯一的数字和对应的频数
unique_numbers <- as.numeric(names(freq_list))
frequencies <- as.numeric(unlist(freq_list))
 
# 构建矩阵
result_matrix <- matrix(c(unique_numbers, frequencies), nrow = 2, byrow = TRUE)
 
 
return(result_matrix)

}

# 示例数据
x <- c(1,3,1,1,4,2,2,4,3,5)


# 调用快速 table 函数
result <- tablef(x)
cat("重新定义的table函数：",'\n')
print(result)


## -----------------------------------------------------------------------------
fast_chi = function(numbers){
  
  inputs = tablef(numbers)
  
  # 获取inputs矩阵的行和与列和
  nrow = nrow(inputs) # 行和
  ncol = ncol(inputs) # 列和
  
  #获取总样本量
  N = sum(inputs)
  
  #获取每行与每列之和
  sumrow = rowSums(inputs)
  sumcol = colSums(inputs)
  
  #计算期望函数
  exp = outer(sumrow,sumcol)/N
  
  #计算卡方统计量
  chi_square = (inputs-exp)^2/exp
  return(sum(chi_square))
}

observed = matrix(c(1,3,1,1,4,2,2,4,3,5), nrow = 2, byrow = TRUE)

cat("加速卡方检验的统计量：",'\n')
fast_chi(observed)

