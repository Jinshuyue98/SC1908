## ------------------------------------------------------------------------
library(xtable)
x<-c(1,2,3,4)
y<-c(1,2,3,4)
x1<-matrix(c(1,2,3,4,1,2,3,4),nrow=2,ncol=4,byrow=TRUE,dimnames=list(c("x","y"),c("1","2","3","4")))
print(xtable(x1),type="html")

x<-c(1,2,3,4)
y<-c(1,2,3,4)
plot(x,y)

x<-c(1,2,3,4)
y<-c(1,2,3,4)
x1<-matrix(c(1,2,3,4,1,2,3,4),nrow=2,ncol=4,byrow=TRUE,dimnames=list(c("x","y"),c("1","2","3","4")))
knitr::kable(x1,format="markdown")
plot(x,y,type="b",col="blue")

knitr::kable(head(cars),format="markdown")
plot(head(cars),type="p",col="red",main="cars")

## ------------------------------------------------------------------------
n<-1000
sigma<-0.5#给σ赋值
u<-runif(n)#F(x)服从(0,1)上的均匀分布
x<-sqrt(-2*sigma^2*log(1-u,base=exp(1)))#对x取样
hist(x,prob=TRUE,breaks=25,main=expression(f(x)==(x/sigma^2)*e^(-x^2/2*sigma^2)))#画x的直方图
y<-seq(0,2,0.01)
lines(y,(y/sigma^2)*exp(-y^2/(2*sigma^2)))#画密度函数图像

## ------------------------------------------------------------------------
n<-1000
sigma<-1#给σ赋值
u<-runif(n)#F(x)服从(0,1)上的均匀分布
x<-sqrt(-2*sigma^2*log(1-u,base=exp(1)))#对x取样
hist(x,prob=TRUE,breaks=25,main=expression(f(x)==(x/sigma^2)*e^(-x^2/2*sigma^2)))#画x的直方图
y<-seq(0,4,0.01)
lines(y,(y/sigma^2)*exp(-y^2/(2*sigma^2)))#画密度函数图像

## ------------------------------------------------------------------------
n<-1000
sigma<-1.2#给σ赋值
u<-runif(n)#F(x)服从(0,1)上的均匀分布
x<-sqrt(-2*sigma^2*log(1-u,base=exp(1)))#对x取样
hist(x,prob=TRUE,breaks=25,main=expression(f(x)==(x/sigma^2)*e^(-x^2/2*sigma^2)))#画x的直方图
y<-seq(0,6,0.01)
lines(y,(y/sigma^2)*exp(-y^2/(2*sigma^2)))#画密度函数图像

## ------------------------------------------------------------------------
n<-1000#样本容量
p<-c(0.25,0.75)#p1=0.75
x1<-rnorm(n,mean=0,sd=1)#对x1，x2取样
x2<-rnorm(n,mean=3,sd=1)
r<-sample(c(0,1),n,prob=p,replace=TRUE)#对混合概率取样
z<-r*x1+(1-r)*x2#混合分布
hist(z,prob=TRUE,breaks=50,main=expression(0.75*x[1]+0.25*x[2]))#画直方图
lines(density(z))#画密度图像

## ------------------------------------------------------------------------
n<-1000
p<-c(0.3,0.7)#p1=0.7
x1<-rnorm(n,mean=0,sd=1)
x2<-rnorm(n,mean=3,sd=1)
r<-sample(c(0,1),n,prob=p,replace=TRUE)
z<-r*x1+(1-r)*x2
hist(z,prob=TRUE,breaks=50,main=expression(0.7*x[1]+0.3*x[2]))
lines(density(z))

## ------------------------------------------------------------------------
n<-1000
p<-c(0.4,0.6)#p1=0.6
x1<-rnorm(n,mean=0,sd=1)
x2<-rnorm(n,mean=3,sd=1)
r<-sample(c(0,1),n,prob=p,replace=TRUE)
z<-r*x1+(1-r)*x2
hist(z,prob=TRUE,breaks=50,main=expression(0.6*x[1]+0.4*x[2]))
lines(density(z))

## ------------------------------------------------------------------------
n<-1000
p<-c(0.5,0.5)#p1=0.5
x1<-rnorm(n,mean=0,sd=1)
x2<-rnorm(n,mean=3,sd=1)
r<-sample(c(0,1),n,prob=p,replace=TRUE)
z<-r*x1+(1-r)*x2
hist(z,prob=TRUE,breaks=50,main=expression(0.5*x[1]+0.5*x[2]))
lines(density(z))

## ------------------------------------------------------------------------
n<-1000
p<-c(0.6,0.4)#p1=0.4
x1<-rnorm(n,mean=0,sd=1)
x2<-rnorm(n,mean=3,sd=1)
r<-sample(c(0,1),n,prob=p,replace=TRUE)
z<-r*x1+(1-r)*x2
hist(z,prob=TRUE,breaks=50,main=expression(0.4*x[1]+0.6*x[2]))
lines(density(z))

## ------------------------------------------------------------------------
sample.Wishart<-function(sigma,nf,d,N){
  #根据Bartlett分解对Wishart分布进行取样的函数
  #Σ是Wishart分布的协方差矩阵
  #nf是Wishart分布的自由度
  #d是正态分布的维度
  #N是样本容量
  L<-chol(sigma)# 计算Σ的Cholesky矩阵
  n<-1
  X<-as.list(numeric(N))
  for(n in 1:N){
    #对Wishart分布取样N次
    i<-1
    A<-matrix(numeric(d^2),nrow=d)
    c2<-numeric(d)
    for(i in 1:d){
      #对A矩阵对角元取样
      c2[i]<-rchisq(n=1,df=nf-i+1)
      A[i,i]<-sqrt(c2[i])
    }
    A
    i<-1
    j<-1
    for(i in 1:d){
      #对A其他元素取样
      for(j in 1:i){
      A[i,j]<rnorm(n=1)
      }
    }
    X[[n]]<-L%*%A%*%t(L)%*%t(A)#得到一个Wishart分布样本
  }
  return(X)#返回取样结果
}
#试运行Wishart取样函数
si<-matrix(c(1,2,3,2,1,2,3,2,1),nrow=3,byrow=TRUE)#任取3X3矩阵
sig<-si%*%t(si)#保证Σ是对称正定矩阵
#取自由度为5，维度为3，样本容量为5
n<-5
d<-3
N<-5
X<-sample.Wishart(sigma=sig,nf=n,d=d,N=N)#用Wishart取样函数取样
X

## ------------------------------------------------------------------------
M<-10000#总样本容量
k<-10#分层数
r<-M/k#每层取样数
N<-50#重复估计的次数
T<-numeric(k)
est<-numeric(N)
g<-function(x)sin(x)*(x>0)*(x<(pi/3))
for(i in 1:N){
  for(j in 1:k)T[j]<-mean(g(runif(r,(j-1)*pi/(k*3),j*pi/(k*3))))
  est[i]<-mean(T)
}
mean<-mean(est)
s<-sd(est)
mean
s

## ------------------------------------------------------------------------
M<-100000#总样本容量
k<-100#分层数
r<-M/k#每层取样数
N<-50#重复估计的次数
T<-numeric(k)
est<-numeric(N)
g<-function(x)sin(x)*(x>0)*(x<(pi/3))
for(i in 1:N){
  for(j in 1:k)T[j]<-mean(g(runif(r,(j-1)*pi/(k*3),j*pi/(k*3))))
  est[i]<-mean(T)
}
mean<-mean(est)
s<-sd(est)
mean
s

## ------------------------------------------------------------------------
k<-100#普通蒙特卡洛估计的样本容量
N<-50#重复估计的次数
T1<-numeric(k)
T2<-numeric(k/2)
est<-matrix(0,N,2)
g1<-function(x)exp(-x)/(1+x^2)*(x>0)*(x<1)#普通蒙特卡洛估计
g2<-function(x)(exp(-x)/(1+x^2)+exp(-(1-x))/(1+(1-x)^2))*(x>0)*(x<1)#对偶变量蒙特卡洛估计
for(i in 1:N){
  for(j in 1:k)T1[j]<-g1(runif(1,0,1))
  est[i,1]<-mean(T1)
  for(a in 1:(k/2))T2[a]<-g2(runif(1,0,1))
  est[i,2]<-sum(T2)/k
}
c(sd(est[,1]),sd(est[,2]),sd(est[,2])/sd(est[,1]))

## ------------------------------------------------------------------------
p<-c(0,0.2,0.4,0.6,0.8,1)
q<-numeric(6)
for(a in 1:6){
  #求f(x)5分位点
  q[a]<--log(1-p[a]*(1-exp(-1)))
}
M<-100#总样本容量
k<-5#分层数
r<-M/k#每层取样数
N<-50#重复估计的次数
T<-numeric(k)
est<-numeric(N)
for(i in 1:N){
  for(j in 1:5){
    u<-runif(r)#f3, inverse transform method
    x<--log(exp(-q[j])-u*((1-exp(-1))/5))
    gq<-function(x)exp(-x-log(1+x^2))
    fq<-function(x)(exp(-x)/(1-exp(-1)))
    fg<-gq(x)/(5*fq(x))
    T[j]<-mean(fg)#计算每个区间的积分估计值
  }
  est[i]<-sum(T)#对所有区间求和，得到整个区间的积分估计值
}
mean(est)
sd(est)

## ------------------------------------------------------------------------
m<-1000#Monte Carlo estimation times
n<-20#sample size
up<-numeric(m)
down<-numeric(m)
mean<-numeric(m)
s<-numeric(m)
t<-qt(0.975,df=(n-1))
for(i in 1:m){
  x <- rchisq(n, df=2)
  mean[i]<-mean(x)
  s[i]<-var(x)
  up[i]<-(mean[i]+t*s[i]/sqrt(n))#up line of CI
  down[i]<-(mean[i]-t*s[i]/sqrt(n))#down line of CI
}
p<-mean(up>2&down<2)#compute the mean to get the confidence level
p
se<-sqrt((1-p)*p/m)#The standard error of the estimate
se

## ------------------------------------------------------------------------
m<-100#Monte Carlo estimation times
n<-200#sample size
N<-100#Number of repeated sampling
skew<-numeric(N)
p<-c(0.025,0.05,0.95,0.975)
q<-matrix(m*4,nrow=m,ncol=4)
for(i in 1:m){
  for(j in 1:N){
    x<-rnorm(n,mean=0,sd=1)#generate n samples from normal distribution
    skew[j]<-mean((x-mean(x))*3/(var(x))*(3/2))#compute skew of n samples
  }
  q[i,]<-quantile(skew,probs=p)#compute quantiles of N skew
}
quantile<-apply(q,2,mean)#the estimates of quantiles
quantile
ql<-qnorm(p,mean=0,sd=sqrt(6/100000))#the quantiles of large sample approximation
ql
f<-dnorm(quantile,mean=0,sd=sqrt(6/n))#the density of the normal distribution
v<-numeric(4)
for(a in 1:4){
  v[a]<-p[a]*(1-p[a])/(m*(f[a])^2)#The standard error of the estimate
}
v

## ------------------------------------------------------------------------
set.seed(1234)
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
alpha <- .05
n <- 30
m <- 2500
beta <- c(seq(0, 10, .5), seq(10, 100, 2))
N <- length(beta)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
  b <- beta[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    x <- rbeta(n, b, b)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs beta
plot(beta, pwr, type = "b",pch=19,cex=0.5,xlab = bquote(beta), ylim = c(0,1))
abline(h = alpha, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(beta, pwr+se, lty = 3)
lines(beta, pwr-se, lty = 3)

## ------------------------------------------------------------------------
set.seed(1234)
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
alpha <- .05
n <- 30
m <- 2500
degree <- c(seq(1, 100, 1))
N <- length(degree)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon
  df <- degree[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    x <- rt(n,df)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs degree
plot(degree, pwr, type = "b",pch=19,cex=0.5,xlab = bquote(degree), ylim = c(0,1))
abline(h = alpha, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(degree, pwr+se, lty = 3)
lines(degree, pwr-se, lty = 3)

## ------------------------------------------------------------------------
set.seed(1234)
alpha <- 0.05 
miu <- 1 # Null hypothesis!
m <- 1e4
n <- 60
set.seed(123)
p.val1 <- numeric(m)
for(i in 1:m){
x<-rchisq(n,1)
hat<-sqrt(n)*(mean(x)-miu)
se<-var(x)
p.val1[i] <- 2*(1-pt(abs(hat/se),n-1)) # t-test p-value
}
print(c(mean(p.val1<=alpha),alpha))

## ------------------------------------------------------------------------
set.seed(1234)
alpha <- 0.05 
miu <- 1 # Null hypothesis!
m <- 1e4
n <- 60
set.seed(123)
p.val1 <- numeric(m)
for(i in 1:m){
x<-runif(n,min=0,max=2)
hat<-sqrt(n)*(mean(x)-miu)
se<-var(x)
p.val1[i] <- 2*(1-pt(abs(hat/se),n-1)) # t-test p-value
}
print(c(mean(p.val1<=alpha),alpha))

## ------------------------------------------------------------------------
set.seed(1234)
alpha <- 0.05 
miu <- 1 # Null hypothesis!
m <- 1e4
n <- 60
set.seed(123)
p.val1 <- numeric(m)
for(i in 1:m){
x<-rexp(n,rate=1)
hat<-sqrt(n)*(mean(x)-miu)
se<-var(x)
p.val1[i] <- 2*(1-pt(abs(hat/se),n-1)) # t-test p-value
}
print(c(mean(p.val1<=alpha),alpha))

## ------------------------------------------------------------------------
set.seed(1234)
library(bootstrap)
score<-as.matrix(scor)
#display the scatter plots for each pair of test scores
for(i in 1:5){
  for(j in 1:5){
    if(j>i) plot(score[,i],score[,j],xlab=colnames(score)[i],ylab=colnames(score)[j])
  }
}
#compute the sample correlation matrix
cor<-as.matrix(cor(score))
cor
#function to obtain bootstrap estimates of correlation
cor.b<-function(x,y,B){
  n<-length(x)
  corstar<-numeric(B)
  for(i in 1:B){
    xstar<-sample(x,n,replace=TRUE)
    ystar<-sample(y,n,replace=TRUE)
    corstar[i]<-cor(xstar,ystar)
  }
  return(corstar)
}
B<-100
cor12.b<-cor.b(score[,1],score[,2],B)
cor34.b<-cor.b(score[,3],score[,4],B)
cor35.b<-cor.b(score[,3],score[,5],B)
cor45.b<-cor.b(score[,4],score[,5],B)
SE<-matrix(c(sd(cor12.b),sd(cor34.b),sd(cor35.b),sd(cor45.b)),nrow=1,ncol=4,byrow=TRUE,dimnames=list(c("SE"),c("SE12.b","SE34.b","SE35.b","SE45.b")))
SE

## ------------------------------------------------------------------------
set.seed(1234)
#computes the sample skewness coeff.
sk<-function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
#normal populations
alpha<-0.05
M<-100
B<-100
n<-2000
sk.star<-numeric(B)
se<-up<-down<-numeric(M)

#estimate the coverage probabilities of the standard normal bootstrap confidence interval
for(i in 1:M){
  xhat<-rnorm(n,mean=0,sd=1)
  sk.hat<-sk(xhat)
  for(j in 1:B){
    xstar<-sample(xhat,n,replace=TRUE)
    sk.star[j]<-sk(xstar)
  }
  se[i]<-sd(sk.star)
  up[i]<-(sk.hat-qnorm(alpha/2,mean=0,sd=1)*se[i])
  down[i]<-(sk.hat-qnorm(1-alpha/2,mean=0,sd=1)*se[i])
}
CP.s<-mean(down<0&up>0)
s.left<-mean(down>0)
s.right<-mean(up<0)
#estimate the coverage probabilities of the basic bootstrap confidence interval
for(i in 1:M){
  xhat<-rnorm(n,mean=0,sd=1)
  sk.hat<-sk(xhat)
  for(j in 1:B){
    xstar<-sample(xhat,n,replace=TRUE)
    sk.star[j]<-sk(xstar)
  }
  up[i]<-(2*sk.hat-quantile(sk.star,probs=alpha/2))
  down[i]<-(2*sk.hat-quantile(sk.star,probs=(1-alpha/2)))
}
CP.b<-mean(down<0&up>0)
b.left<-mean(down>0)
b.right<-mean(up<0)
#estimate the coverage probabilities of the percentile confidence interval
for(i in 1:M){
  xhat<-rnorm(n,mean=0,sd=1)
  sk.hat<-sk(xhat)
  for(j in 1:B){
    xstar<-sample(xhat,n,replace=TRUE)
    sk.star[j]<-sk(xstar)
  }
  up[i]<-quantile(sk.star,probs=(1-alpha/2))
  down[i]<-quantile(sk.star,probs=alpha/2)
}
CP.p<-mean(down<0&up>0)
p.left<-mean(down>0)
p.right<-mean(up<0)

norm<-matrix(c(CP.s,s.left,s.right,CP.b,b.left,b.right,CP.p,p.left,p.right),nrow=1,ncol=9)

#χ2(5) distributions
sk1<-sqrt(8/5)
sk.star<-numeric(B)
se<-up<-down<-numeric(M)

#estimate the coverage probabilities of the standard normal bootstrap confidence interval
for(i in 1:M){
  xhat<-rchisq(n,df=5)
  sk.hat<-sk(xhat)
  for(j in 1:B){
    xstar<-sample(xhat,n,replace=TRUE)
    sk.star[j]<-sk(xstar)
  }
  se[i]<-sd(sk.star)
  up[i]<-(sk.hat-qnorm(alpha/2,mean=0,sd=1)*se[i])
  down[i]<-(sk.hat-qnorm(1-alpha/2,mean=0,sd=1)*se[i])
}
CP.s<-mean(down<sk1&up>sk1)
s.left<-mean(down>sk1)
s.right<-mean(up<sk1)

#estimate the coverage probabilities of the basic bootstrap confidence interval
for(i in 1:M){
  xhat<-rchisq(n,df=5)
  sk.hat<-sk(xhat)
  for(j in 1:B){
    xstar<-sample(xhat,n,replace=TRUE)
    sk.star[j]<-sk(xstar)
  }
  up[i]<-(2*sk.hat-quantile(sk.star,probs=alpha/2))
  down[i]<-(2*sk.hat-quantile(sk.star,probs=(1-alpha/2)))
}
CP.b<-mean(down<sk1&up>sk1)
b.left<-mean(down>sk1)
b.right<-mean(up<sk1)

#estimate the coverage probabilities of the percentile confidence interval
for(i in 1:M){
  xhat<-rchisq(n,df=5)
  sk.hat<-sk(xhat)
  for(j in 1:B){
    xstar<-sample(xhat,n,replace=TRUE)
    sk.star[j]<-sk(xstar)
  }
  up[i]<-quantile(sk.star,probs=(1-alpha/2))
  down[i]<-quantile(sk.star,probs=alpha/2)
}
CP.p<-mean(down<sk1&up>sk1)
p.left<-mean(down>sk1)
p.right<-mean(up<sk1)

chisq<-matrix(c(CP.s,s.left,s.right,CP.b,b.left,b.right,CP.p,p.left,p.right),nrow=1,ncol=9)
compare<-matrix(c(norm,chisq),nrow=2,ncol=9,byrow=TRUE,dimnames=list(c("norm","chisq"),c("CP.standard","standard.left","standard.right","CP.basic","basic.left","basic.right","CP.percentile","percentile.left","percentile.right")))
compare


## ------------------------------------------------------------------------
set.seed(1234)
library(bootstrap)
scor<-as.matrix(scor)
lamda<-eigen(cor(scor))$values
theta.hat<-lamda[1]/sum(lamda)
J<-nrow(scor)
theta.jack<-numeric(J)
for(i in 1:J){
  #obtain the jackknife estimates of theta
  scor.hat<-scor[-i,]
  n<-nrow(scor.hat)
  sigma.MLE<-((n-1)/n)*cor(scor.hat)
  lamda.hat<-eigen(sigma.MLE)$values
  theta.jack[i]<-lamda.hat[1]/sum(lamda.hat)
}
#Obtain the jackknife estimates of bias and standard error
bias.jack<-(J-1)*(mean(theta.jack)-theta.hat)
sd.jack<-sqrt((J-1)*mean((theta.jack-theta.hat)^2))
round(c(bias.jack=bias.jack,sd.jack=sd.jack),2)

## ------------------------------------------------------------------------
set.seed(1234)
library(DAAG)
attach(ironslag)
ybar<-mean(magnetic)
n<-length(magnetic) #in DAAG ironslag
e1<-e2<-e3<-e4<-numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for(k in 1:n){
  y<-magnetic[-k]
  x<-chemical[-k]

  J1<-lm(y~x)
  yhat1<-J1$coef[1]+J1$coef[2]*chemical[k]
  e1[k]<-magnetic[k]-yhat1

  J2<-lm(y~x+I(x^2))
  yhat2<-J2$coef[1]+J2$coef[2]*chemical[k]+J2$coef[3]*chemical[k]^2
  e2[k]<-magnetic[k]-yhat2

  J3<-lm(log(y)~x)
  logyhat3<-J3$coef[1]+J3$coef[2]*chemical[k]
  yhat3<-exp(logyhat3)
  e3[k]<-magnetic[k]-yhat3

  J4<-lm(y~x+I(x^2)+I(x^3))
  yhat4<-J4$coef[1]+J4$coef[2]*chemical[k]+J4$coef[3]*chemical[k]^2+J4$coef[4]*chemical[k]^3
  e4[k]<-magnetic[k]-yhat4
}
round(c(MSE1=mean(e1^2),MSE2=mean(e2^2),MSE3=mean(e3^2),MSE4=mean(e4^2)),4)
#adjusted R^2
SSE1<-SSE2<-SSE3<-SSE4<-numeric(n)
SST<-sum((magnetic-ybar)^2)

M1<-lm(magnetic~chemical)
yhat1<-M1$coef[1]+M1$coef[2]*chemical
SSE1<-sum((yhat1-magnetic)^2)
R1<-(1-SSE1*(n-1)/(SST*(n-1-1)))

M2<-lm(magnetic~chemical+I(chemical^2))
yhat2<-M2$coef[1]+M2$coef[2]*chemical+M2$coef[3]*chemical^2
SSE2<-sum((yhat2-magnetic)^2)
R2<-(1-SSE2*(n-1)/(SST*(n-2-1)))

M3<-lm(log(magnetic)~chemical)
logyhat3<-M3$coef[1]+M3$coef[2]*chemical
yhat3<-exp(logyhat3)
SSE3<-sum((yhat3-magnetic)^2)
R3<-(1-SSE3*(n-1)/(SST*(n-1-1)))

M4<-lm(magnetic~chemical+I(chemical^2)+I(chemical^3))
yhat4<-M4$coef[1]+M4$coef[2]*chemical+M4$coef[3]*chemical^2+M4$coef[4]*chemical^3
SSE4<-sum((yhat4-magnetic)^2)
R4<-(1-SSE4*(n-1)/(SST*(n-3-1)))

round(c(adjusted.R1=R1,adjusted.R2=R2,adjusted.R3=R3,adjusted.R4=R4),4)

## ------------------------------------------------------------------------
set.seed(1234)
count5<-function(x1,x2){
  R<-1000
  T1<-T2<-numeric(R)
  z<-c(x1,x2)
  n<-length(z)
  K<-1:n
  n1<-length(x1)
  for(i in 1:R){
    k<-sample(K,n1,replace=FALSE)
    x1.star<-z[k]
    x2.star<-z[-k]
    T1[i]<-sum(x1.star>max(x1))+sum(x1.star<min(x1))
    T2[i]<-sum(x2.star>max(x2))+sum(x2.star<min(x2))
  }
  p1<-(1+sum(T1>=5))/(1+R)
  p2<-(1+sum(T2>=5))/(1+R)
  return(round(c(p1=p1,p2=p2),3))
}
x1<-rnorm(10,0,3)
x2<-rnorm(20,0,3)
alpha<-0.05
p<-max(count5(x1,x2))
if(p>=alpha) sprintf("p-value is %.3f.At the %.2f confidence level, the variances of two samples are different",p,alpha)
if(p<alpha) sprintf("p-value is %.3f.At the %.2f confidence level, the variances of two samples are equal",p,alpha)
x1<-rnorm(10,0,1)
x2<-rnorm(20,0,3)
alpha<-0.05
p<-max(count5(x1,x2))
if(p>=alpha) sprintf("p-value is %.3f.At the %.2f confidence level, the variances of two samples are different",p,alpha)
if(p<alpha) sprintf("p-value is %.3f.At the %.2f confidence level, the variances of two samples are equal",p,alpha)

## ------------------------------------------------------------------------
dCov<-function(x, y){
  x<-as.matrix(x)
  y<-as.matrix(y)
  n<-nrow(x)
  m<-nrow(y)
  if(n!=m||n<2) stop("Sample sizes must agree")
  if(!(all(is.finite(c(x, y))))) stop("Data contains missing or infinite values")
  Akl<-function(x){
    d<-as.matrix(dist(x))
    m<-rowMeans(d)
    M<-mean(d)
    a<-sweep(d,1,m)
    b<-sweep(a,2,m)
    return(b+M)
  }
  A<-Akl(x)
  B<-Akl(y)
  return(sqrt(mean(A*B)))
}
ndCov2<-function(z,ix,dims){
#dims contains dimensions of x and y
  p<-dims[1]
  q<-dims[2]
  d<-p+q
  x<-z[,1:p]#leave x as is
  y<-z[ix,-(1:p)]#permute rows of y
  return(nrow(z)*dCov(x, y)^2)
}
library("boot")
library("Ball")
set.seed(12345)
m<-10
R<-100
alpha<-0.05
size<-c(10,20,30,40,50,60,70,80,90,100,200)
b<-length(size)
pow1<-pow2<-matrix(NA,b,2)
for(j in 1:b){
  p.values1<-p.values2<-matrix(NA,m,2)
  for(i in 1:m){
    x<-matrix(rnorm(size[j]*2),ncol=2)
    e<-matrix(rnorm(size[j]*2),ncol=2)
    y1<-x/4+e
    y2<-(x/4)*e
    z1<-cbind(x,y1)
    z2<-cbind(x,y2)
    boot.obj1<-boot(data=z1,statistic=ndCov2,R=R,sim="permutation",dims=c(2,2))
    #permutatin:resampling without replacement
    p.values1[i,1]<-mean(boot.obj1$t>=boot.obj1$t0)
    p.values1[i,2]<-bd.test(x=x,y=y1,R=R,seed=j*i*123)$p.value
    boot.obj2<-boot(data=z2,statistic=ndCov2,R=R,sim="permutation",dims=c(2,2))
    #permutatin:resampling without replacement
    p.values2[i,1]<-mean(boot.obj2$t>=boot.obj2$t0)
    p.values2[i,2]<-bd.test(x=x,y=y2,R=R,seed=j*i*123)$p.value
  }
  pow1[j,1]<-mean(p.values1[,1]<alpha)
  pow1[j,2]<-mean(p.values1[,2]<alpha)
  pow2[j,1]<-mean(p.values2[,1]<alpha)
  pow2[j,2]<-mean(p.values2[,2]<alpha)
}
plot(size,pow1[,1],sub="Model1",xlab="sample size",ylab="power",type="l",col=1)
lines(size,pow1[,2],type="l",col=2)
legend('topright',legend=c('distance correlation','ball covariance'),col=1:2,lty=c(1,1))
plot(size,pow2[,1],sub="Model2",xlab="sample size",ylab="power",type="l",col=1)
lines(size,pow2[,2],type="l",col=2)
legend('topright',legend=c('distance correlation','ball covariance'),col=1:2,lty=c(1,1))


## ------------------------------------------------------------------------
L.d<-function(x){
  f<-(exp(-abs(x)))/2
  return(f)
}
L.Metropolis<-function(sigma,N){
  x<-numeric(N)
  x[1]<-rnorm(1,0,sigma)
  u<-runif(N)
  r<-0
  a<-0
  for(i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    fy<-L.d(y)
    fx<-L.d(x[i-1])
    if (u[i]<=fy/fx){
      x[i]<-y
      a<-a+1
    }
    else {
      x[i]<-x[i-1]
      r<-r+1
    }
  }
  a.rate<-a/N
return(list(x=x,r=r,a.rate=a.rate))
}
library("GeneralizedHyperbolic")
set.seed(123)
N<-2000
sigma<-c(5, 10, 15, 20)
L1<-L.Metropolis(sigma[1],N)
L2<-L.Metropolis(sigma[2],N)
L3<-L.Metropolis(sigma[3],N)
L4<-L.Metropolis(sigma[4],N)
#number of candidate points rejected
round(c(reject.1=L1$r,reject.2=L2$r,reject.3=L3$r,reject.4=L4$r))
#acceptance rate
round(c(acceptance.rate1=L1$a.rate,acceptance.rate2=L2$a.rate,acceptance.rate3=L3$a.rate,acceptance.rate4=L4$a.rate),4)

refline<-qskewlap(c(0.025,0.975))
L<-cbind(L1$x,L2$x,L3$x,L4$x)
for (j in 1:4) {
  plot(L[,j],type="l",xlab=bquote(sigma==.(round(sigma[j],3))),ylab="X",ylim=range(L[,j]))
  abline(h=refline)
}
par(mfrow=c(1,1)) #reset to default

p<-c(.05,seq(.1, .9, .1),.95)
Q<-qskewlap(p)
L<-cbind(L1$x,L2$x,L3$x,L4$x)
mc<-L[501:N, ]
QL<-apply(mc,2,function(x) quantile(x, p))


Qc<-matrix(round(cbind(Q,QL),3),nrow=11,ncol=5,dimnames=list(c("5%","10%","20%","30%","40%","50%","60%","70%","80%","90%","95%"),c("origin","5","10","15","20")))
knitr::kable(Qc,format="markdown")

## ------------------------------------------------------------------------
x<-0.1
isTRUE(log(exp(x))==exp(log(x)))
isTRUE(all.equal(log(exp(x)),exp(log(x))))
identical(log(exp(x)),exp(log(x)))

## ------------------------------------------------------------------------
f<-function(k,a){
  ck<-sqrt((a^2*k)/(k+1-a^2))
  f1<-2*gamma((k+1)/2)/(sqrt(pi*k)*gamma(k/2))
  f2<-function(u){
    (1+u^2/k)^(-(k+1)/2)
  }
  inte<-integrate(f2,lower=0,upper=ck,rel.tol=.Machine$double.eps^0.25)
  return(f1*inte$value)
}
findroot<-function(k){
  a<-seq(0,sqrt(k),length.out=50*k)
  n<-length(a)
  b<-numeric(n-2)
  for(i in 2:(n-1)){
    b[i]<-abs(f(k-1,a[i])-f(k,a[i]))
  }
  if(min(b)<=0.001){
    i<-which(b==min(b),arr.ind=TRUE)
    return(a[i+1])
  }
  else return(a[1])
}
k<-c(4:25,100)
N<-length(k)
a.root<-numeric(N)
for(j in 1:N){
  a.root[j]<-findroot(k[j])
}
a.root

#points A(k) in 11.4
s<-function(k,a){
  ck<-sqrt((a^2*k)/(k+1-a^2))
  return(1-pt(ck,k))
}
sroot<-function(k){
  a<-seq(0,sqrt(k),length.out=50*k)
  n<-length(a)
  b<-numeric(n-2)
  for(i in 2:(n-1)){
    b[i]<-abs(s(k-1,a[i])-s(k,a[i]))
  }
  if(min(b)<=0.001){
    i<-which(b==min(b),arr.ind=TRUE)
    return(a[i+1])
  }
  else return(a[1])
}
s.root<-numeric(N)
j<-1
for(j in 1:N){
  s.root[j]<-sroot(k[j])
}
s.root

## ------------------------------------------------------------------------

library(nloptr)
# Mle 
eval_f0 = function(x,x1,n.A=28,n.B=24,nOO=41,nAB=70) {
  
  r1 = 1-sum(x1)
  nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r = 1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}


# constraint
eval_g0 = function(x,x1,n.A=28,n.B=24,nOO=41,nAB=70) {
  return(sum(x)-0.999999)
}

opts = list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-8)
mle = NULL
r = matrix(0,1,2)
r = rbind(r,c(0.2,0.35))# the beginning value of p0 and q0
j = 2
while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
res = nloptr( x0=c(0.3,0.25),
               eval_f=eval_f0,
               lb = c(0,0), ub = c(1,1), 
               eval_g_ineq = eval_g0, 
               opts = opts, x1=r[j,],n.A=28,n.B=24,nOO=41,nAB=70 )
j = j+1
r = rbind(r,res$solution)
mle = c(mle,eval_f0(x=r[j,],x1=r[j-1,]))
}
#the result of EM algorithm
r 
#the max likelihood values
plot(mle,type = 'l')


## ------------------------------------------------------------------------
formulas<-list(mpg~disp,mpg~I(1/disp),mpg~disp+wt,mpg~I(1/disp)+wt)
fun<-function(i){
  lm(formulas[[i]],data=mtcars)
}
lapply(1:4,fun)
lm<-list(NULL)
length(lm)<-4
for(i in 1:4){
  lm[[i]]<-fun(i)
}
lm

## ------------------------------------------------------------------------
boot<-function(x) x[sample(nrow(x),replace=TRUE),]
boot_lm<-function(i){
  dat<-boot(mtcars)
  lm(mpg~disp,data=dat)
}
lapply(1:10,boot_lm)
b<-list(NULL)
length(b)<-10
for(j in 1:10){
  b[[j]]<-boot_lm(j)
}
b

## ------------------------------------------------------------------------
rsq<-function(i) summary(fun(i))$r.squared
lapply(1:4,rsq)
boot_rsq<-function(i){
  summary(boot_lm(i))$r.squared
}
lapply(1:10,boot_rsq)

## ------------------------------------------------------------------------
trials<-replicate(100,t.test(rpois(10,10),rpois(7,10)),simplify=FALSE)
p<-function(i) trials[[i]]$p.value
sapply(1:100,p)

## ------------------------------------------------------------------------
library(parallel)
boot.rsq<-function(i){
  dat<-mtcars[sample(nrow(mtcars),replace=TRUE),]
  lm<-lm(mpg~disp,data=dat)
  summary(lm)$r.squared
}
no_cores<-7
cl<-makeCluster(no_cores)
system.time(sapply(1:100,boot.rsq))
system.time(parSapply(cl,1:100,boot.rsq))
system.time(vapply(1:100,boot.rsq,0))
stopCluster(cl)

## ------------------------------------------------------------------------

set.seed(12345)

#the function using R
lap_f = function(x) exp(-abs(x))

rw.Metropolis = function(sigma, x0, N){
 x = numeric(N)
 x[1] = x0
 u = runif(N)
 k = 0
 for (i in 2:N) {
  y = rnorm(1, x[i-1], sigma)
  if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
  else {
  x[i] = x[i-1]
  k = k+1
  }
 }
 return(list(x = x, k = k))
}

#the function using Rcpp
library(Rcpp)
cppFunction('NumericVector cpp_Metropolis(double sigma, double x0, int N){
NumericVector x(N);
x[0]=x0;
for (int i=1;i<N;i++) {
double e=runif(1)[0];
double z=rnorm(1,x[i-1],sigma)[0];
if (e<=exp(abs(x[i-1])-abs(z))) x[i]=z;
else {
x[i]=x[i-1];
}
}
return x;
}')

#generate random numbers
N<-1000
sigma<-1
x0<-0
R<-rw.Metropolis(sigma,x0,N)$x
cpp<-cpp_Metropolis(sigma,x0,N)

#qqplot
par(mfrow=c(1,2))
qqnorm(R)
qqline(R)
qqnorm(cpp,col="red")
qqline(cpp,col="red")

#computation time
library(microbenchmark)
time<-microbenchmark(R=rw.Metropolis(sigma,x0,N),Rcpp=cpp_Metropolis(sigma,x0,N))
summary(time)[,c(1,3,5,6)]

