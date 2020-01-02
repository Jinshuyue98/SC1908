#' @title Estimating BCa confidence interval with Bootstrap using R
#' @description Estimating BCa confidence interval with Bootstrap using R
#' @param x the sample
#' @param conf confidence level
#' @return a list including the sample estimate and BCa confidence interval
#' @examples
#' \dontrun{
#' library("bootstrap")
#' data(scor,package="bootstrap")
#' CI_BCa_Bootstrap(scor,conf=0.95)
#' }
#' @export
CI_BCa_Bootstrap <-function(x, conf) {
  x <- as.matrix(x)
  n <- nrow(x) 
  N <- 1:n
  alpha <- (1 + c(-conf, conf))/2
  zalpha <- qnorm(alpha)
  cov.hat<-cov(x)
  lamda.hat<-eigen(cov.hat)$values
  lamda.sum<-sum(lamda.hat)
  th0<-lamda.hat[1]/lamda.sum
  B<-2000
  th<-numeric(B)
  for(b in 1:B){
    i <- sample(1:n, size = n, replace = TRUE)
    scor.boot<-scor[i,]
    cov.boot<-cov(scor.boot)
    lamda.boot<-eigen(cov.boot)$values
    th[b]<-lamda.boot[1]/sum(lamda.boot)
  }
  z0 <- qnorm(sum(th < th0) / length(th))
  scor.jack<-x[-1,]
  theta.jack<-numeric(n)
  for(i in 1:n){
    scor.jack<-scor[-i,]
    cov.jack<-cov(scor.jack)
    lamda.jack<-eigen(cov.jack)$values
    theta.jack[i]<-lamda.jack[1]/sum(lamda.jack)
  }
  L <- mean(theta.jack) - theta.jack
  a <- sum(L^3)/(6 * sum(L^2)^1.5)
  adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
  limits <- quantile(th, adj.alpha, type=6)
  return(list("est"=th0, "BCa"=limits))
}