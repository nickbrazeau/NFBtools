#' @title playingwithcpp
#'
#' @export


library(Rcpp)


rnbincpp <- cppFunction("
  int rnbinom1(double lambda, double gamma) {
    double ret = R::rnbinom(lambda/(gamma-1), 1/gamma);
    return(ret);
  }")


rpoiscpp <- cppFunction("
  int rpois1(double lambda) {
    double ret = R::rpois(lambda);
    return(ret);
  }")



rnbincppwrapper <- function(n, lambda, gamma){
  lambda <- rep(lambda, n)
  gamma <- rep(gamma, n)
  ret <- rep(NA, n)
  for(i in 1:n){
   ret[i] <- rnbincpp(lambda[i], gamma[i])
  }
  return(ret)
}


rpoiscppwrapper <- function(n, lambda){
  lambda <- rep(lambda, n)
  ret <- rep(NA, n)
  for(i in 1:n){
    ret[i] <- rpoiscpp(lambda[i])
  }
  return(ret)
}



set.seed(44)

lambda <- -0.5
gamma <- 5

cppNBret <- rnbincppwrapper(1e4, lambda, gamma)
hist(cppNBret)
mean(cppNBret)
sd(cppNBret)
summary(cppNBret)

rNBret <- rnbinom(n = 1e4, mu = lambda, size = gamma)
hist(rNBret)
mean(rNBret)
sd(rNBret)
summary(rNBret)

cpppoisret <- rpoiscppwrapper(1e4, lambda)
hist(cpppoisret)
mean(cpppoisret)
sd(cpppoisret)

rpoisret <- rpois(1e4, lambda)
hist(rpoisret)
mean(rpoisret)
sd(rpoisret)





