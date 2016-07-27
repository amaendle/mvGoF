#' mvGoF: A package for computating several multivariate EDF goodness of fit statistics.
#'
#' Multivariate versions of the popular EDF statistics by \emph{Cramér-von Mises}, \emph{Anderson-Darling} and \emph{Kolmogorov-Smirnov} are supported.
#'
#' @section List of public functions:
#' \itemize{
#'  \item{"\code{\link{AD.stat}}"}{ Anderson-Darling}
#'  \item{"\code{\link{CvM.stat}}"}{ Cramér-von Mises}
#'  \item{"\code{\link{biKS.stat}}"}{ bivariate Kolmogorov-Smirnov}
#'  \item{"\code{\link{rosenblatt.norm}}"}{ Rosenblatt's transformation for normally distributed data}
#' }
#'
#' @docType package
#' @name mvGoF
#' @import copula
#' @examples
#' # Evaluate percentiles of the distribution by Monte-Carlo simulation
#' # sampling from independent uniforms 0-1 for F_0 specified completely by H_0
#' n <- 10 # e.g. sample size 10
#' set.seed(14072016)
#' sims <- stat.MC(n,100,biKS.stat) # 100 MC simulations (usually larger)
#' # now quantiles can be estimated
#' quantile(sims)
#' # critival value for sample size 48 and alpha=0.05
#' KS.cval(48,0.05)
#'
#' # for the Lilliefors statistic:
#' n <- 10 # e.g. sample size 10
#' set.seed(14072016)
#' sims <- stat.MC(n,100,biKS.stat) # 100 MC simulations (usually larger)
#'
#' #evaluate empirical power of statistics
#' #set parameters of true distribution
#' (mu0 <- rep(0,2))
#' (Sigma0 <- diag(0.5,2,2)+matrix(0.5,2,2))
#' #set parameters for alternative distribution
#' nMC <- 1000 # number of MC simulations
#' n <- 15 # sample size
#' eps <- 0.4 #c(0.1,0.2,0.4)
#' mu <- c(3,3) #mu <- c(3,-1)
#' tmp.KS.n15.sims <- rep(NA, nMC)
#' tmp.AD.n15.sims <- rep(NA, nMC)
#' tmp.CvM.n15.sims <- rep(NA, nMC)
#' set.seed(15072016)
#' for (i in 1:nMC) {
#'   rbts <- rosenblatt.norm( ((1-eps)*MASS::mvrnorm(n, mu0, Sigma0) + eps*MASS::mvrnorm(n, mu, Sigma0))  ,mu0,Sigma0)
#'   tmp.KS.n15.sims[i] <- biKS.stat( rbts  )
#'   tmp.AD.n15.sims[i] <- AD.stat( rbts  )
#'   tmp.CvM.n15.sims[i] <- CvM.stat( rbts  )
#' }
#' #length(t1.KS.n15.sims[t1.KS.n15.sims>KS.cval(15,0.05)])/length(t1.KS.n15.sims)
#' #length(t2.KS.n15.sims[t2.KS.n15.sims>KS.cval(15,0.05)])/length(t2.KS.n15.sims)
#' #length(t3.KS.n15.sims[t3.KS.n15.sims>KS.cval(15,0.05)])/length(t3.KS.n15.sims)
#' #t1: mu <- c(3,3), eps <- 0.1
#' #t2: mu <- c(3,3), eps <- 0.2
#' #t3: mu <- c(3,3), eps <- 0.4
#'
NULL

#' Compute Monte-Carlo simulations of test statistic \code{func}
#'
#' \code{stat.MC} returns the simulated values of function \code{func}.
#'
#' @param n Sample size.
#' @param rp Number of Monte Carlo simulations.
#' @param func Function which evaluates the test statistic.
#'
#' @examples
#' KS.n12test.sims <- stat.MC(12,100,biKS.stat)
#' # now quantiles can be estimated
#' quantile(KS.n12test.sims,probs=(1-c(0.01,0.05,0.1,0.25)))
#' @export
stat.MC <- function(n,rp,func) {
  # Evaluate percentiles of the distribution by Monte-Carlo simulation
  # sampling from independent uniforms 0-1 for F_0 specified completely by H_0
  sims <- rep(NA,rp) #matrix(NA, 10, 2, dimnames = list(c(1:10),c("X", "Y"))) #
 # set.seed(14072016)
  #system.time(
  for (i in 1:rp) {
    rbts <- matrix(runif(2*n), n, 2, dimnames = list(c(1:n),c("X", "Y")))
    sims[i] <- func( rbts  )
  }#)
  return(sims)
}
stat.MC.clust <- function(n,rp,func){
  #create rb.array
  rm(rb.array)
  #rp<-10000
  #n<-50
  rb.array <- array(NA, dim=c(rp,n,2))
  #set.seed(14072016)
  for (i in 1:rp) {
    rb.array[i,,] <- matrix(runif(2*n), n, 2, dimnames = list(c(1:n),c("X", "Y")))
  }
  #apply(rb.array[1:10,,],c(1),biKS.stat)

  stopCluster(cl)
  cl <- makeCluster(15, outfile="")
  registerDoParallel(cl)
  getDoParWorkers()
  tmp = foreach(i = 1:rp, .packages = c("np","MASS","copula","fBasics","doParallel","foreach"), .combine = rbind) %dopar% {
    #print(i)
    return(c(i,biKS.stat(rb.array[i,,])))
  }  ; return(tmp);
}

#Lilliefors
stat.MC.Lilliefors <- function(n,rp,func) {
  ##library(mclust)
  # Evaluate percentiles of the distribution by Monte-Carlo simulation
  # sampling from .. for F_0 specified completely by H_0
  sims <- rep(NA,rp)
  # set.seed(14072016)
  #system.time(
  for (i in 1:rp) {
    mvns <- MASS::mvrnorm(n,rep(0,2),diag(1,2))
    tmpfit <- mvn(modelName="Ellipsoidal",mvns)
    (tmpmu<-tmpfit$parameters$mean)
    (tmpSigma<-tmpfit$parameters$variance$Sigma)
    (rbts <- rosenblatt.norm(mvns,tmpmu,tmpSigma))
    sims[i] <- func( rbts  )
  }#)
  return(sims)
}

#' Quantiles of the Kolmogorov–Smirnov test
#'
#' \code{KS.quantiles} returns quantiles for the Kolmogorov–Smirnov test.
#'
#' @param alphas Alpha levels for which the quantiles should be returned.
#'
#' @examples
#' x <- matrix(rnorm(50), nrow=25, ncol=2)
#' @seealso \code{\link[copula]{rtrafo}} for Rosenblatt's transformation on elliptical and Archimedean copulas.
#' @export
KS.quantiles <- function(alphas=c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)){
  n<-c(10,11,12,13,14,15,20,25,30,40,50,60,80,100,150,200,300)
  return(cbind(n,rbind(
    quantile(KS.n10.sims,probs=(1-alphas)),
    quantile(KS.n11.sims,probs=(1-alphas)),
    quantile(KS.n12.sims,probs=(1-alphas)),
    quantile(KS.n13.sims,probs=(1-alphas)),
    quantile(KS.n14.sims,probs=(1-alphas)),
    quantile(KS.n15.sims,probs=(1-alphas)),
    quantile(KS.n20.sims,probs=(1-alphas)),
    quantile(KS.n25.sims,probs=(1-alphas)),
    quantile(KS.n30.sims,probs=(1-alphas)),
    quantile(KS.n40.sims,probs=(1-alphas)),
    quantile(KS.n50.sims,probs=(1-alphas)),
    quantile(KS.n60.sims,probs=(1-alphas)),
    quantile(KS.n80.sims,probs=(1-alphas)),
    quantile(KS.n100.sims,probs=(1-alphas)),
    quantile(KS.n150.sims,probs=(1-alphas)),
    quantile(KS.n200.sims,probs=(1-alphas)),
    quantile(KS.n300.sims,probs=(1-alphas)) )))
}
#' Quantiles of the Anderson-Darling test
#'
#' \code{KS.quantiles} returns quantiles for the Anderson-Darling test.
#'
#' @param alphas Alpha levels for which the quantiles should be returned.
#'
#' @examples
#' x <- matrix(rnorm(50), nrow=25, ncol=2)
#' @seealso \code{\link[copula]{rtrafo}} for Rosenblatt's transformation on elliptical and Archimedean copulas.
#' @export
AD.quantiles <- function(alphas=c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)){
  n<-c(10,11,12,13,14,15,20,25,30,40,50,60,80,100,150,200,300)
  return(cbind(n,rbind(
    quantile(AD.n10.sims,probs=(1-alphas)),
    quantile(AD.n11.sims,probs=(1-alphas)),
    quantile(AD.n12.sims,probs=(1-alphas)),
    quantile(AD.n13.sims,probs=(1-alphas)),
    quantile(AD.n14.sims,probs=(1-alphas)),
    quantile(AD.n15.sims,probs=(1-alphas)),
    quantile(AD.n20.sims,probs=(1-alphas)),
    quantile(AD.n25.sims,probs=(1-alphas)),
    quantile(AD.n30.sims,probs=(1-alphas)),
    quantile(AD.n40.sims,probs=(1-alphas)),
    quantile(AD.n50.sims,probs=(1-alphas)),
    quantile(AD.n60.sims,probs=(1-alphas)),
    quantile(AD.n80.sims,probs=(1-alphas)),
    quantile(AD.n100.sims,probs=(1-alphas)),
    quantile(AD.n150.sims,probs=(1-alphas)),
    quantile(AD.n200.sims,probs=(1-alphas)),
    quantile(AD.n300.sims,probs=(1-alphas)) )))
}
#' Quantiles of the Cramér-von Mises test
#'
#' \code{KS.quantiles} returns quantiles for the Cramér-von Mises test.
#'
#' @param alphas Alpha levels for which the quantiles should be returned.
#'
#' @examples
#' x <- matrix(rnorm(50), nrow=25, ncol=2)
#' @seealso \code{\link[copula]{rtrafo}} for Rosenblatt's transformation on elliptical and Archimedean copulas.
#' @export
CvM.quantiles <- function(alphas=c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)){
  n<-c(10,11,12,13,14,15,20,25,30,40,50,60,80,100,150,200,300)
  return(cbind(n,rbind(
    quantile(CvM.n10.sims,probs=(1-alphas)),
    quantile(CvM.n11.sims,probs=(1-alphas)),
    quantile(CvM.n12.sims,probs=(1-alphas)),
    quantile(CvM.n13.sims,probs=(1-alphas)),
    quantile(CvM.n14.sims,probs=(1-alphas)),
    quantile(CvM.n15.sims,probs=(1-alphas)),
    quantile(CvM.n20.sims,probs=(1-alphas)),
    quantile(CvM.n25.sims,probs=(1-alphas)),
    quantile(CvM.n30.sims,probs=(1-alphas)),
    quantile(CvM.n40.sims,probs=(1-alphas)),
    quantile(CvM.n50.sims,probs=(1-alphas)),
    quantile(CvM.n60.sims,probs=(1-alphas)),
    quantile(CvM.n80.sims,probs=(1-alphas)),
    quantile(CvM.n100.sims,probs=(1-alphas)),
    quantile(CvM.n150.sims,probs=(1-alphas)),
    quantile(CvM.n200.sims,probs=(1-alphas)),
    quantile(CvM.n300.sims,probs=(1-alphas)) )))
}
KSL.quantiles <- function(alphas=c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)){
  n<-c(10,11,12,13,14,15,20,25,30,40,50,60,80,100,150,200,300)
  return(cbind(n,rbind(
    quantile(KS.n10.sims,probs=(1-alphas)),
    quantile(KS.n11.sims,probs=(1-alphas)),
    quantile(KS.n12.sims,probs=(1-alphas))#,
    # quantile(KS.n13.sims,probs=(1-alphas)),
    # quantile(KS.n14.sims,probs=(1-alphas)),
    # quantile(KS.n15.sims,probs=(1-alphas)),
    # quantile(KS.n20.sims,probs=(1-alphas)),
    # quantile(KS.n25.sims,probs=(1-alphas)),
    # quantile(KS.n30.sims,probs=(1-alphas)),
    # quantile(KS.n40.sims,probs=(1-alphas)),
    # quantile(KS.n50.sims,probs=(1-alphas)),
    # quantile(KS.n60.sims,probs=(1-alphas)),
    # quantile(KS.n80.sims,probs=(1-alphas)),
    # quantile(KS.n100.sims,probs=(1-alphas)),
    # quantile(KS.n150.sims,probs=(1-alphas)),
    # quantile(KS.n200.sims,probs=(1-alphas)),
    # quantile(KS.n300.sims,probs=(1-alphas))
    )))
}

#' Compute critical value for the Kolmogorov–Smirnov test
#'
#' \code{KS.cval} returns the critical values of for the Kolmogorov–Smirnov test statistic.
#'
#' @param n Sample size.
#' @param prob Quantile.
#' @export
KS.cval <- function(n,prob) {stat.cval(n,prob,KS.quantiles)}
#' Compute critical value for the Anderson-Darling test
#'
#' \code{KS.cval} returns the critical values of for the Anderson-Darling test statistic.
#'
#' @param n Sample size.
#' @param prob Quantile.
#' @export
AD.cval <- function(n,prob) {stat.cval(n,prob,AD.quantiles)}
#' Compute critical value for the Cramér-von Mises test
#'
#' \code{KS.cval} returns the critical values of for the Cramér-von Mises test statistic.
#'
#' @param n Sample size.
#' @param prob Quantile.
#' @export
CvM.cval <- function(n,prob) {stat.cval(n,prob,CvM.quantiles)}

stat.cval <- function(n,prob,qfun){
  if (n>350) {warning("critical value due to extrapolation probably not reliable")}

  xy <- qfun(prob)
  #if (n<=80) {
    return(spline(xy[,1],xy[,2],method="natural",xout=n)$y)  #alternative without extrapolation:   return(approx(xy[,1],xy[,2],method="linear",xout=n)$y)
  # }
  # else {
  #   t15 <- as.numeric(KS.quantiles(0.05)[6,2])
  #   t50 <- as.numeric(KS.quantiles(0.05)[11,2])
  #   t60 <- as.numeric(KS.quantiles(0.05)[12,2])
  #   t80 <- as.numeric(KS.quantiles(0.05)[13,2])
  #   dslnex <- function(x) {
  #     y <- numeric(1)
  #     y[1] <- x[1]*(x[2]^60-x[2]^80)/(1-x[2]) +t80-t60
  #     #y[2] <- x[1]*(x[2]^50-x[2]^60)/(1-x[2]) +t60-t50
  #     #y[2] <- x[1]*(x[2]^15-x[2]^650)/(1-x[2]) +t50-t15
  #     return(y)
  #   }
  #   pars<-optim(c(0.5,0.99999992),dslnex,method = "Nelder-Mead")$par
  #   #require("nleqslv")
  #   #nleqslv(x=c(0.5,0.1),dslnex,method="Newton",global="none")
  #   return(t60-pars[1]*((pars[2]^60-pars[2]^n)/(1-pars[2])))#return(t80-pars[1]*((pars[2]^100-pars[2]^n)/(1-pars[2])))
  # }
}

KS.pval <- function(x,n){
  optfun <- function(p) {return(abs(KS.cval(n,p)-x))}
  p<-optim(par=0.43, optfun, method="Brent",lower=0,upper=1)$par
 # if (p>=0.99) {warning("p-value larger than 0.99")}
 # if (p<=0.01) {warning("p-value smaller than 0.01")}
  return(p)
}
AD.pval <- function(x,n){
  optfun <- function(p) {return(abs(AD.cval(n,p)-x))}
  return(optim(par=0.43, optfun, method="Brent",lower=0,upper=1)$par)
}
CvM.pval <- function(x,n){
  optfun <- function(p) {return(abs(CvM.cval(n,p)-x))}
  return(optim(par=0.43, optfun, method="Brent",lower=0,upper=1)$par)
}


#' Multivariate Cramér-von Mises statistic
#'
#' \code{CvM.stat} returns the value of the Cramér-von Mises test statistic for testing uniformity.
#'
#' The CvM test statistic is in fact the same as the squared L2 star discrepancy from \emph{Warnock (1972)}.
#' \code{u} has couloumn size equalto the dimension of the samples and line number equal to the number of observations.
#'
#' @param u Data set. Coloumn size equals the dimension of the samples and the line number equals the number of observations.
#' @examples
#' set.seed(102015)
#' u <- matrix(runif(50), ncol=2)
#' plot(u)
#' CvM.stat(u)
#' #[1] 0.1061887
#'
#' x <- MASS::mvrnorm(n = 200, rep(0, 2), diag(1,2,2))
#' plot(x)
#' u <- rosenblatt.norm(x, rep(0,2), diag(1,2,2))
#' plot(u)
#' CvM.stat(u)
#'
#' u <- 1/array(9:4, c(6,2))
#' CvM.stat(u)
#' @family "test statistics"
#' @seealso
#' \itemize{
#'  \item{Warnock, T.T., 1972. Computational investigations of low discrepancy point sets. In: Zaremba, S.K. (Ed.), Applications of Number Theory to Numerical Analysis. Academic Press, pp. 319-343.}
#'  \item{Chiu, S.N. and K.I. Liu. Generalized Cramér--von Mises goodness-of-fit tests for multivariate distributions. Computational Statistics & Data Analysis, 53(11):3817-3834, 2009.}
#' }
#' @export
CvM.stat <- function(u) {
  n <- dim(u)[1]
  s <- dim(u)[2]

  t3 <- 0
  for (i in 1:n) {
    ti<-0
    imax <- function(x){ pmax(x,u[i,])}
    ti <- t(apply(u,1,imax))
    ti <- sum(apply((1-ti),1,prod))
    t3 <- t3+ti
  }
  t3 <- t3/(n^2)

  t2 <- 0.5-((u*u)/2)
  t2 <- sum(apply(t2,1,prod))
  t2 <- 2*t2/(n)

  t1 <- (1/3)^s
  return(t1-t2+t3)
}

# superior distance D_n^+
biKS.supdist <- function(u,dat) {
  return(mecdf(u,dat) - prod(u))
}

# inferiour distance D_n^-
biKS.infdist <- function(u,dat) {
  return((-1)*biKS.supdist(u,dat))
}

# intersection points
biKS.isecpts <- function(dat) {
  n <- dim(dat)[1]
  cmbns <- t(combn(1:n, 2)) # all? combinations (only unique ones!!)
  cmbnsXcheckA <- cmbns[dat[cmbns[,2],1]>dat[cmbns[,1],1],] # combinations for which dat_2 has greater x value than dat_1
  cmbnsXcheckB <- cmbns[dat[cmbns[,2],1]<dat[cmbns[,1],1],c(2,1)]
  cmbnsXcheck <- rbind(cmbnsXcheckA,cmbnsXcheckB)
  cmbnsYcheck <- cmbnsXcheck[dat[cmbnsXcheck[,2],2]<dat[cmbnsXcheck[,1],2],]
  return(cbind(dat[cmbnsYcheck[,2],1],dat[cmbnsYcheck[,1],2]))
}



#' Kolmogorov-Smirnov test statistic
#'
#' \code{CvM.stat} returns the value of the bivariate Kolmogorov-Smirnov test statistic
#'
#' The Kolmogorov-Smirnov test statistic is in fact the same as the squared L_infinity star discrepancy.
#'
#' @param u Data set. Coloumn size equals the dimension of the samples and the line number equals the number of observations.
#' @examples
#' set.seed(102015)
#' u <- matrix(runif(50), ncol=2)
#' plot(u)
#' CvM.stat(u)
#' #[1] 0.006029274
#' biKS.stat(u)
#' #[1] 0.252292
#'
#' x <- MASS::mvrnorm(n = 200, rep(0, 2), diag(1,2,2))
#' plot(x)
#' u <- rosenblatt.norm(x, rep(0,2), diag(1,2,2))
#' plot(u)
#' CvM.stat(u)
#' #[1] 0.002172406
#' biKS.stat(u)
#' #[1] 0.1058533
#'
#' u <- 1/array(9:4, c(6,2))
#' CvM.stat(u)
#' biKS.stat(u)
#' @family "test statistics"
#' @seealso
#' \itemize{
#'  \item{Justel, A, D. Pe~na and R. Zamar, 1997. A multivariate Kolmogorov-Smirnov test of goodness of fit. Statistics & Probability Letters 35, 251-259.}
#'  \item{Chiu, S.N. and K.I. Liu. Generalized Cramér--von Mises goodness-of-fit tests for multivariate distributions. Computational Statistics & Data Analysis, 53(11):3817-3834, 2009.}
#' }
#' @export
biKS.stat <- function(u) {
  n <- dim(u)[1]
  s <- dim(u)[2]

  # according to justel pena zamar:
  # (1) compute the maximum distance in the  observed points, D_n^1 = max D_n^+(u_i)
  dn1 <- max(apply(u,1,biKS.supdist,dat=u))
  # (2) compute the maximum and minimum distances in the intersection points,
  ispts <- biKS.isecpts(u)
  dn2 <- max(apply(ispts,1,biKS.supdist,dat=u))
  dn3 <- (2/n)-min(apply(ispts,1,biKS.supdist,dat=u))
  #     D_n^2 = max_i,j {D^+_n (x_j, y_i)       | x_j > x_i, y_j < y_i} and
  #     D_n^3 = 2/n - min_i,j {D^+_n (x_j, y_i) | x_j > x_i, y_j < y_i }
  # (3) compute the maximum distance among the projections of the observed points on the right unit-square
  #     border, D_n^4 = 1/n- min D^+_n (1,y_i)
  dn4 <- (1/n) - min(apply(cbind(rep(1,n),u[,2]),1,biKS.supdist,dat=u))
  # (4) compute the maximum distance among the projections of the observed points on the top unit-square border,
  #     D_n^5 = 1/n - min D^+_n (X_i, 1)
  dn5 <- (1/n) - min(apply(cbind(u[,1],rep(1,n)),1,biKS.supdist,dat=u))
  # (5) compute the maximum D_n = max{D_n^1, D_n^2, D_^3n, D_n^4, D_n^5}.
  return(max(dn1,dn2,dn3,dn4,dn5))
}

#' approximative Kolmogorov-Smirnov test statistic
#'
#' \code{KSa.stat} returns the value of the approximative Kolmogorov-Smirnov test statistic
#'
#' The statistic used by this function \eqn{D_n^*} from
#' Justel, A.; Peña, D.; Zamar, R. (1997). "A multivariate Kolmogorov–Smirnov test of goodness of fit". Statistics & Probability Letters 35 (3): 251–259. doi:10.1016/S0167-7152(97)00020-5.
#' as an approximation of the Kolmogorov-Smirnov test statistic.
#'
#' @param u Data set. Coloumn size equals the dimension of the samples and the line number equals the number of observations.
#' @examples
#' set.seed(102015)
#' u <- matrix(runif(50), ncol=2)
#' plot(u)
#' CvM.stat(u)
#' #[1] 0.006029274
#' KSa.stat(u)
#' #[1] 0.252292
#'
#' x <- MASS::mvrnorm(n = 200, rep(0, 2), diag(1,2,2))
#' plot(x)
#' u <- rosenblatt.norm(x, rep(0,2), diag(1,2,2))
#' plot(u)
#' CvM.stat(u)
#' #[1] 0.002172406
#' biKS.stat(u)
#' #[1] 0.1058533
#'
#' u <- 1/array(9:4, c(6,2))
#' CvM.stat(u)
#' KSa.stat(u)
#' @family "test statistics"
#' @seealso
#' \itemize{
#'  \item{Justel, A, D. Pe~na and R. Zamar, 1997. A multivariate Kolmogorov-Smirnov test of goodness of fit. Statistics & Probability Letters 35, 251-259.}
#'  \item{Chiu, S.N. and K.I. Liu. Generalized Cramér--von Mises goodness-of-fit tests for multivariate distributions. Computational Statistics & Data Analysis, 53(11):3817-3834, 2009.}
#' }
#' @export
KSa.stat <- function(u) {
  n <- dim(u)[1]
  s <- dim(u)[2]

  # according to justel pena zamar:
  # (1) compute the maximum distance in the  observed points, D_n^1 = max D_n^+(u_i)
  dn1 <- max(apply(u,1,biKS.supdist,dat=u))
  # # (2) compute the maximum and minimum distances in the intersection points,
  # ispts <- biKS.isecpts(u)
  # dn2 <- max(apply(ispts,1,biKS.supdist,dat=u))
  # dn3 <- (2/n)-min(apply(ispts,1,biKS.supdist,dat=u))
  # #     D_n^2 = max_i,j {D^+_n (x_j, y_i)       | x_j > x_i, y_j < y_i} and
  # #     D_n^3 = 2/n - min_i,j {D^+_n (x_j, y_i) | x_j > x_i, y_j < y_i }
  # # (3) compute the maximum distance among the projections of the observed points on the right unit-square
  # #     border, D_n^4 = 1/n- min D^+_n (1,y_i)
  # dn4 <- (1/n) - min(apply(cbind(rep(1,n),u[,2]),1,biKS.supdist,dat=u))
  # # (4) compute the maximum distance among the projections of the observed points on the top unit-square border,
  # #     D_n^5 = 1/n - min D^+_n (X_i, 1)
  # dn5 <- (1/n) - min(apply(cbind(u[,1],rep(1,n)),1,biKS.supdist,dat=u))
  # # (5) compute the maximum D_n = max{D_n^1, D_n^2, D_^3n, D_n^4, D_n^5}.
  return(dn1) #return(max(dn1,dn2,dn3,dn4,dn5))
}

# t2 <- 0.5-(u*u)/2
# t2 <- apply(t2,1,prod)
# t2 <- sum(t2)


#roxygen2::roxygenise()

# dummy Kolmogorov Smirnov
mult5KS <- function(data) {
  n<-dim(data)[1]
  maxret <-0
  for (i1 in 1:n){print(i1)
    for (i2 in 1:n){print(i2)
      for (i3 in 1:n){#print(i3)
        for (i4 in 1:n){#print(i4)
          for (i5 in 1:n){
            maxret<-max(maxret,multKSdist(c(data[i1,1],data[i2,2],data[i3,3],data[i4,4],data[i5,5]),data))
          }
          #print(maxret)
        }
      }
      print(maxret)
    }
  }
  return(maxret)
}
multKSdist <- function(u,data) {
  myprod<-prod(u)
  return(max(mecdf0(u,data)-myprod,mecdf1(u,data)-myprod))
}
mecdf <- function(u,data){
  d<- length(u)
  if (d!=dim(data)[2]) {
    print("mecdf0 error, different dimension of u and data")
    return(NULL)
  }
  n <- dim(data)[1]

  data <- (-1)*data
  data <- sweep(data,2,matrix(u,nrow=1),"+")
  return(sum( apply(apply(data, c(1,2), Heaviside1),1,prod) )/n)
}
Heaviside0 <- function(x) {
  if (x > 0) return(1)
  else return(0)
}
Heaviside1 <- function(x) {
  if (x < 0) return(0)
  else return(1)
}


### now content of mADpack.R follows

#aus Markov ADessohneskript.R

#library(copula)
#library(MASS)
# library("parallel")
# library("foreach")
# library("doParallel")

########################################################
# rosenblatt.norm <- function(x, Mu, Sigma) 		##
# fillCol <- function(mat, col, x)				##
# getV <- function(u,m)						##
# StatC1 <- function(x, nsum)					## die Komponente hinter 1/n
# StatC2 <- function(x, nsum)					## die Komponente hinter 1/n-2
# rekStatC1 <- function(xalt, xneu, nsum)			## gibt an wie viel die Komponente hinter 1/n durch die neuen Daten waechst
# rekStatC2 <- function(xalt, xneu, nsum)			## gibt an wie viel die Komponente hinter 1/n durch die neuen Daten waechst
# buildADStat <- function(c1, c2, n, dim,nsum)		## berechnet die Statistik auf Basis der beiden genannten Komponenten
# AD.stat <- function(x, nsum)					## wie oben, berechnet die Statistik aber direkt mit Hilfe von StatC1, StatC2
# summand12nuevo <- function(x, nsum)			## computes the whole statistic directly
# summand2tails <- function(x, nsum) 			## sum of the statistic for both tails
# summand2tmax <- function(x, nsum) 			## max of the statistic for both tails
### Eine Art Backtesting:					##
# simNTmc <- function(n,MuF,SigmaF,k,Mu,Sigma,l,nsum)	## Simulation NV(MuF,SigmaF) gegen NV(Mu,Sigma), auf summand12nuevo(dataRoseSim, nsum)
# simNTmc2 <- function(n,MuF,SigmaF,k,Mu,Sigma,l,nsum)## basierend auf summand2tails(dataRoseSim, nsum)
### Kritische Werte bestimmen, auf Basis von runif	##
# cvalMC <- function(n,dim,k,nsum)				## critical value distribution free
# cvalMCZ <- function(n,dim,k,rn,nsum)			## random numbers given
# cvalMCndf <- function(n,MuF,SigmaF,k,nsum)		## critical value based on true (normal) distribution
#									##
########################################################

#' Rosenblatt's transformation for normally distributed data
#'
#' \code{rosenblatt.norm} returns Rosenblatt's transformation of normally distributed vectors which are the lines of the array \code{x}, with mean vector \code{Mu} and correlation matrix \code{Sigma}.
#'
#' @param x Array containing the data which have to be transformed. Each line represents a data vector with the coloumn number as dimension.
#' @param Mu Mean vector.
#' @param Sigma Correlation matrix.
#'
#' @examples
#' x <- matrix(rnorm(50), nrow=25, ncol=2)
#' plot(x)
#' u <- rosenblatt.norm(x, rep(0,2), diag(1,2,2))
#' plot(u)
#'
#' x <- MASS::mvrnorm(n = 500, rep(0, 2), diag(1,2,2))
#' plot(x)
#' u <- rosenblatt.norm(x, rep(0,2), diag(1,2,2))
#' plot(u)
#'
#' set.seed(5)
#' library("gumbel")
#' gdata <- qlnorm(rgumbel(12,alpha=2,dim=5))
#' pairs(gdata)
#' sigma2 <- matrix((1/sqrt(2)),nrow=5,ncol=5)-diag(((1/sqrt(2))-1),5)
#' sigma2
#' ndata <- qnorm(plnorm(gdata))
#' pairs(ndata)
#' rbdata <- rosenblatt.norm(ndata, rep(0,5),sigma2)
#' pairs(rbdata)
#' @seealso \code{\link[copula]{rtrafo}} for Rosenblatt's transformation on elliptical and Archimedean copulas.
#' @export
rosenblatt.norm <- function(x, Mu, Sigma) {
  if (class(x)=="numeric"){
    x <- matrix(x,1,length(x))
  }
  d <- dim(x)[2]
  n <- dim(x)[1]
  u <- matrix(NaN,dim(x)[1],dim(x)[2])

  ### case k <- 1 ###
  ### marginal distribution ###
  k <- 1
  mu1 <- Mu[k]
  sigma1 <- Sigma[k,k]
  tempfun <- function(x)  { pnorm(x, mean = mu1, sd = sqrt(sigma1), lower.tail = TRUE, log.p = FALSE) }    #!!!!!!!!
  u[,k] <- sapply( x[,k] , tempfun )

  #print(x)
  for (k in 2:d) {
    # make x matrix
    tMu1 <- matrix(Mu[1:(k-1)],k-1,1)
    tMu2 <- matrix(Mu[k:d],d-k+1,1)
    tSigma11 <- matrix(Sigma[1:(k-1),1:(k-1)],k-1,k-1)
    tSigma21 <- matrix(Sigma[k:d,1:(k-1)],d-k+1,k-1)
    tSigma12 <- matrix(Sigma[1:(k-1),k:d],k-1,d-k+1)
    tSigma22 <- matrix(Sigma[k:d,k:d],d-k+1,d-k+1)
    tSigma11I <- solve(tSigma11)  # inverse matrix
    muk <- NULL
    sigmak <- NULL

    muk <- tMu2[1]%*%matrix(1,1,n) + (tSigma21[1,] %*% tSigma11I) %*% (t(x[,1:(k-1)]) - tMu1%*%matrix(1,1,n))
    #  sigmak <- tSigma22 - (tSigma21 %*% tSigma11I %*% tSigma12)
    sigmak <- tSigma22[1,1] - (tSigma21[1,] %*% tSigma11I %*% tSigma12[,1])
    for (l in 1:n) {
      u[l,k] <- pnorm(x[l,k], mean = muk[1,l], sd = sqrt(sigmak[1,1]), lower.tail = TRUE, log.p = FALSE)    #!!!!!!!!
    }
  }
  return(u)
}

fillCol <- function(mat, col, x) {
  mat[,col]<-x
  return(mat)
}

getV <- function(u,m) {
  d <-  dim(m)[2]
  vplus <- matrix(1,0,d)
  vminus <- matrix(1,0,d)
  # k <- 1
  vplus <- rbind(vplus,u)
  vminus <- rbind(vminus,fillCol(vplus,1,m[1]))
  for (k in 2:d) {
    vplusback <- vplus
    vplus <- rbind(vplus,fillCol(vminus,k,m[k]))
    vminus <- rbind(vminus,fillCol(vplusback,k,m[k]))
  }
  erg=list(vplus,vminus)
  names(erg) = c("vplus","vminus")
  return(erg)
}

StatC1 <- function(x, nsum) {
  d <- dim(x)[2]
  n <- dim(x)[1]
  u <- rep(1,d)
  v <- c()
  s<-0
  prodsum <- 0
  plogx <- function(x) {
    return(polylog(z = x, s = d, n.sum = nsum)) #d dimensionen
  }
  for (k in 1:n) {
    # Produkt
    proddlvn <- prod(-log(x[k,]))
    prodsum  <- prodsum +proddlvn
    for (j in (k+1):n) {
      if (k+1>n)
        break
      v <- apply(rbind(x[k,],x[j,]),2,max)
      erg<-getV(matrix(u,1,d),matrix(v,1,d))
      # Produkt
      proddlvn <- prod(-log(v))
      s <- s+sum(plogx(apply(erg$vplus,1,prod)))-sum(plogx(apply(erg$vminus,1,prod)))+proddlvn
    }
  }
  s <- 2*s + prodsum
  return(s)
}
StatC2 <- function(x, nsum) {
  d <- dim(x)[2]
  n <- dim(x)[1]
  u <- rep(1,d)
  v <- c()
  s2<-0
  plogx <- function(x) {
    return(polylog(z = x, s = d, n.sum = nsum)) #d dimensionen
  }
  for (k in 1:n) {
    erg<-getV(matrix(u,1,d),matrix(x[k,],1,d))
    s2 <- s2 + sum(plogx(apply(erg$vplus,1,prod)))-sum(plogx(apply(erg$vminus,1,prod)))
  }
  return(s2)
}
StatC3 <- function(n,d,nsum) {
  n*(polylog(z = 1, s = d, n.sum = nsum)-1)
}


rekStatC1 <- function(xalt, xneu, nsum) {
  d <- dim(xalt)[2]
  if (dim(xneu)[2]!=d) print("oh, das geht aber nicht!")
  nalt <- dim(xalt)[1]
  nneu <- dim(xneu)[1]
  n <- nalt+nneu
  x <- rbind(xalt,xneu)
  u <- rep(1,d)
  v <- c()
  s<-0
  prodsum <- 0
  plogx <- function(x) {
    return(polylog(z = x, s = d, n.sum = nsum)) #d dimensionen
  }
  for (k in (nalt+1):n) {
    if ((nalt+1)>n)
      break
    # Produkt
    prodsum  <- prodsum + prod(-log(x[k,]))
    for (j in 1:(k-1)) {
      if (k-1<1)
        break
      v <- apply(rbind(x[k,],x[j,]),2,max)
      erg<-getV(matrix(u,1,d),matrix(v,1,d))
      # Produkt
      proddlvn <- prod(-log(v))
      s <- s+sum(plogx(apply(erg$vplus,1,prod)))-sum(plogx(apply(erg$vminus,1,prod)))+proddlvn
    }
  }
  s <- 2*s + prodsum
  return(s)
}
rekStatC2 <- function(xalt, xneu, nsum) {
  d <- dim(xalt)[2]
  if (dim(xneu)[2]!=d) print("oh, das geht aber nicht!")
  nalt <- dim(xalt)[1]
  nneu <- dim(xneu)[1]
  n <- nalt+nneu
  x <- rbind(xalt,xneu)
  u <- rep(1,d)
  v <- c()
  s2<-0
  plogx <- function(x) {
    return(polylog(z = x, s = d, n.sum = nsum)) #d dimensionen
  }
  for (k in (nalt+1):n) {
    if ((nalt+1)>n) break
    erg<-getV(matrix(u,1,d),matrix(x[k,],1,d))
    s2 <- s2 + sum(plogx(apply(erg$vplus,1,prod)))-sum(plogx(apply(erg$vminus,1,prod)))
  }
  return(s2)
}

buildADStat <- function(c1, c2, n, dim,nsum) {
  return( (c1/n) + ((1/n)-2)*c2 + n*(polylog(z = 1, s = dim, n.sum = nsum)-1) )
}

#' Multivariate Anderson-Darling test statistic
#'
#' \code{AD.stat} returns the multivariate Anderson-Darling test statistic for uniformly distributed data.
#'
#' \code{AD.stat} uses \code{polylog} from the \code{copula} package.
#'
#' @param x The sample.
#' @param nsum See \code{copula::polylog}.
#'
#' @examples
#' set.seed(102015)
#' u <- matrix(runif(50), ncol=2)
#' plot(u)
#' AD.stat(u)
#' #[1] 1.570951
#'
#' x <- matrix( seq(0.01, 0.99, length=100), ncol=2)
#' plot(x)
#' AD.stat(x)
#'
#' x <- matrix(c(0.5,0.1), ncol=2)
#' plot(x)
#' AD.stat(x)
#'
#' x <- MASS::mvrnorm(n = 200, rep(0, 2), diag(1,2,2))
#' plot(x)
#' u <- rosenblatt.norm(x, rep(0,2), diag(1,2,2))
#' plot(u)
#' AD.stat(u)
#'
#' # true situation: Data from Gumbel copula with known corellation.
#' # false assumption: Gauss copula
#' # Rosenblatt#s transformation under wrong assumptions and computation of AD.stat
#' set.seed(5)
#' library("gumbel")
#' gdata <- qlnorm(rgumbel(50,alpha=2,dim=5))
#' pairs(gdata)
#' sigma2 <- matrix((1/sqrt(2)),nrow=5,ncol=5)-diag(((1/sqrt(2))-1),5)
#' sigma2
#' ndata <- qnorm(plnorm(gdata))
#' pairs(ndata)
#' rbdata <- rosenblatt.norm(ndata, rep(0,5),sigma2)
#' pairs(rbdata)
#' AD.stat(rbdata)
#' @family "test statistics"
#' @export
AD.stat <- function(x, nsum=1e3) {
  d <- dim(x)[2]
  n <- dim(x)[1]
  # print(StatC1(x,nsum)/n) # wird manchmal unendlich
  return(StatC1(x,nsum)/n + ((1/n)-2) * StatC2(x,nsum)  + n*(copula::polylog(z = 1, s = d, n.sum = nsum)-1))
}
summand12nuevo <- function(x, nsum) {
  d <- dim(x)[2]
  n <- dim(x)[1]
  u <- rep(1,d)
  v <- c()
  s<-0
  plogx <- function(x) {
    return(copula::polylog(z = x, s = d, n.sum = nsum)) #d dimensionen
  }
  for (k in 1:n) {
    for (j in (k+1):n) {
      if (k+1>n)
        break
      v <- apply(rbind(x[k,],x[j,]),2,max)
      erg<-getV(matrix(u,1,d),matrix(v,1,d))
      # Produkt
      proddlvn <- prod(-log(v))
      s <- s+sum(plogx(apply(erg$vplus,1,prod)))-sum(plogx(apply(erg$vminus,1,prod)))+proddlvn
    }
  }
  s <- 2*s
  s2<-0
  for (k in 1:n) {
    # Produkt
    proddlvn <- prod(-log(x[k,]))
    s <- s+proddlvn #sum(plogx(apply(erg$vplus,1,prod)))-sum(plogx(apply(erg$vminus,1,prod)))+proddlvn
    erg<-getV(matrix(u,1,d),matrix(x[k,],1,d))
    s2 <- s2 + sum(plogx(apply(erg$vplus,1,prod)))-sum(plogx(apply(erg$vminus,1,prod)))
  }
  #  print(n*(copula::polylog(z = 1, s = d, n.sum = nsum)-1))
  #  print(s/n) # wird manchmal unendlich
  #  print(((1/n)-2) * sum(s2))
  return(s/n + ((1/n)-2) * sum(s2)  + n*(copula::polylog(z = 1, s = d, n.sum = nsum)-1))
}

summand2tails <- function(x, nsum) {
  return(summand12nuevo(x, nsum)+summand12nuevo(1-x, nsum))
}

summand2tmax <- function(x, nsum) {
  return(max(summand12nuevo(x, nsum),summand12nuevo(1-x, nsum)))
}

simNTmc <- function(n,MuF,SigmaF,k,Mu,Sigma,l,nsum) {
  SumSim = foreach(i = 1:k, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
    # print("Runde")
    # print(i)
    dataSim <- MASS::mvrnorm(n=n, MuF, SigmaF)
    dataRoseSim <-rosenblatt.norm(dataSim,MuF,SigmaF)
    return(summand12nuevo(dataRoseSim, nsum))  ###########!!!!!!!
  }
  print(SumSim)
  print(quantile(SumSim, probs = seq(0.95, 1, 0.05)))
  ## fehlspezifiziert

  SumSimST = foreach(i = 1:l, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
    # print("Runde")
    # print(i)
    dataSimST <- MASS::mvrnorm(n=n, Mu, Sigma)
    dataRoseSimST <-rosenblatt.norm(dataSimST,MuF,SigmaF)
    return(summand12nuevo(dataRoseSimST, nsum))
  }
  Fn <- ecdf(SumSim)
  Fn(SumSimST)
  print(SumSimST)
  print("Anzahl Ablehnungen:")
  print(length(Fn(SumSimST)[Fn(SumSimST)>0.95]))
}

simNTmc2 <- function(n,MuF,SigmaF,k,Mu,Sigma,l,nsum) {
  SumSim = foreach(i = 1:k, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
    dataSim <- MASS::mvrnorm(n=n, MuF, SigmaF)
    dataRoseSim <-rosenblatt.norm(dataSim,MuF,SigmaF)
    return(summand2tails(dataRoseSim, nsum))  ############!!!!!!!!!!!!
  }
  print(SumSim)
  print(quantile(SumSim, probs = seq(0.95, 1, 0.05)))

  SumSimST = foreach(i = 1:l, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
    dataSimST <- MASS::mvrnorm(n=n, Mu, Sigma)
    dataRoseSimST <-rosenblatt.norm(dataSimST,MuF,SigmaF)
    return(summand2tails(dataRoseSimST, nsum)) ###############!!!!!!!!
  }
  Fn <- ecdf(SumSim)
  Fn(SumSimST)
  print(SumSimST)
  print("Anzahl Ablehnungen:")
  print(length(Fn(SumSimST)[Fn(SumSimST)>0.95]))
}

# critical value distribution free
cvalMC <- function(n,dim,k,nsum) {
  SumSim = foreach(i = 1:k, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
    return(summand12nuevo( matrix(runif(n*dim),n,dim), nsum)) # runif(dim)
  }
  probabs <- 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)
  #print(quantile(SumSim, probs = 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001), na.rm = FALSE, names = TRUE, type = 1,  ))
  nt(quantile(SumSim, probs = 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001), na.rm = FALSE, names = TRUE, type = 2,  ))
  #print(quantile(SumSim, probs = 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001), na.rm = FALSE, names = TRUE, type = 3,  ))
  print(quantile(SumSim, probs = 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)  ))  #seq(0.95, 1, 0.05)))
}

# ---
#   k<-1e6
# nsum<-1e3
# n<-10
# d<-5
# cvaluesADn10d5 = foreach(i = 1:k, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
#   print(i)
#   return(summand12nuevo( matrix(runif(n*d),n,d), nsum)) # runif(dim)
# }
# save(cvaluesADn10d5, file = "cvaluesADn10d5.RData")
# save.image()
# unlink("cvaluesADn10d5.RData")
# print(quantile(cvaluesADn10d5, probs = c(0.85,0.9,0.95,0.975,0.99)))
# ---
#   cvaluesADn103d7c3 <- NULL
# n<-130
# d <- 7
# sims<-10e4
# musim <- rep(0,d)
# sigmasim <- diag(1,d)
# cvaluesADn103d7c3 = foreach(i = 1:sims, .packages = c("MASS","copula","fBasics","mclust"), .combine = c) %dopar% {
#   tmpsim <- MASS::mvrnorm(130,musim ,sigmasim )
#   tmpfit <- mvn(modelName = "Ellipsoidal", tmpsim)
#   tmpmu <- tmpfit$parameters$mean
#   tmpsigma <- tmpfit$parameters$variance$Sigma
#   print(i)
#   return(summand12nuevo(rosenblatt.norm(tmpsim,tmpmu ,tmpsigma ),nsum))
# }
# save(cvaluesADn103d7c3, file = "cvaluesADn103d7c3.RData")
# save.image()
# unlink("cvaluesADn103d7c3.RData")
# print(quantile(cvaluesADn103d7c3, probs = c(0.85,0.9,0.95,0.975,0.99)))
# ----
#
#
#   # wie oben mit uebergebenen Zufallszahlen
#   # k Anzahl Stichproben, z.B: 10000
#   # dim Dimensionen, z.B. 2
#   # n Stichprobengroesse, z.B. 10:15,20,25,...,300
#   cvalMCZ <- function(n,dim,k,rn,nsum) {
#     #- n dim k muessen auch > 1
#     if ((dim)>dim(rn)[2]) print("oh, das geht aber nicht! Wrong number of random number dimensions.")
#     if ((n)>dim(rn)[1]) print("oh, das geht aber nicht! Wrong number of sample size.")
#     if ((k)>dim(rn)[3]) print("oh, das geht aber nicht! Wrong number of number of samples.")
#     if (length(dim(rn))!=3) print("oh, das geht aber nicht! falsche Dimension von zf.")
#     SumSim = foreach(i = 1:k, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
#       return(summand12nuevo( matrix(rn[(1:n),(1:dim),i],n,dim), nsum)) # runif(dim)
#     }
#     probabs <- 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)
#     print(quantile(SumSim, probs = 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)  ))  #seq(0.95, 1, 0.05)))
#   }
#
#
# # not dist free??
# cvalMCndf <- function(n,MuF,SigmaF,k,nsum) {
#   SumSim = foreach(i = 1:k, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
#     dataSim <- MASS::mvrnorm(n=n, MuF, SigmaF)
#     dataRoseSim <-rosenblatt.norm(dataSim,MuF,SigmaF)
#     return(summand12nuevo(dataRoseSim, nsum))
#   }
#   probabs <- 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)
#   print(quantile(SumSim, probs = 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)  ))  #seq(0.95, 1, 0.05)))
# }
#


