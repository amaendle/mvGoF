#aus Markov ADessohneskript.R

library(copula)
library(MASS)
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
#' \code{rosenblatt.norm} returns Rosenblatt's transformation of a normally distributed vector \code{x} with mean vector \code{Mu} and correlation matrix \code{Sigma}.
#'
#' @examples
#' x <- matrix(rnorm(50), nrow=25, ncol=2)
#' plot(x)
#' u <- rosenblatt.norm(x, rep(0,2), diag(1,2,2))
#' plot(u)
#'
#' x <- mvrnorm(n = 500, rep(0, 2), diag(1,2,2))
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
#' @seealso \code{\link[copula]{gtrafo}} for Rosenblatt's transformation on elliptical and Archimedean copulas.
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
# rosenblatt.norm <- function(x, Mu, Sigma) {
#   if (class(x)=="numeric"){
#     x <- matrix(x,1,length(x))
#   }
#   d <- dim(x)[2]
#   n <- dim(x)[1]
#   u <- matrix(NaN,dim(x)[1],dim(x)[2])
#
#   ### case k <- 1 ###
#   ### marginal distribution ###
#   k <- 1
#   mu1 <- Mu[k]
#   sigma1 <- Sigma[k,k]
#   #u[1] <- pnorm(x[1], mean = mu1, sd = sigma1, lower.tail = TRUE, log.p = FALSE)
#   tempfun <- function(x)  { pnorm(x, mean = mu1, sd = sigma1, lower.tail = TRUE, log.p = FALSE) }
#   u[,k] <- sapply( x[,k] , tempfun )
#
#   #print(x)
#   for (k in 2:d) {
#      # make x matrix
#     tMu1 <- matrix(Mu[1:(k-1)],k-1,1)
#     tMu2 <- matrix(Mu[k:d],d-k+1,1)
#     tSigma11 <- matrix(Sigma[1:(k-1),1:(k-1)],k-1,k-1)
#     tSigma21 <- matrix(Sigma[k:d,1:(k-1)],d-k+1,k-1)
#     tSigma12 <- matrix(Sigma[1:(k-1),k:d],k-1,d-k+1)
#     tSigma22 <- matrix(Sigma[k:d,k:d],d-k+1,d-k+1)
#     tSigma11I <- solve(tSigma11)
#     muk <- NULL #remove(muk)
#     sigmak <- NULL #remove(sigmak)
#
#     muk <- tMu2[1]%*%matrix(1,1,n) + (tSigma21[1,] %*% tSigma11I) %*% (t(x[,1:(k-1)]) - tMu1%*%matrix(1,1,n))
#     sigmak <- tSigma22 - (tSigma21 %*% tSigma11I %*% tSigma12)
#
#     for (l in 1:n) {
#       u[l,k] <- pnorm(x[l,k], mean = muk[1,l], sd = sigmak[1,1], lower.tail = TRUE, log.p = FALSE)
#     }
#   }
#     return(u)
# }

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
#' x <- mvrnorm(n = 200, rep(0, 2), diag(1,2,2))
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
#' @family test statistics
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
    dataSim <- mvrnorm(n=n, MuF, SigmaF)
    dataRoseSim <-rosenblatt.norm(dataSim,MuF,SigmaF)
    return(summand12nuevo(dataRoseSim, nsum))  ###########!!!!!!!
  }
  print(SumSim)
  print(quantile(SumSim, probs = seq(0.95, 1, 0.05)))
  ## fehlspezifiziert

  SumSimST = foreach(i = 1:l, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
    # print("Runde")
    # print(i)
    dataSimST <- mvrnorm(n=n, Mu, Sigma)
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
    dataSim <- mvrnorm(n=n, MuF, SigmaF)
    dataRoseSim <-rosenblatt.norm(dataSim,MuF,SigmaF)
    return(summand2tails(dataRoseSim, nsum))  ############!!!!!!!!!!!!
  }
  print(SumSim)
  print(quantile(SumSim, probs = seq(0.95, 1, 0.05)))

  SumSimST = foreach(i = 1:l, .packages = c("MASS","copula","fBasics"), .combine = c) %dopar% {
    dataSimST <- mvrnorm(n=n, Mu, Sigma)
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
#   tmpsim <- mvrnorm(130,musim ,sigmasim )
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
#     dataSim <- mvrnorm(n=n, MuF, SigmaF)
#     dataRoseSim <-rosenblatt.norm(dataSim,MuF,SigmaF)
#     return(summand12nuevo(dataRoseSim, nsum))
#   }
#   probabs <- 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)
#   print(quantile(SumSim, probs = 1-c(0.25,0.2,0.15,0.1,0.05,0.025,0.01,0.005,0.0025,0.001)  ))  #seq(0.95, 1, 0.05)))
# }
#

