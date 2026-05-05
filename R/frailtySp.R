## usethis namespace: start
#' @importFrom stats runif aggregate stepfun quantile model.matrix optim
#' @importFrom graphics abline title
#' @importFrom utils flush.console
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom Rcpp sourceCpp
#' @useDynLib SpaCP, .registration = TRUE
## usethis namespace: end
NULL
## usethis namespace: start
#' @import RcppEigen
## usethis namespace: end
NULL
#' @import survival
#' @import Matrix
#' @import nleqslv
#' @import doParallel
#' @import foreach
#' @import kernlab
GenRanEf <- function(M = NULL, eta, sigma2, regular = FALSE, Plot = FALSE, m1 = NULL, m2 = NULL){
  if (is.null(m1) || is.null(m2)) {
    if (is.null(M)) stop("M or m1, m2")
    factors <- which(M %% 1:floor(sqrt(M)) == 0)
    m1 <- max(factors)
    m2 <- M / m1
    if (m1 * m2 != M) stop("M?")
  } else {
    if (m1 * m2 != M) stop("m1 * m2 = M")
  }
  grid_size1 <- 1 / m1
  grid_size2 <- 1 / m2
  if (regular == FALSE) {
    cells <- expand.grid(ix = 0:(m1-1), iy = 0:(m2-1))
    k <- nrow(cells)
    x = (cells$ix + runif(k)) / m1
    y = (cells$iy + runif(k)) / m2
    centers <- cbind(x, y)
  } else {
    x <- seq(grid_size1 / 2, 1 - grid_size1 / 2, by = grid_size1)
    y <- seq(grid_size2 / 2, 1 - grid_size2 / 2, by = grid_size2)
    centers <- expand.grid(x = x, y = y)
  }
  if (Plot == TRUE) {
    plot(centers, xlim = c(0, 1), ylim = c(0, 1), pch = 19, xlab = "x", ylab = "y")
    abline(v = seq(0, 1, length.out = m1 + 1), h = seq(0, 1, length.out = m2 + 1), lty = 2, col = "grey70")
    title(main = paste("Number of Random Effects: ", m1, "x", m2, "\n Regular: ", regular))
  }
  rbf <- kernlab::rbfdot(sigma = (1/eta)^2 )
  rho <- kernlab::kernelMatrix(rbf, as.matrix(centers) )
  Sigma <- sigma2 * rho
  d2 <- -log(rho)*eta^2
  rpi <- MASS::mvrnorm(n = 1, mu = rep(0, m1 * m2), Sigma = Sigma)
  return(list(centers = centers, rpi = rpi, d2 = d2, rho = rho))
}
#' @export
GenData<-function(n, M, beta, lambda0, sigma2, eta, Random = F, cc){
  if (length(n)==1){
    N <- n*M
    ni <- rep(n, M)
  }else{
    N <- sum(n)
    ni <- n
  }
  q <- dim(beta)[1]
  Xij <- MASS::mvrnorm(n = N, mu=rep(0, q), Sigma = diag(q) )
  Mi <- rep(1:M, times=ni)
  u <- runif(N)
  ri <- GenRanEf(M=M, eta = eta, sigma2 = sigma2, regular = T)
  Loc <- ri$centers
  d2 <- ri$d2
  allri <- rep(ri$rpi, times = ni)
  if (Random) {
    Tm <- -log(1-u)/(lambda0*exp(Xij%*%beta+allri))
  }else{
    Tm <- -log(1-u)/(lambda0*exp(Xij%*%beta))
  }
  cn <- runif(N, 0, cc)
  delta <- ifelse(Tm <= cn, 1, 0)
  Yij <- pmin(Tm, cn)
  return( list(Tij = Tm, Yij = Yij, Xij = Xij, delta=delta, ri=ri, Loc=Loc, Mi=Mi, d2=d2) )
}
#' @export
MyLambda0 <- function(dt, beta, rhat, Random = T){
  myyij <- c(dt$Yij)
  Yij <- outer(myyij, myyij, FUN = ">=")
  myxij <- dt$Xij
  #allr <- rhat[dt$Mi]
  if (Random){
    allr <- rhat[dt$Mi]
    b <- exp(myxij %*% beta + allr)
    sr1 <- apply(sweep(Yij, 1, b, FUN = "*"), 2,
                 FUN = function(x) aggregate(x ~ dt$Mi, FUN = "sum")$x)
    sr0 <- colSums(sr1)
    dt0 <- data.frame(Yij=dt$Yij, Sr0=sr0, delta=dt$delta)
    dt1 <- dt0[dt0$delta==1,]
    dtNew <- dt1[order(dt1$Yij), ]
    Cumsr0 <- cumsum(1/dtNew$Sr0)
    mystepfun <- stepfun(x = dtNew$Yij, y=c(0, Cumsr0) )
  }else{
    cox1 <- survival::coxph(survival::Surv(dt$Yij, dt$delta) ~ dt$Xij)
    bh <- survival::basehaz(cox1)
    mystepfun <- stepfun(bh$time, c(0, bh$hazard))
  }
  return(Lambda0 = mystepfun)
}
#' @export
GetWeight <- function(dt){
  dat <- data.frame(y=c(dt$Yij), delta = c(dt$delta))
  cox.censor <- survival::survfit(Surv(y, (1-delta))~1, data=dat)
  id0 <- which(dat$delta==1)
  n1 <- length(id0)
  dat.event.Yij <- dt$Yij[id0, ]
  dat.event.Xij <- dt$Xij[id0, ]
  dat.event.Mij <- dt$Mi[id0]
  tto <- dat.event.Yij
  p.tto <- summary(cox.censor, tto)$surv
  mass.tto <- 1/p.tto
  mass.tto <- mass.tto / sum(mass.tto)
  return(mass.tto)
}
#' Parameters Estimation
#'
#' This function estimates the model parameters for spatial survival data
#'
#' @param dt A list object containing the training data, with components
#'   \code{Yij} (survival time), \code{Xij} (covariate matrix),
#'   \code{delta} (event indicator), \code{Mi} (spatial location),
#'   and \code{d2} (distance matrix between locations).
#' @param diff.tol A numeric value specifying the convergence tolerance.
#'   The algorithm stops when the maximum parameter change is less than
#'   \code{diff.tol}. Default is \code{1e-2}.
#' @param max_iter A positive integer specifying the maximum number of
#'   iterations. Default is \code{100}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{parms}{A numeric vector of the estimated parameters,
#'       including regression coefficients, variance component
#'       \code{sigma2}, and spatial parameter \code{eta}.}
#'     \item{rhat}{A numeric vector of the estimated random effects.}
#'     \item{iter}{The number of iterations until convergence.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Load built-in dataset
#' data(train.dt)
#'
#' # Estimate model parameters
#' parms <- Findparms(dt = train.dt)
#'
#' # View results
#' parms
#' }
#'
#' @export
Findparms <- function(dt, diff.tol = 1e-2, max_iter = 100){
  ###
  #library(Matrix)
  q <- dim(dt$Xij)[2]
  Yij.t <- dt$Yij[order(dt$Yij, decreasing = F)]
  nrisk <- sapply(Yij.t, FUN = function(x)  Yij.t >= x )
  #nriskSp<- as(nrisk, "dgCMatrix")
  nriskSp <- Matrix::Matrix(nrisk, sparse = T)
  nriskSp <- methods::as(nrisk, Class = "dgCMatrix")
  X <- model.matrix(~-1+dt$Xij+as.factor(dt$Mi))
  X.t <- X[order(dt$Yij, decreasing = F), ]
  Mi.t <- dt$Mi[order(dt$Yij, decreasing = F) ]
  delta.t <- dt$delta[order(dt$Yij, decreasing = F) ]
  ###
  initBeta <- matrix(coef(survival::coxph(survival::Surv(dt$Yij, dt$delta) ~ dt$Xij)), ncol = 1)
  initsigma2 <- 0.5
  initeta <- 0.5
  #initr <- rep(0, length(dt$ri$rpi))
  initr <- rep(0, length(unique(dt$Mi)))
  rhat <- nleqslv::nleqslv(x = initr, fn = FindrCpp, jac = MyjacCpp,
                           d2 = dt$d2, Mi=Mi.t, delta=delta.t, nrisk = nriskSp, X = X.t,
                           beta=initBeta, eta=initeta, sigma2=initsigma2, control = list(ftol=1e-2))
  initTheta <- c(log(initsigma2), initeta)
  initparms <- c(initBeta, initTheta)
  res <- optim(initparms, l21Cpp,
               rhat = rhat$x, Mi=Mi.t, delta=delta.t, nrisk = nriskSp, X = X.t, d2=dt$d2,
               method = "L-BFGS-B", control = list(pgtol = 1e-2))
  iter <- 0
  diff <- 1
  while (diff > diff.tol && iter < max_iter) {
    iter <- iter + 1
    #initBeta <- matrix(res$par[1:2], ncol = 1)
    initBeta <- matrix(res$par[1:q], ncol = 1)
    #initsigma2 <- exp(res$par[3])
    initsigma2 <- exp(tail(res$par, 2)[1])
    #initeta <- exp(res$par[4])
    #initeta <- res$par[4]
    initeta <- tail(res$par, 2)[2]
    #initparms <- c(initBeta, initsigma2, initeta)
    initparms <- c(initBeta, initsigma2, res$par[4])
    rhat <- nleqslv::nleqslv(x = initr, fn = FindrCpp, jac = MyjacCpp,
                             d2 = dt$d2, Mi=Mi.t, delta=delta.t, nrisk = nriskSp, X = X.t,
                             beta=initBeta, eta=initeta, sigma2=initsigma2, control = list(ftol=1e-2))
    initTheta <- c(log(initsigma2), initeta)
    initparms <- c(initBeta, initTheta)
    res <- optim(initparms, l21Cpp,
                 rhat = rhat$x, Mi=Mi.t, delta=delta.t, nrisk = nriskSp, X = X.t, d2=dt$d2,
                 method = "L-BFGS-B", control = list(pgtol = 1e-2))
    updateparms <- res$par
    #diff <- max(abs(updateparms - initparms))
    diff <- max(abs( c(initBeta - matrix(res$par[1:q], ncol = 1),
                       initsigma2 - exp(tail(res$par, 2)[1]),
                       initeta - tail(res$par, 2)[2] ) ))
    #cat("Iter:", iter, "Diff:", diff, "\n")
  }
  finalparms <- c(updateparms[1:q], exp(tail(updateparms, 2)[1]), tail(updateparms, 2)[2] )
  #print(round(finalparms, 4) )
  return(list(parms = finalparms, rhat=rhat$x, iter=iter) )
}
#' Prediction
#'
#' This function computes conformal prediction intervals for spatial survival outcomes.
#'
#' @param train_dt A list object containing the training data, with components
#'   \code{Yij} (survival time), \code{Xij} (covariate matrix),
#'   \code{delta} (event indicator), \code{Mi} (spatial location),
#'   and \code{d2} (distance matrix between locations).
#' @param test_dt A list object containing the testing data, with the same
#'   structure as \code{train_dt}.
#' @param B A positive integer specifying the number of bootstrap samples.
#' @param alpha A numeric value in (0, 1) specifying the miscoverage level,
#'   yielding 100(1 - alpha) percent conformal prediction intervals.
#'
#' @return A matrix with the following columns:
#'   \describe{
#'     \item{Mi}{Spatial location of each test subject.}
#'     \item{RE}{The estimated random effects.}
#'     \item{Up}{Upper bounds of the conformal prediction intervals.}
#'     \item{Up.m}{Modified upper bounds of the conformal prediction intervals.}
#'     \item{Lp}{Lower bounds of the conformal prediction intervals.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Load built-in datasets
#' data(train.dt)
#' data(test.dt)
#'
#' # Compute 95% conformal prediction intervals
#' result <- Predict.Sp(
#'   train_dt = train.dt,
#'   test_dt  = test.dt,
#'   B        = 10,
#'   alpha    = 0.05
#' )
#'
#' # View results
#' head(result)
#' }
#'
#' @export
Predict.Sp <- function(train_dt, test_dt, B = 250, alpha = 0.05){
  fit.Sp <- Findparms(dt = train_dt)
  n.parms <- length(fit.Sp$parms)
  N <- length(train_dt$Yij)
  n1 <- sum(train_dt$delta)
  mass.tto <- GetWeight(dt = train_dt)
  id0 <- which(train_dt$delta==1)
  dat.event.Yij <- train_dt$Yij[id0, ]
  dat.event.Xij <- train_dt$Xij[id0, ]
  dat.event.Mij <- train_dt$Mi[id0]
  beta.1 <- matrix(fit.Sp$parms[1:(n.parms-2)], ncol = 1)
  rhat <- fit.Sp$rhat
  sfun0 <- MyLambda0(dt = train_dt, beta = beta.1, rhat = rhat, Random = T)
  library(parallel)
  library(doSNOW)
  library(foreach)
  cat("Starting parallel bootstrapping...\n")
  ncore <- max(1, parallel::detectCores() - 2)
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  #pb <- txtProgressBar(min = 0, max = B, style = 3)
  progress <- function(n) {
    #  setTxtProgressBar(pb, n)
    cat(sprintf("Bootstrapping progress: %d / %d", n, B), "\r")
    flush.console()
  }
  opts <- list(progress = progress)
  z <- foreach(i = 1:B, .options.snow = opts, 
               .combine = rbind,
               .packages = c("SpaF","survival")
  ) %dopar% {
    while (TRUE) {
      boot.loc <- sample(1:N, size = N, replace = T)
      boot.dat <- list(Yij = train_dt$Yij[boot.loc,],
                       Xij = train_dt$Xij[boot.loc, ],
                       delta = train_dt$delta[boot.loc],
                       ri = NA,
                       Loc = NA,
                       Mi = train_dt$Mi[boot.loc], d2 = train_dt$d2)
      if (min(table(train_dt$Mi[boot.loc]))>=2) break
    }
    boot.cox1 <- Findparms(dt = boot.dat)
    boot.coef <- matrix(boot.cox1$parms[1:(n.parms-2)], ncol = 1)
    boot.sfun0 <- MyLambda0(dt = boot.dat, beta = boot.coef, rhat = boot.cox1$rhat, Random = T)
    ##
    tid <- sample(1:n1, size = 1, replace = T, prob = mass.tto)
    new.x <- dat.event.Xij[tid, ]
    new.y <- dat.event.Yij[tid]
    loc.ri <- dat.event.Mij[tid]
    new.xb.r <- as.numeric( t(boot.coef) %*% new.x + boot.cox1$rhat[loc.ri] )
    z1 <- exp(-boot.sfun0(new.y)*exp(new.xb.r) )
    c(z1)
    ##
    #i^2
  }
  #close(pb)
  stopCluster(cl)
  z1 <- sort(z[,1])
  test.x <- test_dt$Xij
  test.ri <- rhat
  RE <- test.ri[test_dt$Mi]
  test.xb <- as.numeric( t(beta.1) %*% t(test.x)  + test.ri[test_dt$Mi] )
  u <- quantile(z1, 1-alpha/2)
  l <- quantile(z1, alpha/2)
  Up <- -log(l)*exp(-test.xb) / mean(sfun0(train_dt$Yij)/train_dt$Yij)
  Lp <- -log(u)*exp(-test.xb) / mean(sfun0(train_dt$Yij)/train_dt$Yij)
  ## with modify
  id1 <- which.max(dat.event.Yij)
  Y.event.max <- dat.event.Yij[id1]
  X.event <- dat.event.Xij[id1,]
  temp.u <- -log(l)*exp(-test.xb) / mean(sfun0(train_dt$Yij)/train_dt$Yij)
  Up.M <- pmin(Y.event.max, temp.u)
  test.x <- cbind (test.x, Mi=test_dt$Mi)
  pred <- cbind(test.x, RE, Up, Up.M, Lp)
  return(pred)
}
