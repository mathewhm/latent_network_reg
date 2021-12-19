# Code for Submission
# Function Code
########################################### Standard AMEN models with diff. Multi. Effects
# standard AME model but multiplicative effects are U\LambdaV^T where U/V are 
# binary (this is more efficient for sampling multi effects if underlying structure
# is believed to come from community structure)

# most from ame function from 'amen' package, initialU/initialV are initializations of
# U and V which shouold be n x k matrices
library(latex2exp);library(irlba);library(pbivnorm); library(TruncatedNormal)
library(amen); library(dplyr); library(ggplot2); library(reshape); library(bayesplot);
library(R.utils);library(latex2exp)

# functions from AMEN version 1.3 in case 
# there are any version problems
# X: n x n x p arryay
# beta: p x 1 vector
Xbeta <- function (X, beta) 
{
  XB <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
  for (k in seq(1, length(beta), length = length(beta))) {
    XB <- XB + beta[k] * X[, , k]
  }
  XB
}

# Function for simulating Y for censored binary likelihood
#' @param EZ: n x n matrix of expected Z values
#' @param rho: current rho
#' @param odmax: n x1 vector of max censoring
simY_cbin <- function (EZ, rho, odmax, YO = NULL) 
{
  if (length(odmax) == 1) {
    odmax <- rep(odmax, nrow(EZ))
  }
  ZS <- simZ(EZ, rho)
  diag(ZS) <- -Inf
  if (!is.null(YO)) {
    ZS[is.na(YO)] <- -Inf
  }
  YS <- ZS * 0
  for (i in 1:nrow(EZ)) {
    rs <- rank(ZS[i, ]) - (nrow(EZ) - odmax[i])
    YS[i, ] <- rs * (rs > 0) * (ZS[i, ] > 0)
    YS[i, YS[i, ] > 0] <- match(YS[i, YS[i, ] > 0], sort(unique(YS[i, 
                                                                   YS[i, ] > 0])))
  }
  diag(YS) <- NA
  if (!is.null(YO)) {
    YS[is.na(YO)] <- NA
  }
  YS
}


# Function for simulating Z for censored binary likelihood
#' @param Z: n x n matrix of current Z values
#' @param EZ: n x n matrix of expected Z value
#' @param rho: current rho
#' @param Y: n x n network
#' @param odmax: n x1 vector of max censoring
#' @param odobs: out deegree
rZ_cbin_fc <- function (Z, EZ, rho, Y, odmax, odobs) 
{
  sz <- sqrt(1 - rho^2)
  ut <- upper.tri(EZ)
  lt <- lower.tri(EZ)
  for (y in sample(0:1)) {
    if (y == 1) {
      ub <- Inf
      lbm <- matrix(pmax(apply(Z - (Y != 0) * (Inf^(Y != 
                                                      0)), 1, max, na.rm = TRUE), 0), nrow(Z), nrow(Z))
    }
    if (y == 0) {
      lb <- -Inf
      ubm <- matrix(apply(Z + (Y != 1) * (Inf^(Y != 1)), 
                          1, min, na.rm = TRUE), nrow(Z), nrow(Z))
      ubm[odobs < odmax] <- 0
    }
    up <- ut & Y == y
    if (y == 0) {
      ub <- ubm[up]
    }
    if (y == 1) {
      lb <- lbm[up]
    }
    ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
    Z[up] <- ez + sz * qnorm(runif(sum(up), pnorm((lb - 
                                                     ez)/sz), pnorm((ub - ez)/sz)))
    up <- lt & Y == y
    if (y == 0) {
      ub <- ubm[up]
    }
    if (y == 1) {
      lb <- lbm[up]
    }
    ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
    Z[up] <- ez + sz * qnorm(runif(sum(up), pnorm((lb - 
                                                     ez)/sz), pnorm((ub - ez)/sz)))
  }
  Z[is.na(Y)] <- rnorm(sum(is.na(Y)), EZ[is.na(Y)], 1)
  Z
}



# Simulate from Wishart
#' @param S0: pos. definite matrix
#' @param nu: pos. integer
rwish <- function (S0, nu = dim(S0)[1] + 1) 
{
  sS0 <- chol(S0)
  Z <- matrix(rnorm(nu * dim(S0)[1]), nu, dim(S0)[1]) %*% 
    sS0
  t(Z) %*% Z
}


# Updating multi. effects in standard amen way for sym
#' @param E: Residual matrix
#' @param U: Current U
#' @param V: Current v
#' @param s2: dyadic variace
#' @param shrink: shrink with prior or not

rUV_sym_fc<- function (E, U, V, s2 = 1, shrink = TRUE) 
{
  R <- ncol(U)
  n <- nrow(U)
  L <- diag((V[1, ]/U[1, ]), nrow = R)
  L[is.na(L)] <- 1
  if (shrink) {
    ivU <- diag(rgamma(R, (2 + n)/2, (1 + apply(U^2, 2, 
                                                sum))/2), nrow = R)
  }
  if (!shrink) {
    ivU <- diag(1/n, nrow = R)
  }
  for (i in rep(sample(1:n), 4)) {
    l <- L %*% (apply(U * E[i, ], 2, sum) - U[i, ] * E[i, 
                                                       i])/s2
    iQ <- solve((ivU + L %*% (crossprod(U) - U[i, ] %*% 
                                t(U[i, ])) %*% L/s2))
    U[i, ] <- iQ %*% l + t(chol(iQ)) %*% rnorm(R)
  }
  for (r in 1:R) {
    Er <- E - U[, -r, drop = FALSE] %*% L[-r, -r, drop = FALSE] %*% 
      t(U[, -r, drop = FALSE])
    l <- sum((Er * (U[, r] %*% t(U[, r])))[upper.tri(Er)])/s2
    iq <- 1/(1 + sum(((U[, r] %*% t(U[, r]))^2)[upper.tri(Er)])/s2)
    L[r, r] <- rnorm(1, iq * l, sqrt(iq))
  }
  list(U = U, V = U %*% L)
}


# Updating multi. effects in standard amen way assuming
# same across replicates
#' @param E.T: Square residual relational matrix
#' @param U: Current U
#' @param V: Current v
#' @param rho: current rho
#' @param s2: dyadic variace
#' @param shrink: shrink with prior or not
rUV_rep_fc <- function (E.T, U, V, rho, s2 = 1, shrink = TRUE) 
{
  Time <- dim(E.T)[3]
  R <- ncol(U)
  n <- nrow(U)
  UV <- cbind(U, V)
  if (shrink) {
    Suv <- solve(rwish(solve(diag(nrow = 2 * R) + t(UV) %*% 
                               UV), n + R + 2))
  }
  if (!shrink) {
    Suv <- diag(n, nrow = 2 * R)
  }
  Se <- matrix(c(1, rho, rho, 1), 2, 2) * s2
  iSe2 <- mhalf(solve(Se))
  g <- iSe2[1, 1]
  d <- iSe2[1, 2]
  get.Er <- function(E, UVmr) {
    return(E - UVmr)
  }
  get.Es <- function(Er, g, d) {
    n <- sqrt(length(Er))
    return((g^2 + d^2) * matrix(Er, n) + 2 * g * d * matrix(Er, 
                                                            n, byrow = T))
  }
  for (r in sample(1:R)) {
    UVmr <- tcrossprod(U[, -r], V[, -r])
    Er.t <- apply(E.T, 3, get.Er, UVmr)
    Es.t <- apply(Er.t, 2, get.Es, g, d)
    vr <- V[, r]
    b0 <- c(Suv[r, -r] %*% solve(Suv[-r, -r]))
    v0 <- c(Suv[r, r] - b0 %*% Suv[-r, r])
    m0 <- cbind(U[, -r], V) %*% b0
    ssv <- max(sum(vr^2), 1e-06)
    a <- Time * (g^2 + d^2) * ssv + 1/v0
    c <- -2 * Time * g * d/(a^2 + a * 2 * Time * g * d * 
                              ssv)
    Esv.vec <- rowSums(Es.t)
    nEsv <- sqrt(length(Esv.vec))
    Esv <- matrix(Esv.vec, nEsv) %*% vr
    m1 <- Esv/a + c * vr * sum((Esv + m0/v0) * vr) + m0/(a * 
                                                           v0)
    ah <- sqrt(1/a)
    bh <- (sqrt(1/a + ssv * c) - sqrt(1/a))/ssv
    e <- rnorm(nrow(E.T[, , 1]))
    U[, r] <- m1 + ah * e + bh * vr * sum(vr * e)
    ur <- U[, r]
    rv <- R + r
    b0 <- c(Suv[rv, -rv] %*% solve(Suv[-rv, -rv]))
    v0 <- c(Suv[rv, rv] - b0 %*% Suv[-rv, rv])
    m0 <- cbind(U, V[, -r]) %*% b0
    ssu <- max(sum(ur^2), 1e-06)
    a <- Time * (g^2 + d^2) * ssu + 1/v0
    c <- -2 * Time * g * d/(a^2 + a * 2 * Time * g * d * 
                              ssu)
    tEsu <- matrix(Esv.vec, nEsv, byrow = T) %*% ur
    m1 <- tEsu/a + c * ur * sum((tEsu + m0/v0) * ur) + m0/(a * 
                                                             v0)
    ah <- sqrt(1/a)
    bh <- (sqrt(1/a + ssu * c) - sqrt(1/a))/ssu
    e <- rnorm(nrow(E.T[, , 1]))
    V[, r] <- m1 + ah * e + bh * ur * sum(ur * e)
  }
  list(U = U, V = V)
}

# Updating multi. effects in standard amen way 
#' @param E: Residual matrix
#' @param U: Current U
#' @param V: Current v
#' @param s2: dyadic variace
#' @param shrink: shrink with prior or not

rUV_fc <- function (E, U, V, rho, s2 = 1, shrink = TRUE) 
{
  R <- ncol(U)
  n <- nrow(U)
  UV <- cbind(U, V)
  if (shrink) {
    Suv <- solve(rwish(solve(diag(nrow = 2 * R) + t(UV) %*% 
                               UV), n + R + 2))
  }
  if (!shrink) {
    Suv <- diag(n, nrow = 2 * R)
  }
  Se <- matrix(c(1, rho, rho, 1), 2, 2) * s2
  iSe2 <- mhalf(solve(Se))
  g <- iSe2[1, 1]
  d <- iSe2[1, 2]
  for (r in sample(1:R)) {
    Er <- E - U[, -r] %*% t(V[, -r])
    Es <- (g^2 + d^2) * Er + 2 * g * d * t(Er)
    vr <- V[, r]
    b0 <- c(Suv[r, -r] %*% solve(Suv[-r, -r]))
    v0 <- c(Suv[r, r] - b0 %*% Suv[-r, r])
    m0 <- cbind(U[, -r], V) %*% b0
    ssv <- max(sum(vr^2), 1e-06)
    a <- (g^2 + d^2) * ssv + 1/v0
    c <- -2 * g * d/(a^2 + a * 2 * g * d * ssv)
    Esv <- Es %*% vr
    m1 <- Esv/a + c * vr * sum((Esv + m0/v0) * vr) + m0/(a * 
                                                           v0)
    ah <- sqrt(1/a)
    bh <- (sqrt(1/a + ssv * c) - sqrt(1/a))/ssv
    e <- rnorm(nrow(E))
    U[, r] <- m1 + ah * e + bh * vr * sum(vr * e)
    ur <- U[, r]
    rv <- R + r
    b0 <- c(Suv[rv, -rv] %*% solve(Suv[-rv, -rv]))
    v0 <- c(Suv[rv, rv] - b0 %*% Suv[-rv, rv])
    m0 <- cbind(U, V[, -r]) %*% b0
    ssu <- max(sum(ur^2), 1e-06)
    a <- (g^2 + d^2) * ssu + 1/v0
    c <- -2 * g * d/(a^2 + a * 2 * g * d * ssu)
    tEsu <- t(Es) %*% ur
    m1 <- tEsu/a + c * ur * sum((tEsu + m0/v0) * ur) + m0/(a * 
                                                             v0)
    ah <- sqrt(1/a)
    bh <- (sqrt(1/a + ssu * c) - sqrt(1/a))/ssu
    e <- rnorm(nrow(E))
    V[, r] <- m1 + ah * e + bh * ur * sum(ur * e)
  }
  list(U = U, V = V)
}




myAME_estimateLambda <-function (Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL, rvar = !(model == 
                                                                                       "rrl"), cvar = TRUE, dcor = !symmetric, nvar = TRUE, R = 0, 
                                 model = "nrm", intercept = !is.element(model, c("rrl", "ord")), 
                                 symmetric = FALSE, odmax = rep(max(apply(Y > 0, 1, sum, na.rm = TRUE)), 
                                                                nrow(Y)), seed = 1, nscan = 10000, burn = 500, odens = 25, 
                                 plot = TRUE, print = TRUE, gof = TRUE, initialU, initialV) 
{
  set.seed(seed)
  Y_binary <- 1 * (Y > 0)
  
  counter <- 0
  MH <-0
  
  U <-  initialU
  V <-  initialV
  prior.lambda.mu  <- rep(1, dim(U)[2]^2)
  prior.lambda.var <- prior.lambda.prec <- diag(2, dim(U)[2]^2) ## have had .2
  
  #prior.lambda.prec <- 1/prior.lambda.var
  
  lambda <- matrix(prior.lambda.mu, nrow = dim(U)[2])
  uvl <- U %*% lambda %*% t(V)
  
  
  diag(Y) <- NA
  Y_binary <- Y
  diag(Y_binary) <- 0
  if (is.element(model, c("bin", "cbin"))) {
    Y <- 1 * (Y > 0)
  }
  if (is.element(model, c("cbin", "frn", "rrl"))) {
    odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
    if (length(odmax) == 1) {
      odmax <- rep(odmax, nrow(Y))
    }
  }
  if (symmetric) {
    Xcol <- Xrow
    rvar <- cvar <- nvar
  }
  n <- nrow(Y)
  pr <- length(Xrow)/n
  pc <- length(Xcol)/n
  pd <- length(Xdyad)/n^2
  X <- design_array(Xrow, Xcol, Xdyad, intercept, nrow(Y))
  if (model == "rrl" & any(apply(apply(X, c(1, 3), var), 2, 
                                 sum) == 0) & !any(apply(X, 3, function(x) {
                                   var(c(x))
                                 }) == 0)) {
    cat("WARNING: row effects are not estimable using this procedure ", 
        "\n")
  }
  if (is.element(model, c("ord", "rrl")) & any(apply(X, 3, 
                                                     function(x) {
                                                       var(c(x))
                                                     }) == 0)) {
    cat("WARNING: an intercept is not estimable using this procedure ", 
        "\n")
  }
  if (is.element(model, c("frn", "rrl"))) {
    ymx <- max(apply(1 * (Y > 0), 1, sum, na.rm = TRUE))
    YL <- NULL
    warn <- FALSE
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(yi[rnkd]) > length(unique(yi[rnkd]))) {
        warn <- TRUE
      }
      yi[rnkd] <- rank(yi[rnkd], ties.method = "random")
      Y[i, ] <- yi
      YL <- rbind(YL, match(1:ymx, yi))
    }
    if (warn) {
      cat("WARNING: Random reordering used to break ties in ranks\n")
    }
  }
  if (model == "nrm") {
    Z <- Y
  }
  if (model == "ord") {
    Z <- matrix(zscores(Y), nrow(Y), ncol(Y))
  }
  if (model == "rrl") {
    Z <- matrix(t(apply(Y, 1, zscores)), nrow(Y), ncol(Y))
  }
  if (model == "bin") {
    Z <- matrix(zscores(Y), nrow(Y), nrow(Y))
    z01 <- 0.5 * (max(Z[Y == 0], na.rm = TRUE) + min(Z[Y == 
                                                         1], na.rm = TRUE))
    Z <- Z - z01
  }
  if (is.element(model, c("cbin", "frn"))) {
    Z <- Y
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      zi <- zscores(yi)
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
        zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 0.001
      }
      if (length(rnkd) < odmax[i]) {
        urnkd <- which(!is.na(yi) & yi == 0)
        if (max(zi[urnkd]) > 0) {
          zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 0.001
        }
      }
      Z[i, ] <- zi
    }
  }
  mu <- mean(Z, na.rm = TRUE)
  a <- rowMeans(Z, na.rm = TRUE)
  b <- colMeans(Z, na.rm = TRUE)
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  ZA <- mu + outer(a, b, "+")
  Z[is.na(Z)] <- ZA[is.na(Z)]
  beta <- rep(2, dim(X)[3])
  s2 <- 1
  rho <- 0
  Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar))
  # U <- V <- matrix(0, nrow(Y), R)
  BETA <- matrix(nrow = 0, ncol = dim(X)[3] - pr * symmetric)
  LAMBDA <- matrix(nrow = 0, ncol = dim(U)[2]^2)
  U_vect <- matrix(nrow = 0, ncol = dim(U)[1]*dim(U)[2])
  
  
  VC <- matrix(nrow = 0, ncol = 5 - 3 * symmetric)
  UVPS <- U %*% t(V) * 0
  APS <- BPS <- rep(0, nrow(Y))
  YPS <- matrix(0, nrow(Y), ncol(Y))
  dimnames(YPS) <- dimnames(Y)
  GOF <- matrix(gofstats(Y), 1, 4)
  rownames(GOF) <- "obs"
  colnames(GOF) <- c("sd.rowmean", "sd.colmean", "dyad.dep", 
                     "triad.dep")
  names(APS) <- names(BPS) <- rownames(U) <- rownames(V) <- rownames(Y)
  if (!symmetric) {
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
    colnames(BETA) <- dimnames(X)[[3]]
  }
  if (symmetric) {
    colnames(VC) <- c("va", "ve")
    rb <- intercept + seq(1, pr, length = pr)
    cb <- intercept + pr + seq(1, pr, length = pr)
    bnames <- dimnames(X)[[3]]
    bni <- bnames[1 * intercept]
    bnn <- gsub("row", bnames[rb], replacement = "node")
    bnd <- bnames[-c(1 * intercept, rb, cb)]
    colnames(BETA) <- c(bni, bnn, bnd)
  }
  
  Xr <- apply(X, c(1, 3), sum)
  Xc <- apply(X, c(2, 3), sum)
  mX <- apply(X, 3, c)
  mXt <- apply(aperm(X, c(2, 1, 3)), 3, c)
  XX <- t(mX) %*% mX
  XXt <- t(mX) %*% mXt
  have_coda <- suppressWarnings(try(requireNamespace("coda", 
                                                     quietly = TRUE), silent = TRUE))
  print("start")
  for (s in 1:(nscan + burn)) {
    EZ <- Xbeta(X, beta) + outer(a, b, "+") + U %*% lambda %*% t(V)
    if (model == "nrm") {
      Z <- rZ_nrm_fc(Z, EZ, rho, s2, Y)
    }
    if (model == "bin") {
      Z <- rZ_bin_fc(Z, EZ, rho, Y)
    }
    if (model == "ord") {
      Z <- rZ_ord_fc(Z, EZ, rho, Y)
    }
    if (model == "cbin") {
      Z <- rZ_cbin_fc(Z, EZ, rho, Y, odmax, odobs)
    }
    if (model == "frn") {
      Z <- rZ_frn_fc(Z, EZ, rho, Y, YL, odmax, odobs)
    }
    if (model == "rrl") {
      Z <- rZ_rrl_fc(Z, EZ, rho, Y, YL)
    }
    if (model == "nrm") {
      s2 <- rs2_fc(Z - EZ, rho)
      
    }
    
    
    print("did z")
    
    tmp <- rbeta_ab_fc(Z - uvl, Sab, rho, X, mX, mXt, 
                       XX, XXt, Xr, Xc, s2)
    print("Ab")
    beta <- tmp$beta
    a <- tmp$a * rvar
    b <- tmp$b * cvar
    
    if (symmetric) {
      a <- b <- (a + b)/2
    }
    
    ## AS LONG AS WE STAY SYMMETRIc WE GOOD
    if (rvar & cvar & !symmetric) {
      if (is.element(model, c("nrm", "ord"))) {
        Sab <- solve(rwish(solve(diag(2) + crossprod(cbind(a, 
                                                           b))), 3 + nrow(Z)))
      }
      if (model == "bin") {
        tmp <- raSab_bin_fc_new(Z, Y, a, b, Sab)
        print("sab")
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
      if (model == "cbin") {
        tmp <- raSab_cbin_fc(Z, Y, a, b, Sab, odmax, 
                             odobs)
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
      if (model == "frn") {
        tmp <- raSab_frn_fc(Z, Y, YL, a, b, Sab, odmax, 
                            odobs)
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
    }
    if (rvar & !cvar & !symmetric) {
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(a^2))/2)
    }
    if (!rvar & cvar & !symmetric) {
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(b^2))/2)
    }
    if (symmetric & nvar) {
      Sab[1, 1] <- Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, 
                                         (1 + sum(a^2))/2)
      Sab[1, 2] <- Sab[2, 1] <- 0.999 * Sab[1, 1]
    }
    if (dcor) {
      rho <- rrho_mh(Z - (Xbeta(X, beta) + outer(a, b, 
                                                 "+") + uvl), rho, s2)
    }
    if (symmetric) {
      rho <- min(0.9999, 1 - 1/sqrt(s))
    }
    Se   <-matrix(c(1,rho,rho,1),2,2)*1
    iSe2 <-mhalf(solve(Se))
    td   <-iSe2[1,1] ; to<-iSe2[1,2]
    XB <-  Xbeta(X, beta)
    lambda <- drawLambda(Z, XB + outer(a,b, "+"), U,U, prior.lambda.var,
                         prior.lambda.mu, td, to, TRUE)
    if (R > 0) {
      E <- Z - (Xbeta(X, beta) + outer(a, b, "+"))
      if (symmetric) {
        E <- 0.5 * (E + t(E))
      }
      shrink <- (s > 0.5 * burn)
      
      
      if (symmetric) {
        XB <- Xbeta(X, beta)
        # UV <- myrUV_sym_fc_membership(lambda, U, V, s2, shrink, X,beta, Y,a,b)
        # UV <- drawU(lambda, U, U,  XB + outer(a,b, "+"), Y_binary, symmetric, counter, MH)
        # U_info       <- drawU(lambda, U, U,  XB + outer(a,b, "+"), Y_binary, symmetric,counter, MH)
        # U            <- U_info$U
        
        UV_info      <- drawUV_probitlLikelihood_amen_lambda(lambda, U, U,  XB ,a,b,
                                                             Y_binary, UandV=FALSE,counter, MH,
                                                             h_r_past, h_c_past,Z,td,to,Sab, rho)
        
        
        
      }
      if (!symmetric) {
        XB <- Xbeta(X, beta)
        #UV <- myrUV_fc_membership(lambda, U, V,  s2, shrink, X, beta, Y,a,b)
        # UV <- drawU(lambda, U, U,  XB + outer(a,b, "+"), Y_binary, symmetric, counter, MH)
        
        
        # U_info       <- drawU(lambda, U, U,  XB + outer(a,b, "+"), Y_binary, symmetric,counter, MH)
        UV_info      <- drawUV_probitlLikelihood_amen_lambda(lambda, U, U,  XB ,a,b,
                                                             Y_binary, UandV=FALSE,counter, MH,
                                                             h_r_past, h_c_past,Z,td,to,Sab, rho)
        
        
      }
      U <- UV_info$U
      V <- UV_info$U
      updatedU <- UV_info$updateU
      swapPeep     <- UV_info$swapPeep
      counter      <- UV_info$counter
      MH           <- UV_info$MH
      
      print(beta)
      
      
      
      uvl <- U %*% lambda %*% t(V)
      
      EZ <- XB + uvl + outer(a, b, "+")
      
      if(updatedU == TRUE){
        for(h in 1:length(swapPeep)){
          Z <-     update_Z_bin(swapPeep[h], Y_binary, EZ, rho,  n,Z)
        }
        
      }
      
      
      #membership <- UV_info$membership
      
    }
    if (s%%odens == 0 & s <= burn & print) {
      cat(round(100 * s/burn, 2), " pct burnin complete \n")
    }
    if (s%%odens == 0 & s > burn) {
      if (symmetric) {
        br <- beta[rb]
        bc <- beta[cb]
        bn <- (br + bc)/2
        sbeta <- c(beta[1 * intercept], bn, beta[-c(1 * 
                                                      intercept, rb, cb)])
        BETA <- rbind(BETA, sbeta)
        LAMBDA <- rbind(LAMBDA, c(lambda))
        VC <- rbind(VC, c(Sab[1, 1], s2))
        U_vect <- rbind(U_vect, c(U))
      }
      if (!symmetric) {
        BETA <- rbind(BETA, beta)
        
        LAMBDA <- rbind(LAMBDA, c(lambda))
        VC <- rbind(VC, c(Sab[upper.tri(Sab, diag = T)], 
                          rho, s2))
        U_vect <- rbind(U_vect, c(U))
        
      }
      UVPS <- UVPS + uvl
      APS <- APS + a
      BPS <- BPS + b
      EZ <- Xbeta(X, beta) + outer(a, b, "+") + uvl
      if (symmetric) {
        EZ <- (EZ + t(EZ))/2
      }
      if (model == "bin") {
        Ys <- simY_bin(EZ, rho)
      }
      if (model == "cbin") {
        Ys <- 1 * (simY_frn(EZ, rho, odmax, YO = Y) > 
                     0)
      }
      if (model == "frn") {
        Ys <- simY_frn(EZ, rho, odmax, YO = Y)
      }
      if (model == "rrl") {
        Ys <- simY_rrl(EZ, rho, odobs, YO = Y)
      }
      if (model == "nrm") {
        Ys <- simY_nrm(EZ, rho, s2)
      }
      if (model == "ord") {
        Ys <- simY_ord(EZ, rho, Y)
      }
      if (symmetric) {
        Ys[lower.tri(Ys)] <- 0
        Ys <- Ys + t(Ys)
      }
      YPS <- YPS + Ys
      if (gof) {
        Ys[is.na(Y)] <- NA
        GOF <- rbind(GOF, gofstats(Ys))
      }
      if (print) {
        cat(s, round(apply(BETA, 2, mean), 2), ":", round(apply(VC, 
                                                                2, mean), 2), "\n")
        if (have_coda & nrow(VC) > 3 & length(beta) > 
            0) {
          cat(round(coda::effectiveSize(BETA)), "\n")
        }
      }
      if (plot) {
        par(mfrow = c(1 + 2 * gof, 2), mar = c(3, 3, 
                                               1, 1), mgp = c(1.75, 0.75, 0))
        mVC <- apply(VC, 2, median)
        matplot(VC, type = "l", lty = 1)
        abline(h = mVC, col = 1:length(mVC))
        if (length(beta) > 0) {
          mBETA <- apply(BETA, 2, median)
          matplot(BETA, type = "l", lty = 1, col = 1:length(mBETA))
          abline(h = mBETA, col = 1:length(mBETA))
          abline(h = 0, col = "gray")
        }
        if (gof) {
          for (k in 1:4) {
            hist(GOF[-1, k], xlim = range(GOF[, k]), 
                 main = "", prob = TRUE, xlab = colnames(GOF)[k], 
                 col = "lightblue", ylab = "", yaxt = "n")
            abline(v = GOF[1, k], col = "red")
          }
        }
      }
    }
  }
  APM <- APS/nrow(VC)
  BPM <- BPS/nrow(VC)
  UVPM <- UVPS/nrow(VC)
  YPM <- YPS/nrow(VC)
  EZ <- Xbeta(X, apply(BETA, 2, mean)) + outer(APM, BPM, "+") + 
    UVPM
  names(APM) <- names(BPM) <- rownames(UVPM) <- colnames(UVPM) <- dimnames(Y)[[1]]
  dimnames(YPM) <- dimnames(EZ) <- dimnames(Y)
  rownames(BETA) <- NULL
  rownames(LAMBDA) <- NULL
  if (!symmetric) {
    # UDV <- svd(UVPM)
    # U <- UDV$u[, seq(1, R, length = R)] %*% diag(sqrt(UDV$d)[seq(1, 
    #                                                              R, length = R)], nrow = R)
    # V <- UDV$v[, seq(1, R, length = R)] %*% diag(sqrt(UDV$d)[seq(1, 
    #                                                              R, length = R)], nrow = R)
    rownames(U) <- rownames(V) <- rownames(Y)
    fit <- list(BETA = BETA, VC = VC, APM = APM, BPM = BPM, 
                U = U_vect, V = V, UVPM = UVPM, EZ = EZ, YPM = YPM, GOF = GOF, LAMBDA = LAMBDA)
  }
  if (symmetric) {
    ULUPM <- UVPM
    #  eULU <- eigen(ULUPM)
    # eR <- which(rank(-abs(eULU$val), ties.method = "first") <= 
    #               R)
    #  U <- eULU$vec[, seq(1, R, length = R), drop = FALSE]
    # L <- eULU$val[eR]
    rownames(U) <- rownames(ULUPM) <- colnames(ULUPM) <- rownames(Y)
    EZ <- 0.5 * (EZ + t(EZ))
    YPM <- 0.5 * (YPM + t(YPM))
    fit <- list(BETA = BETA, VC = VC, APM = APM, U = U_vect, V = V, 
                ULUPM = ULUPM, EZ = EZ, YPM = YPM, GOF = GOF, LAMBDA = LAMBDA)
  }
  class(fit) <- "ame"
  fit
}


# AME function that allows to estimate multiplicative effects where U is binary and 
# written as UU^T
# most from ame function from 'amen' package, initialU/initialV are initializations of
# U and V which shouold be n x k

myAME <-function (Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL, rvar = !(model == 
                                                                        "rrl"), cvar = TRUE, dcor = !symmetric, nvar = TRUE, R = 0, 
                  model = "nrm", intercept = !is.element(model, c("rrl", "ord")), 
                  symmetric = FALSE, odmax = rep(max(apply(Y > 0, 1, sum, na.rm = TRUE)), 
                                                 nrow(Y)), seed = 1, nscan = 10000, burn = 500, odens = 25, 
                  plot = TRUE, print = TRUE, gof = TRUE, initialU, initialV) 
{
  Y_binary <- 1 * (Y > 0)
  counter <- 0
  MH <-0
  
  set.seed(seed)
  U <-  initialU
  V <-  initialV
  k <-dim(U)[2]
  
  diag(Y) <- NA
  if (is.element(model, c("bin", "cbin"))) {
    Y <- 1 * (Y > 0)
  }
  if (is.element(model, c("cbin", "frn", "rrl"))) {
    odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
    if (length(odmax) == 1) {
      odmax <- rep(odmax, nrow(Y))
    }
  }
  if (symmetric) {
    Xcol <- Xrow
    rvar <- cvar <- nvar
  }
  n <- nrow(Y)
  pr <- length(Xrow)/n
  pc <- length(Xcol)/n
  pd <- length(Xdyad)/n^2
  X <- design_array(Xrow, Xcol, Xdyad, intercept, nrow(Y))
  if (model == "rrl" & any(apply(apply(X, c(1, 3), var), 2, 
                                 sum) == 0) & !any(apply(X, 3, function(x) {
                                   var(c(x))
                                 }) == 0)) {
    cat("WARNING: row effects are not estimable using this procedure ", 
        "\n")
  }
  if (is.element(model, c("ord", "rrl")) & any(apply(X, 3, 
                                                     function(x) {
                                                       var(c(x))
                                                     }) == 0)) {
    cat("WARNING: an intercept is not estimable using this procedure ", 
        "\n")
  }
  if (is.element(model, c("frn", "rrl"))) {
    ymx <- max(apply(1 * (Y > 0), 1, sum, na.rm = TRUE))
    YL <- NULL
    warn <- FALSE
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(yi[rnkd]) > length(unique(yi[rnkd]))) {
        warn <- TRUE
      }
      yi[rnkd] <- rank(yi[rnkd], ties.method = "random")
      Y[i, ] <- yi
      YL <- rbind(YL, match(1:ymx, yi))
    }
    if (warn) {
      cat("WARNING: Random reordering used to break ties in ranks\n")
    }
  }
  if (model == "nrm") {
    Z <- Y
  }
  if (model == "ord") {
    Z <- matrix(zscores(Y), nrow(Y), ncol(Y))
  }
  if (model == "rrl") {
    Z <- matrix(t(apply(Y, 1, zscores)), nrow(Y), ncol(Y))
  }
  if (model == "bin") {
    Z <- matrix(zscores(Y), nrow(Y), nrow(Y))
    z01 <- 0.5 * (max(Z[Y == 0], na.rm = TRUE) + min(Z[Y == 
                                                         1], na.rm = TRUE))
    Z <- Z - z01
  }
  if (is.element(model, c("cbin", "frn"))) {
    Z <- Y
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      zi <- zscores(yi)
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
        zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 0.001
      }
      if (length(rnkd) < odmax[i]) {
        urnkd <- which(!is.na(yi) & yi == 0)
        if (max(zi[urnkd]) > 0) {
          zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 0.001
        }
      }
      Z[i, ] <- zi
    }
  }
  mu <- mean(Z, na.rm = TRUE)
  a <- rowMeans(Z, na.rm = TRUE)
  b <- colMeans(Z, na.rm = TRUE)
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  ZA <- mu + outer(a, b, "+")
  Z[is.na(Z)] <- ZA[is.na(Z)]
  beta <- rep(0, dim(X)[3])
  s2 <- 1
  rho <- 0
  Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar))
  # U <- V <- matrix(0, nrow(Y), R)
  U_vect <- matrix(nrow = 0, ncol = dim(U)[1]*dim(U)[2])
  BETA <- matrix(nrow = 0, ncol = dim(X)[3] - pr * symmetric)
  VC <- matrix(nrow = 0, ncol = 5 - 3 * symmetric)
  UVPS <- U %*% t(V) * 0
  APS <- BPS <- rep(0, nrow(Y))
  YPS <- matrix(0, nrow(Y), ncol(Y))
  dimnames(YPS) <- dimnames(Y)
  GOF <- matrix(gofstats(Y), 1, 4)
  rownames(GOF) <- "obs"
  colnames(GOF) <- c("sd.rowmean", "sd.colmean", "dyad.dep", 
                     "triad.dep")
  names(APS) <- names(BPS) <- rownames(U) <- rownames(V) <- rownames(Y)
  if (!symmetric) {
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
    colnames(BETA) <- dimnames(X)[[3]]
  }
  if (symmetric) {
    colnames(VC) <- c("va", "ve")
    rb <- intercept + seq(1, pr, length = pr)
    cb <- intercept + pr + seq(1, pr, length = pr)
    bnames <- dimnames(X)[[3]]
    bni <- bnames[1 * intercept]
    bnn <- gsub("row", bnames[rb], replacement = "node")
    bnd <- bnames[-c(1 * intercept, rb, cb)]
    colnames(BETA) <- c(bni, bnn, bnd)
  }
  Xr <- apply(X, c(1, 3), sum)
  Xc <- apply(X, c(2, 3), sum)
  mX <- apply(X, 3, c)
  mXt <- apply(aperm(X, c(2, 1, 3)), 3, c)
  XX <- t(mX) %*% mX
  XXt <- t(mX) %*% mXt
  lambda <- diag(k)
  have_coda <- suppressWarnings(try(requireNamespace("coda", 
                                                     quietly = TRUE), silent = TRUE))
  
  
  diag(Y_binary) <- 0
  for (s in 1:(nscan + burn)) {
    EZ <- Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)
    if (model == "nrm") {
      Z <- rZ_nrm_fc(Z, EZ, rho, s2, Y)
    }
    if (model == "bin") {
      Z <- rZ_bin_fc(Z, EZ, rho, Y)
    }
    if (model == "ord") {
      Z <- rZ_ord_fc(Z, EZ, rho, Y)
    }
    if (model == "cbin") {
      Z <- rZ_cbin_fc(Z, EZ, rho, Y, odmax, odobs)
    }
    if (model == "frn") {
      Z <- rZ_frn_fc(Z, EZ, rho, Y, YL, odmax, odobs)
    }
    if (model == "rrl") {
      Z <- rZ_rrl_fc(Z, EZ, rho, Y, YL)
    }
    if (model == "nrm") 
      s2 <- rs2_fc(Z - EZ, rho)
    
    
    tmp <- rbeta_ab_fc(Z - U %*% t(V), Sab, rho, X, mX, mXt, 
                       XX, XXt, Xr, Xc, s2)
    beta <- tmp$beta
    XB <- 
      a <- tmp$a * rvar
    b <- tmp$b * cvar
    if (symmetric) {
      a <- b <- (a + b)/2
    }
    if (rvar & cvar & !symmetric) {
      if (is.element(model, c("nrm", "ord"))) {
        Sab <- solve(rwish(solve(diag(2) + crossprod(cbind(a, 
                                                           b))), 3 + nrow(Z)))
      }
      if (model == "bin") {
        tmp <- raSab_bin_fc(Z, Y, a, b, Sab)
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
      if (model == "cbin") {
        tmp <- raSab_cbin_fc(Z, Y, a, b, Sab, odmax, 
                             odobs)
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
      if (model == "frn") {
        tmp <- raSab_frn_fc(Z, Y, YL, a, b, Sab, odmax, 
                            odobs)
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
    }
    if (rvar & !cvar & !symmetric) {
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(a^2))/2)
    }
    if (!rvar & cvar & !symmetric) {
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(b^2))/2)
    }
    if (symmetric & nvar) {
      Sab[1, 1] <- Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, 
                                         (1 + sum(a^2))/2)
      Sab[1, 2] <- Sab[2, 1] <- 0.999 * Sab[1, 1]
    }
    if (dcor) {
      rho <- rrho_mh(Z - (Xbeta(X, beta) + outer(a, b, 
                                                 "+") + U %*% t(V)), rho, s2)
    }
    
    XB <- Xbeta(X, beta)
    if (symmetric) {
      rho <- min(0.9999, 1 - 1/sqrt(s))
    }
    Se   <-matrix(c(1,rho,rho,1),2,2)*1
    iSe2 <-mhalf(solve(Se))
    td   <-iSe2[1,1] ; to<-iSe2[1,2]
    
    if (R > 0) {
      E <- Z - (Xbeta(X, beta) + outer(a, b, "+"))
      if (symmetric) {
        E <- 0.5 * (E + t(E))
      }
      shrink <- (s > 0.5 * burn)
      if (symmetric) {
        UV_info      <- drawUV_probitlLikelihood_amen_lambda(lambda, U, U,  XB ,a,b,
                                                             Y_binary, UandV=FALSE,counter, MH,
                                                             h_r_past, h_c_past,Z,td,to,Sab, rho)
        
        
        U            <- UV_info$U
        updatedU     <- UV_info$updateU
        counter      <- UV_info$counter
        MH           <- UV_info$MH
      }
      if (!symmetric) {
        #UV <- drawU(lambda, U, U,  XB + outer(a,b, "+"), Y_binary, symmetric)
        UV_info      <- drawUV_probitlLikelihood_amen_lambda(lambda, U, U,  XB ,a,b,
                                                             Y_binary, UandV=FALSE,counter, MH,
                                                             h_r_past, h_c_past,Z,td,to,Sab, rho)
        U            <- UV_info$U
        updatedU     <- UV_info$updateU
        counter      <- UV_info$counter
        MH           <- UV_info$MH
      }
      U <- V <- UV_info$U
      
      updatedU <- UV_info$updateU
      swapPeep     <- UV_info$swapPeep
      counter      <- UV_info$counter
      MH           <- UV_info$MH
      
      
      
      
      
      
      EZ <- XB + U%*%t(U) + outer(a, b, "+")
      
      if(updatedU == TRUE){
        for(h in 1:length(swapPeep)){
          Z <-     update_Z_bin(swapPeep[h], Y_binary, EZ, rho,  n,Z)
        }
        
      }
      
    }
    if (s%%odens == 0 & s <= burn & print) {
      cat(round(100 * s/burn, 2), " pct burnin complete \n")
    }
    if (s%%odens == 0 & s > burn) {
      if (symmetric) {
        br <- beta[rb]
        bc <- beta[cb]
        bn <- (br + bc)/2
        sbeta <- c(beta[1 * intercept], bn, beta[-c(1 * 
                                                      intercept, rb, cb)])
        BETA <- rbind(BETA, sbeta)
        VC <- rbind(VC, c(Sab[1, 1], s2))
        U_vect <- rbind(U_vect, c(U))
      }
      if (!symmetric) {
        BETA <- rbind(BETA, beta)
        VC <- rbind(VC, c(Sab[upper.tri(Sab, diag = T)], 
                          rho, s2))
        U_vect <- rbind(U_vect, c(U))
        
      }
      UVPS <- UVPS + U %*% t(V)
      APS <- APS + a
      BPS <- BPS + b
      EZ <- Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)
      if (symmetric) {
        EZ <- (EZ + t(EZ))/2
      }
      if (model == "bin") {
        Ys <- simY_bin(EZ, rho)
      }
      if (model == "cbin") {
        Ys <- 1 * (simY_frn(EZ, rho, odmax, YO = Y) > 
                     0)
      }
      if (model == "frn") {
        Ys <- simY_frn(EZ, rho, odmax, YO = Y)
      }
      if (model == "rrl") {
        Ys <- simY_rrl(EZ, rho, odobs, YO = Y)
      }
      if (model == "nrm") {
        Ys <- simY_nrm(EZ, rho, s2)
      }
      if (model == "ord") {
        Ys <- simY_ord(EZ, rho, Y)
      }
      if (symmetric) {
        Ys[lower.tri(Ys)] <- 0
        Ys <- Ys + t(Ys)
      }
      YPS <- YPS + Ys
      if (gof) {
        Ys[is.na(Y)] <- NA
        GOF <- rbind(GOF, gofstats(Ys))
      }
      if (print) {
        cat(s, round(apply(BETA, 2, mean), 2), ":", round(apply(VC, 
                                                                2, mean), 2), "\n")
        if (have_coda & nrow(VC) > 3 & length(beta) > 
            0) {
          cat(round(coda::effectiveSize(BETA)), "\n")
        }
      }
      if (plot) {
        par(mfrow = c(1 + 2 * gof, 2), mar = c(3, 3, 
                                               1, 1), mgp = c(1.75, 0.75, 0))
        mVC <- apply(VC, 2, median)
        matplot(VC, type = "l", lty = 1)
        abline(h = mVC, col = 1:length(mVC))
        if (length(beta) > 0) {
          mBETA <- apply(BETA, 2, median)
          matplot(BETA, type = "l", lty = 1, col = 1:length(mBETA))
          abline(h = mBETA, col = 1:length(mBETA))
          abline(h = 0, col = "gray")
        }
        if (gof) {
          for (k in 1:4) {
            hist(GOF[-1, k], xlim = range(GOF[, k]), 
                 main = "", prob = TRUE, xlab = colnames(GOF)[k], 
                 col = "lightblue", ylab = "", yaxt = "n")
            abline(v = GOF[1, k], col = "red")
          }
        }
      }
    }
  }
  APM <- APS/nrow(VC)
  BPM <- BPS/nrow(VC)
  UVPM <- UVPS/nrow(VC)
  YPM <- YPS/nrow(VC)
  EZ <- Xbeta(X, apply(BETA, 2, mean)) + outer(APM, BPM, "+") + 
    UVPM
  names(APM) <- names(BPM) <- rownames(UVPM) <- colnames(UVPM) <- dimnames(Y)[[1]]
  dimnames(YPM) <- dimnames(EZ) <- dimnames(Y)
  rownames(BETA) <- NULL
  if (!symmetric) {
    # UDV <- svd(UVPM)
    # U <- UDV$u[, seq(1, R, length = R)] %*% diag(sqrt(UDV$d)[seq(1, 
    #                                                              R, length = R)], nrow = R)
    # V <- UDV$v[, seq(1, R, length = R)] %*% diag(sqrt(UDV$d)[seq(1, 
    #                                                              R, length = R)], nrow = R)
    rownames(U) <- rownames(V) <- rownames(Y)
    fit <- list(BETA = BETA, VC = VC, APM = APM, BPM = BPM, 
                U = U, V = V, UVPM = UVPM, EZ = EZ, YPM = YPM, GOF = GOF)
  }
  if (symmetric) {
    ULUPM <- UVPM
    #  eULU <- eigen(ULUPM)
    # eR <- which(rank(-abs(eULU$val), ties.method = "first") <= 
    #               R)
    #  U <- eULU$vec[, seq(1, R, length = R), drop = FALSE]
    # L <- eULU$val[eR]
    #  rownames(U) <- rownames(ULUPM) <- colnames(ULUPM) <- rownames(Y)
    EZ <- 0.5 * (EZ + t(EZ))
    YPM <- 0.5 * (YPM + t(YPM))
    fit <- list(BETA = BETA, VC = VC, APM = APM, U = U_vect, V = V, 
                ULUPM = ULUPM, EZ = EZ, YPM = YPM, GOF = GOF)
  }
  class(fit) <- "ame"
  fit
}



########################################### Censored Standard Binary AMEN
# this is an updated version of stanndard ammen that uses a censoring parameter
# as such, the model should be specified as 'bin', other arguments are same as AMEN
ame_censored_c_RE <- function (Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL, rvar = !(model == 
                                                                                     "rrl"), cvar = TRUE, dcor = !symmetric, nvar = TRUE, R = 0, 
                               model = "bin", intercept = !is.element(model, c("rrl", "ord")), 
                               symmetric = FALSE, odmax = rep(max(apply(Y > 0, 1, sum, 
                                                                        na.rm = TRUE)), nrow(Y)), seed = 1, nscan = 10000, burn = 500, 
                               odens = 25, plot = TRUE, print = TRUE, gof = TRUE) 
{
  set.seed(seed)
  diag(Y) <- NA
  Y_binary <- Y
  diag(Y_binary) <- 0
  if (is.element(model, c("bin", "cbin"))) {
    Y <- 1 * (Y > 0)
  }
  if (is.element(model, c("cbin", "frn", "rrl"))) {
    odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
    if (length(odmax) == 1) {
      odmax <- rep(odmax, nrow(Y))
    }
  }
  if (symmetric) {
    Xcol <- Xrow
    rvar <- cvar <- nvar
  }
  n <- nrow(Y)
  pr <- length(Xrow)/n
  pc <- length(Xcol)/n
  pd <- length(Xdyad)/n^2
  X <- design_array(Xrow, Xcol, Xdyad, intercept, nrow(Y))
  if (model == "rrl" & any(apply(apply(X, c(1, 3), var), 2, 
                                 sum) == 0) & !any(apply(X, 3, function(x) {
                                   var(c(x))
                                 }) == 0)) {
    cat("WARNING: row effects are not estimable using this procedure ", 
        "\n")
  }
  if (is.element(model, c("ord", "rrl")) & any(apply(X, 3, 
                                                     function(x) {
                                                       var(c(x))
                                                     }) == 0)) {
    cat("WARNING: an intercept is not estimable using this procedure ", 
        "\n")
  }
  if (is.element(model, c("frn", "rrl"))) {
    ymx <- max(apply(1 * (Y > 0), 1, sum, na.rm = TRUE))
    YL <- NULL
    warn <- FALSE
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(yi[rnkd]) > length(unique(yi[rnkd]))) {
        warn <- TRUE
      }
      yi[rnkd] <- rank(yi[rnkd], ties.method = "random")
      Y[i, ] <- yi
      YL <- rbind(YL, match(1:ymx, yi))
    }
    if (warn) {
      cat("WARNING: Random reordering used to break ties in ranks\n")
    }
  }
  if (model == "nrm") {
    Z <- Y
  }
  if (model == "ord") {
    Z <- matrix(zscores(Y), nrow(Y), ncol(Y))
  }
  if (model == "rrl") {
    Z <- matrix(t(apply(Y, 1, zscores)), nrow(Y), ncol(Y))
  }
  if (model == "bin") {
    Z <- matrix(zscores(Y), nrow(Y), nrow(Y))
    z01 <- 0.5 * (max(Z[Y == 0], na.rm = TRUE) + min(Z[Y == 
                                                         1], na.rm = TRUE))
    Z <- Z - z01
  }
  if (is.element(model, c("cbin", "frn"))) {
    Z <- Y
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      zi <- zscores(yi)
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
        zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 0.001
      }
      if (length(rnkd) < odmax[i]) {
        urnkd <- which(!is.na(yi) & yi == 0)
        if (max(zi[urnkd]) > 0) {
          zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 
            0.001
        }
      }
      Z[i, ] <- zi
    }
  }
  mu <- mean(Z, na.rm = TRUE)
  a <- rowMeans(Z, na.rm = TRUE)
  b <- colMeans(Z, na.rm = TRUE)
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  ZA <- mu + outer(a, b, "+")
  Z[is.na(Z)] <- ZA[is.na(Z)]
  beta <- rep(0, dim(X)[3])
  s2 <- 1
  rho <- 0
  Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar))
  U <- V <- matrix(0, nrow(Y), R)
  BETA <- matrix(nrow = 0, ncol = dim(X)[3] - pr * symmetric)
  VC <- matrix(nrow = 0, ncol = 5 - 3 * symmetric)
  UVPS <- U %*% t(V) * 0
  APS <- BPS <- rep(0, nrow(Y))
  YPS <- matrix(0, nrow(Y), ncol(Y))
  dimnames(YPS) <- dimnames(Y)
  GOF <- matrix(gofstats(Y), 1, 4)
  rownames(GOF) <- "obs"
  colnames(GOF) <- c("sd.rowmean", "sd.colmean", "dyad.dep", 
                     "triad.dep")
  names(APS) <- names(BPS) <- rownames(U) <- rownames(V) <- rownames(Y)
  if (!symmetric) {
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
    colnames(BETA) <- dimnames(X)[[3]]
  }
  if (symmetric) {
    colnames(VC) <- c("va", "ve")
    rb <- intercept + seq(1, pr, length = pr)
    cb <- intercept + pr + seq(1, pr, length = pr)
    bnames <- dimnames(X)[[3]]
    bni <- bnames[1 * intercept]
    bnn <- gsub("row", bnames[rb], replacement = "node")
    bnd <- bnames[-c(1 * intercept, rb, cb)]
    colnames(BETA) <- c(bni, bnn, bnd)
  }
  Xr <- apply(X, c(1, 3), sum)
  Xc <- apply(X, c(2, 3), sum)
  mX <- apply(X, 3, c)
  mXt <- apply(aperm(X, c(2, 1, 3)), 3, c)
  XX <- t(mX) %*% mX
  XXt <- t(mX) %*% mXt
  have_coda <- suppressWarnings(try(requireNamespace("coda", 
                                                     quietly = TRUE), silent = TRUE))
  
  censored_IND <- which(rowSums(Y_binary) == odmax)
  h <- rep(0, n)
  for (s in 1:(nscan + burn)) {
    EZ <- h + Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)

    if (model == "bin") {
      Z <- rZ_bin_fc(Z, EZ, rho, Y)
    }
   
    tmp <- rbeta_ab_fc(Z - U %*% t(V)-h, Sab, rho, X, mX, 
                       mXt, XX, XXt, Xr, Xc, s2)
    
    
    
    beta <- tmp$beta
    a <- tmp$a * rvar
    b <- tmp$b * cvar
    if (symmetric) {
      a <- b <- (a + b)/2
    }
    if (rvar & cvar & !symmetric) {
      if (is.element(model, c("nrm", "ord"))) {
        Sab <- solve(rwish(solve(diag(2) + crossprod(cbind(a, 
                                                           b))), 3 + nrow(Z)))
      }
      if (model == "bin") {
        tmp <- raSab_bin_fc(Z, Y, a, b, Sab)
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
      if (model == "cbin") {
        tmp <- raSab_cbin_fc(Z, Y, a, b, Sab, odmax, 
                             odobs)
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
      if (model == "frn") {
        tmp <- raSab_frn_fc(Z, Y, YL, a, b, Sab, odmax, 
                            odobs)
        Z <- tmp$Z
        Sab <- tmp$Sab
        a <- tmp$a
      }
    }
    if (rvar & !cvar & !symmetric) {
      Sab[1, 1] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(a^2))/2)
    }
    if (!rvar & cvar & !symmetric) {
      Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, (1 + sum(b^2))/2)
    }
    if (symmetric & nvar) {
      Sab[1, 1] <- Sab[2, 2] <- 1/rgamma(1, (1 + nrow(Y))/2, 
                                         (1 + sum(a^2))/2)
      Sab[1, 2] <- Sab[2, 1] <- 0.999 * Sab[1, 1]
    }
    
    XB <- Xbeta(X, beta)
    Se   <-matrix(c(1,rho,rho,1),2,2)*1
    iSe2 <-mhalf(solve(Se))
    td   <-iSe2[1,1] ; to<-iSe2[1,2]
    
    sigma_h <- draw_covar_c(100,100, 0, h, n)
  
    h <-  rC(Z, U%*%t(V), a, b, Xbeta(X, beta), td,to,h, censored_IND, rep(0,n), 
             sigma_h)
    
    
    if (dcor) {
      rho <- rrho_mh(Z - (Xbeta(X, beta) + h + outer(a, b, 
                                                     "+") + U %*% t(V)), rho, s2)
    }
    if (symmetric) {
      rho <- min(0.9999, 1 - 1/sqrt(s))
    }
    if (R > 0) {
      E <- Z - (Xbeta(X, beta) +h +  outer(a, b, "+"))
      if (symmetric) {
        E <- 0.5 * (E + t(E))
      }
      shrink <- (s > 0.5 * burn)
      if (symmetric) {
        UV <- rUV_sym_fc(E, U, V, s2, shrink)
      }
      if (!symmetric) {
        UV <- rUV_fc(E, U, V, rho, s2, shrink)
      }
      U <- UV$U
      V <- UV$V
    }
    if (s%%odens == 0 & s <= burn & print) {
      cat(round(100 * s/burn, 2), " pct burnin complete \n")
    }
    if (s%%odens == 0 & s > burn) {
      if (symmetric) {
        br <- beta[rb]
        bc <- beta[cb]
        bn <- (br + bc)/2
        sbeta <- c(beta[1 * intercept], bn, beta[-c(1 * 
                                                      intercept, rb, cb)])
        BETA <- rbind(BETA, sbeta)
        VC <- rbind(VC, c(Sab[1, 1], s2))
      }
      if (!symmetric) {
        BETA <- rbind(BETA, beta)
        VC <- rbind(VC, c(Sab[upper.tri(Sab, diag = T)], 
                          rho, s2))
      }
      UVPS <- UVPS + U %*% t(V)
      APS <- APS + a
      BPS <- BPS + b
      EZ <- Xbeta(X, beta) +h +  outer(a, b, "+") + U %*% 
        t(V)
      if (symmetric) {
        EZ <- (EZ + t(EZ))/2
      }
      if (model == "bin") {
        Ys <- simY_bin(EZ, rho)
      }
      if (model == "cbin") {
        Ys <- 1 * (simY_frn(EZ, rho, odmax, YO = Y) > 
                     0)
      }
      if (model == "frn") {
        Ys <- simY_frn(EZ, rho, odmax, YO = Y)
      }
      if (model == "rrl") {
        Ys <- simY_rrl(EZ, rho, odobs, YO = Y)
      }
      if (model == "nrm") {
        Ys <- simY_nrm(EZ, rho, s2)
      }
      if (model == "ord") {
        Ys <- simY_ord(EZ, rho, Y)
      }
      if (symmetric) {
        Ys[lower.tri(Ys)] <- 0
        Ys <- Ys + t(Ys)
      }
      YPS <- YPS + Ys
      if (gof) {
        Ys[is.na(Y)] <- NA
        GOF <- rbind(GOF, gofstats(Ys))
      }
      if (print) {
        cat(s, round(apply(BETA, 2, mean), 2), ":", 
            round(apply(VC, 2, mean), 2), "\n")
        if (have_coda & nrow(VC) > 3 & length(beta) > 
            0) {
          cat(round(coda::effectiveSize(BETA)), "\n")
        }
      }
      if (plot) {
        par(mfrow = c(1 + 2 * gof, 2), mar = c(3, 3, 
                                               1, 1), mgp = c(1.75, 0.75, 0))
        mVC <- apply(VC, 2, median)
        matplot(VC, type = "l", lty = 1)
        abline(h = mVC, col = 1:length(mVC))
        if (length(beta) > 0) {
          mBETA <- apply(BETA, 2, median)
          matplot(BETA, type = "l", lty = 1, col = 1:length(mBETA))
          abline(h = mBETA, col = 1:length(mBETA))
          abline(h = 0, col = "gray")
        }
        if (gof) {
          for (k in 1:4) {
            hist(GOF[-1, k], xlim = range(GOF[, k]), 
                 main = "", prob = TRUE, xlab = colnames(GOF)[k], 
                 col = "lightblue", ylab = "", yaxt = "n")
            abline(v = GOF[1, k], col = "red")
          }
        }
      }
    }
  }
  APM <- APS/nrow(VC)
  BPM <- BPS/nrow(VC)
  UVPM <- UVPS/nrow(VC)
  YPM <- YPS/nrow(VC)
  EZ <- Xbeta(X, apply(BETA, 2, mean)) + outer(APM, BPM, "+") + 
    UVPM
  names(APM) <- names(BPM) <- rownames(UVPM) <- colnames(UVPM) <- dimnames(Y)[[1]]
  dimnames(YPM) <- dimnames(EZ) <- dimnames(Y)
  rownames(BETA) <- NULL
  if (!symmetric) {
    UDV <- svd(UVPM)
    U <- UDV$u[, seq(1, R, length = R)] %*% diag(sqrt(UDV$d)[seq(1, 
                                                                 R, length = R)], nrow = R)
    V <- UDV$v[, seq(1, R, length = R)] %*% diag(sqrt(UDV$d)[seq(1, 
                                                                 R, length = R)], nrow = R)
    rownames(U) <- rownames(V) <- rownames(Y)
    fit <- list(BETA = BETA, VC = VC, APM = APM, BPM = BPM, 
                U = U, V = V, UVPM = UVPM, EZ = EZ, YPM = YPM, GOF = GOF)
  }
  if (symmetric) {
    ULUPM <- UVPM
    eULU <- eigen(ULUPM)
    eR <- which(rank(-abs(eULU$val), ties.method = "first") <= 
                  R)
    U <- eULU$vec[, seq(1, R, length = R), drop = FALSE]
    L <- eULU$val[eR]
    rownames(U) <- rownames(ULUPM) <- colnames(ULUPM) <- rownames(Y)
    EZ <- 0.5 * (EZ + t(EZ))
    YPM <- 0.5 * (YPM + t(YPM))
    fit <- list(BETA = BETA, VC = VC, APM = APM, U = U, 
                L = L, ULUPM = ULUPM, EZ = EZ, YPM = YPM, GOF = GOF)
  }
  class(fit) <- "ame"
  fit
}


############################################ COMM DEP MODEL


##### Master function for comm. dep model

# Main function to run MCMC, all covariates are assumed to be dep. on community
#' @param X_r: n x n x p_dep array of row covar
#' @param X_c: n x n x p_dep array of col covar
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param Y: Observed network (can be binary or censored binary)
#' @parm iter: number of iterations to run mcmc for 
#' @param numGroup: number of communities (also denote k)
#' @param prior_beta_mu: k dimensional vector with prior mean for beta
#' @param prior_beta_var:k x k matrix for prior variance for beta
#' @param prior_alpha_mu: constant for alpha(intercept) mean
#' @param prior_alpha_var:constant for alpha(intercept) var
#' @param start_beta_c: p_dep x k matrix for initializing beta_c
#' @param start_beta_r:  p_dep x k matrix for initializing beta_c
#' @param start_beta_dyad:  p_dep_dyad x k matrix for initializing beta_c
#' @param start_alpha: constant to initialize alpha
#' @param UandV: Default is FALSE, indicates if U and V should be calculated or just U
#' @param symmetricLambda: Default is TRUE, indicates if lambda should be symmetric
#' @param keep.UV: Default is TRUE, store U/V?
#' @param dcor: Default is TRUE, indicates if correlation in error terms should be estimated
#' @param symmetric: Default is FALSE, indicates that Y is not symmetric
#' @param model:  "cbin" or "bin", indicates what type of network Y is 
#' @param odmax: if model = "cbin", this is the max out degree a node can have
#' @param indepBeta: which covariates are indep of community (does not matter if selectDependentBeta = all)
#' @param selectDependentBeta: "all" "some" all means that all covariates will be comm. dep.
#' @param odens = 10: how often to save output
#' @param burnin: amount of burnin iterations
#' @param badInit: start with completely random initialization for U/V
#' @param prior_group_probs: For future use (if we want to change prior prob on group membership)

# we only allow for cbin and bin: for binary, you can select alll or only some beta
# to be community dependent, for censored it is ony all dep. for now
amen_master <-function(X_r, X_c, X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                       prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r, start_beta_dyad,
                       start_alpha, UandV = TRUE, symmetricLambda = TRUE, keep.UV = TRUE,
                       dcor = TRUE, symmetric = FALSE, 
                       model = "bin",odmax = NULL,
                       indepBeta=1,  selectDependentBeta = "all",
                       odens = 10, burnin = 1000, badInit =FALSE, prior_group_probs = NULL){
  
  if(model == "cbin"){
    results <- subAmen_allDepBeta_censored(X_r, X_c,X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                           prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r,start_beta_dyad,
                                           start_alpha,UandV ,symmetricLambda, keep.UV,
                                           dcor, symmetric, model,odmax, odens, burnin, badInit, prior_group_probs)     }
  
  else if(selectDependentBeta == "all"){
    
    results <- subAmen_allDepBeta(X_r, X_c,X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                  prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r,start_beta_dyad,
                                  start_alpha,UandV ,symmetricLambda, keep.UV,
                                  dcor, symmetric, model,odmax, odens, burnin, badInit, prior_group_probs)
  }else{
    
    results <- subAmen_Dep_And_IndepBeta(X_r, X_c, X_d, Y,iter,indepBeta, numGroup, prior_beta_mu, prior_beta_var,
                                         prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r, 
                                         start_beta_dyad, start_alpha,UandV, symmetricLambda, keep.UV,
                                         dcor, symmetric, model, odmax,odens, burnin)
  }
  
  return(results)
}






# give algo estimated U to stay at (for simulation, would not want to use in practice)
## MASTER FUNCTION, only for binary likelihood


amen_master_give_it_initialU <-function(X_r, X_c, X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                        prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r, start_beta_dyad,
                                        start_alpha, UandV = TRUE, symmetricLambda = TRUE, keep.UV = FALSE,
                                        dcor = FALSE, symmetric = TRUE, 
                                        model = "bin",odmax = NULL,
                                        indepBeta=1,  selectDependentBeta = "all",
                                        odens = 10, burnin = 1000, INITIALU, prior_group_probs = NULL){
  library(pbivnorm)
  
  if(selectDependentBeta == "all"){
    
    if(is.null(prior_group_probs )){
      prior_group_probs <- rep(1/numGroup, numGroup)
    }
    
    results <- subAmen_allDepBeta_stay_at_initialU(X_r, X_c,X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                                   prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r,start_beta_dyad,
                                                   start_alpha,UandV ,symmetricLambda, keep.UV,
                                                   dcor, symmetric, model,odmax, odens, burnin, INITIALU, 
                                                   prior_group_probs)
    
  }else{
    
    results <- subAmen_Dep_And_IndepBeta(X_r, X_c, X_d, Y,iter,indepBeta, numGroup, prior_beta_mu, prior_beta_var,
                                         prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r, 
                                         start_beta_dyad, start_alpha,UandV, symmetricLambda, keep.UV,
                                         dcor, symmetric, model, odmax,odens, burnin)
  }
  
  return(results)
}


# give it true U to stay at
## MASTER FUNCTION, only binary likelihood
amen_master_giveTRUTH <-function(X_r, X_c, X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                 prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r, start_beta_dyad,
                                 start_alpha, UandV = TRUE, symmetricLambda = TRUE, keep.UV = FALSE,
                                 dcor = FALSE, symmetric = TRUE, 
                                 model = "bin",odmax = NULL,
                                 indepBeta=1,  selectDependentBeta = "all",
                                 odens = 10, burnin = 1000, badInit =FALSE, prior_group_probs = NULL){
  library(pbivnorm)
  
  if(selectDependentBeta == "all"){
    
    results <- subAmen_allDepBeta_giveTRUTH(X_r, X_c,X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                            prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r,start_beta_dyad,
                                            start_alpha,UandV ,symmetricLambda, keep.UV,
                                            dcor, symmetric, model,odmax, odens, burnin, badInit)
    
  }else{
    
    results <- subAmen_Dep_And_IndepBeta(X_r, X_c, X_d, Y,iter,indepBeta, numGroup, prior_beta_mu, prior_beta_var,
                                         prior_alpha_mu, prior_alpha_var, start_beta_c, start_beta_r, 
                                         start_beta_dyad, start_alpha,UandV, symmetricLambda, keep.UV,
                                         dcor, symmetric, model, odmax,odens, burnin)
  }
  
  return(results)
}


### HElPER FUNCTIONS
# convert data frame to array of matrices for use in amen_master
#' @param df: is n x p dataframe where p is number of covariates

df_to_matrix <- function(df){
  n <- dim(df)[1]
  p <- dim(df)[2]
  X_r <- array(dim = c(n, n, p))
  X_c <- array(dim = c(n, n, p))
  
  
  for(i in  1:p){
    X_r[,,i] <- matrix(rep(df[,i], n), nrow = n, byrow = FALSE)
    X_c[,,i] <- t(matrix(rep(df[,i], n), nrow = n, byrow = FALSE))
    
  }
  
  return(list(X_r = X_r, X_c = X_c))
}




# Function to draw censoring parameter for censored individuals
# Function to draw alpha
#' @param Z: current Z value
#' @param ULV: current ULV^T
#' @param alpha: current alpha
#' @param a/b: current a/b random row col effects
#' @param XB: current XB value
#' @param td/to: current decorr values
#' @param prior_alpha: prior mean for alpha
#' @param prior_alpha_var: prior var for alpha
rC <- function(Z, ULV, a, b, XB, td,to,c, censored_IND, prior_c = rep(0,n), 
               sigma_c){
  
  n<- dim(Z)[1]
  
  X_cut_bin <- rep(0,n)
  X_cut_bin[censored_IND] <- 1
  
  
  Z <- td*(Z-XB - outer(a, b, "+")-ULV) + to*(t(Z-XB - outer(a , b, "+")-ULV))
  
  ones <- X_cut_bin
  
  
  star <- td*(ones %x% diag(n)) + to*(diag(n)%x% ones)
  
  
  m_c <-  t(c(t(Z))%*%star + prior_c %*% solve(sigma_c*diag(n)))
  
  V_c <-  solve(t(star)%*%star + solve(sigma_c*diag(n)))
  
 
  lb <- rep(-Inf,n); ub <- rep(0,n)
  b_c <-  rtmvnorm(1, V_c%*%m_c, V_c, lb, ub)
  
  c <- b_c*X_cut_bin
  
  return(c)
  
  
}

# draw variance for c
draw_covar_c <- function(alpha, beta, mu, c, n){
  sigma_c <- rinvgamma(1, alpha + n/2, beta + .5*(sum((c-mu)^2)))
  return(sigma_c)
}

# Function to draw random effects a/b
#' @param Z: current Z value
#' @param XB: added xbeta with alpha
#' @param Sab: current Sab value
#' @param rho: current rho
#' @param s2: variance (set to 1)
#' @param ULV: ULambdaV^T (current value)
r_ab_DECORR_DEP <-function(Z,XB, Sab,rho, s2=1,ULV)
  
{
  Z<-Z-ULV-XB
  
  
  n<-nrow(Z)
  ###
  ### decorrelation
  Se<-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2<-mhalf(solve(Se))
  td<-iSe2[1,1] ; to<-iSe2[1,2]
  Sabs<-iSe2%*%Sab%*%iSe2
  tmp<-eigen(Sabs)
  k<-sum(zapsmall(tmp$val)>0 )
  
  
  Zs<-td*Z+to*t(Z)
  zr<-rowSums(Zs) ; zc<-colSums(Zs) ; zs<-sum(zc) 
  ###
  
  
  ### row and column reduction
  ab<-matrix(0,nrow(Z),2)
  if(k>0)
  {
    G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
    K<-matrix(c(0,1,1,0),2,2)
    A<-n*t(G)%*%G + diag(k)
    B<-t(G)%*%K%*%G
    iA0<-solve(A)
    C0<- -solve(A+ n*B)%*%B%*%iA0
    
    iA<-G%*%iA0%*%t(G)
    C<-G%*%C0%*%t(G)
    
    
  }
  
  ### simulate a, b 
  if(k>0) 
  {
    E<- Zs
    er<-rowSums(E) ; ec<-colSums(E) ; es<-sum(ec) 
    m<-t(t(crossprod(rbind(er,ec),t(iA0%*%t(G)))) + rowSums(es*C0%*%t(G)) )
    hiA0<-mhalf(iA0)
    e<-matrix(rnorm(n*k),n,k) 
    w<-m+ t( t(e%*%hiA0) - c(((hiA0-mhalf(iA0+n*C0))/n)%*% colSums(e) ) )
    ab<- w%*%t(G)%*%solve(iSe2) 
  }
  
  list(a=ab[,1],b=ab[,2], to = to, td = td)  
  
}


# Function to draw alpha/intercept term
#' @param Z: current Z value
#' @param ULV: current ULV^T
#' @param alpha: current alpha
#' @param a/b: current a/b random row col effects
#' @param XB: current XB value
#' @param td/to: current decorr values
#' @param prior_alpha: prior mean for alpha
#' @param prior_alpha_var: prior var for alpha
rALPHA <- function(Z, ULV, alpha, a, b, XB, td,to,prior_alpha = 0, 
                   prior_alpha_var = 1){
  
  n<- dim(Z)[1]
  
  X_alpha <- matrix(1,n,n)
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  Z <- td*(Z-XB+alpha - outer(a, b, "+")-ULV) + to*(t(Z-XB+alpha - outer(a, b, "+")-ULV))
  
  X_alpha <- td*X_alpha + to*t(X_alpha)
  
  
  
  
  
  
  m_alpha = t(t(c(Z))%*%c(X_alpha) + prior_alpha %*% solve(prior_alpha_var))
  
  V_alpha = 1/(c(X_alpha) %*% c(X_alpha)+ 1/prior_alpha_var)
  
  b_alpha <-  rnorm(1, V_alpha%*%m_alpha, V_alpha) ## .94, 9.28, 98
  
  alpha <- b_alpha
  
  
  
  
  
  
  
  return(list(alpha = alpha))
  
  
}



# FROM 'amen' package
#' Simulate a and Sab from full conditional distributions under bin likelihood
#' 
#' Simulate a and Sab from full conditional distributions under bin likelihood
#' 
#' 
#' @usage raSab_bin_fc(Z, Y, a, b, Sab, Sab0=NULL, eta0=NULL, SS = round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square binary relational matrix
#' @param a current value of row effects
#' @param b current value of column effects
#' @param Sab current value of Cov(a,b)
#' @param Sab0 prior (inverse) scale matrix for the prior distribution
#' @param eta0 prior degrees of freedom for the prior distribution
#' @param SS number of iterations
#' @return \item{Z}{new value of Z} \item{Sab}{new value of Sab} \item{a}{new
#' value of a}
#' @author Peter Hoff
#' @export raSab_bin_fc
raSab_bin_fc_new <-  function(Z,Y,a,b,Sab,Sab0=NULL,eta0=NULL,SS=round(sqrt(nrow(Z))))
{
  diag(Y)<-NA
  if(is.null(Sab0)){ Sab0<-diag(2) }
  if(is.null(eta0)){ eta0<-4 }
  
  E<-Z-a%*%t(rep(1,nrow(Z))) 
  MEL<-MEU<- -E
  MEL[!is.na(Y) & Y==0]<- -Inf
  MEU[!is.na(Y) & Y==1]<-Inf 
  MEL[is.na(Y)]<- -Inf ; MEU[is.na(Y)]<- Inf
  #  diag(MEU)<-Inf;diag(MEL)<- -Inf
  lba<-apply(MEL,1,max) 
  lba[is.na(lba)]<- -Inf
  uba<-apply(MEU,1,min) 
  uba[is.na(uba)]<- Inf
  
  for(ss in 1:SS)
  {
    ea<-b*Sab[1,2]/Sab[2,2]
    sa<-sqrt(Sab[1,1]-Sab[1,2]^2/Sab[2,2])
    norm <- qnorm(runif(nrow(Z),pnorm((lba-ea)/sa),pnorm((uba-ea)/sa)))
    a<-ea+sa*norm
    # a[is.na(a)] <- -.5
    Sab<-solve(amen::rwish(solve(eta0*Sab0+crossprod(cbind(a,b))),eta0+nrow(Z)))
  }
  list(Z=E+a%*%t(rep(1,nrow(Z))),a=a,Sab=Sab)
}



# from 'ame' package version 1.3 # 
#' Simulate a and Sab from full conditional distributions under the cbin
#' likelihood 
#' Simulate a and Sab from full conditional distributions under the cbin
#' @usage raSab_cbin_fc(Z, Y, a, b, Sab, odmax, odobs, Sab0=NULL, eta0=NULL,SS =
#' round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square matrix of ranked nomination data
#' @param a current value of row effects
#' @param b current value of column effects
#' @param Sab current value of Cov(a,b)
#' @param odmax a scalar or vector giving the maximum number of nominations fo
#' each individual
#' @param odobs observed outdegree
#' @param Sab0 prior (inverse) scale matrix for the prior distribution
#' @param eta0 prior degrees of freedom for the prior distribution
#' @param SS number of iteration
#' @return \item{Z}{new value of Z} \item{Sab}{new value of Sab} \item{a}{new
#' value of a}
#' @author Peter Hoff
#' @export raSab_cbin_fc

raSab_cbin_fc<-
  
  function (Z, Y, a, b, Sab, odmax, odobs, Sab0=NULL, eta0=NULL,SS = round(sqrt(nrow(Z))))
    
  {
    
    if(is.null(Sab0)){ Sab0<-diag(2) }
    
    if(is.null(eta0)){ eta0<-4 }
    
    
    
    
    E <- Z - a %*% t(rep(1, nrow(Z)))
    
    MEL <- MEU <- -E
    
    MEL[!is.na(Y) & Y == 0] <- -Inf
    
    MEU[!is.na(Y) & Y == 1] <- Inf
    
    MEL[is.na(Y)] <- -Inf
    
    MEU[is.na(Y)] <- Inf
    
    lba <- apply(MEL, 1, max)
    
    lba[is.na(lba)] <- -Inf
    
    uba <- apply(MEU, 1, min)
    
    uba[is.na(uba)] <- Inf
    
    uba[odobs == odmax] <- Inf
    
    for (ss in 1:SS) {
      
      ea <- b * Sab[1, 2]/Sab[2, 2]
      
      sa <- sqrt(Sab[1, 1] - Sab[1, 2]^2/Sab[2, 2])
      
      a <- ea + sa * qnorm(runif(nrow(Z), pnorm((lba - ea)/sa),
                                 
                                 pnorm((uba - ea)/sa)))
      
      Sab <- solve(rwish(solve(eta0*Sab0 + crossprod(cbind(a,
                                                           
                                                           b))), eta0 + nrow(Z)))
      
    }
    
    list(Z = E + a %*% t(rep(1, nrow(Z))), a = a, Sab = Sab)
    
  }


# Gibbs step for drawing beta when we have indep betas specified
#' @param Z: current Z
#' @param X_r: row covariates without indep.
#' @param X_c: col covariates without indep.
#' @param U/lambda: current U/lambda values
#' @param p: number of dep covariates
#' @param beta_rowINPUT: current row beta
#' @param beta_colINPUT: current col beta
#' @param alpha: current intercept
#' @param prior_beta_mu: prior mean for beta
#' @param prior_beta_var: prior var for beta
#' @param a/b: current random row/col effects
#' @param h_r_past/h_c_past: past h_r h_c values (can see calculation in code below)
#' @param to/td: decorrelation values
#' @param updatedU: was U updated?
#' @param indepXB: sum of indep XBeta


rbeta_DECORR_DEP_ANDINDEP <- function(Z, X_r, X_c, U,lambda, p, beta_rowINPUT,
                                      beta_colINPUT, alpha, prior_beta_mu=rep(1, dim(U)[2]), 
                                      prior_beta_var = diag(3,dim(U)[2]), a , b , h_r_past, h_c_past, to ,td,updatedU,
                                      indepXB){
  #to = .99; td = 1
  ulu <- U %*% lambda %*% t(U)
  Z <- Z -ulu -alpha  - outer(a,b, "+")-indepXB
  n = dim(X_r)[1]
  ## prior
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  ## get all xrow+xcol
  
  beta_row <- beta_col <- matrix(0, nrow = p, ncol = dim(U)[2] ) # + alpha  + a+ b
  
  for(i in 1:p){
    if(p ==1){
      Ztemp<- Z
    }
    else if(p ==2){
      XB_ind_get    <- XBETA_DEP_partition(X_r, X_c, beta_rowINPUT, beta_colINPUT, U,dim(X_r)[3], updatedU=FALSE, h_r_past, h_c_past)
      XBETA_ROW     <- XB_ind_get$XBETAROW
      XBETA_COL     <- XB_ind_get$XBETACOL
      
      
      Ztemp <- Z - XBETA_COL[,,-i] - XBETA_ROW[,,-i] 
      
    }else{
      Ztemp <- Z - rowSums(XBETA_COL[,,-i], dims = 2) - rowSums(XBETA_ROW[,,-i], dims = 2)
      
    }
    
    Z_d <-td*Ztemp + to*t(Ztemp) 
    
    if(p ==1){
      b_c   <- beta_colINPUT
      b_r   <- beta_rowINPUT
      X_row <- X_r
      X_col <- X_c
    }else{
      b_c   <- beta_colINPUT[i,]
      b_r   <- beta_rowINPUT[i,]
      X_row <- X_r[,,i]
      X_col <- X_c[,,i]
    }
    
    
    
    
    
    id     <- diag(1,n)
    
    h_r <-  td*(ones %x% ((diag(1,n) *X_row)%*%U)) +to*((diag(1,n)*X_row)%*% U) %x% ones
    h_c <-  to*(ones %x% ((diag(1,n) *X_col)%*%U)) +td*((diag(1,n)*X_col)%*% U) %x% ones
    h_rT_h_r <- t(h_r)%*% h_r
    h_cT_h_c <- t(h_c) %*% h_c
    
    c = h_c %*% c(b_c)
    
    
    
    m_r = t(t(c(Z_d)-c)%*%h_r + prior_beta_mu %*% solve(prior_beta_var))
    
    V_r = solve((h_rT_h_r + solve(prior_beta_var)))
    
    b_r <-  rmvnorm(1, V_r%*%m_r, V_r) ## .94, 9.28, 98
    
    if(p  == 1){
      beta_rowINPUT <- b_r
      
    }else{
      beta_rowINPUT[i, ] <- b_r
      
    }   
    print(b_r)
    ## we need to condition on rowXB to get colXB
    c2 = h_r %*% c(b_r)
    
    
    m_c = t(t(c(Z_d)-c2)%*%h_c + prior_beta_mu %*% solve(prior_beta_var))
    
    V_c = solve((h_cT_h_c + solve(prior_beta_var)))
    
    
    b_c <- rmvnorm(1, V_c%*%m_c, V_c)
    
    
    beta_row[i,] <- b_r
    beta_col[i,] <- b_c
    
    if(p  == 1){
      beta_colINPUT <- b_c
      
    }else{
      beta_colINPUT[i, ] <- b_c
      
    }
    
  }
  
  
  return(list(beta_col = beta_col, beta_row = beta_row))
  
  
}

# Gibbs step for drawing beta when we have all dependent betas specified
#' @param Z: current Z value
#' @param X_r: row covariates
#' @param X_c: col covariates
#' @param XBETA_DYAD: added dyadic covariates
#' @param U/V/lambda: current values
#' @param p: number of dep covariates
#' @param beta_rowINPUT: current beta for row
#' @param beta_colINPUT: current beta for col
#' @param alpha: current intercept values
#' @param prior_beta_mu: prior mean for beta
#' @param prior_beta_var: prior beta variance
#' @param a/b: current row/col random effects
#' @param td/to: current decorrelation values
#' @param h_r_past/h_c_past: current values
#' @param updatedUV: does U/V need to be updated?

rbeta_DECORR_DEP <-function(Z, X_r, X_c, XBETA_DYAD, U,V, lambda,
                            p, beta_rowINPUT, beta_colINPUT, alpha, prior_beta_mu=rep(1, dim(U)[2]), 
                            prior_beta_var = diag(3,dim(U)[2]), a , b , h_r_past, h_c_past, to ,td,
                            updatedUV ){
  ulv <- U %*% lambda %*% t(V)
  Z <- Z -ulv -alpha  - outer(a,b, "+") - XBETA_DYAD
  n = dim(X_r)[1]
  ## prior
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  ## get all xrow+xcol
  
  beta_row <- beta_col <- matrix(0, nrow = p, ncol = dim(U)[2] ) # + alpha  + a+ b
  
  for(i in 1:p){
    XB_ind_get    <- XBETA_DEP_UPDATEROWCOLUMN(X_r, X_c, beta_rowINPUT, beta_colINPUT, 
                                               U,V,
                                               updatedUV= FALSE, 
                                               h_r_past, 
                                               h_c_past)
    XBETA_ROW     <- XB_ind_get$XBETAROW
    XBETA_COL     <- XB_ind_get$XBETACOL
    if(p ==1){
      
      Ztemp <- Z 
      
    }else if(p ==2){
      
      
      Ztemp <- Z - XBETA_COL[,,-i] - XBETA_ROW[,,-i]
      
    }else{
      Ztemp <- Z - rowSums(XBETA_COL[,,-i], dims = 2) - rowSums(XBETA_ROW[,,-i], dims = 2)
      
    }
    
    Z_d <-td*Ztemp + to*t(Ztemp) 
    
    b_c   <- beta_colINPUT[i,]
    b_r   <- beta_rowINPUT[i,]
    
    if(p==1){
      X_row <- X_r
      X_col <- X_c
      
    }else{
      X_row <- X_r[,,i]
      X_col <- X_c[,,i]
      
    }
    
    
    
    id     <- diag(1,n)
    
    h_r <-  td*(ones %x% ((diag(1,n) *X_row)%*%U)) +to*((diag(1,n)*X_row)%*% U) %x% ones
    h_c <-  to*(ones %x% ((diag(1,n) *X_col)%*%V)) +td*((diag(1,n)*X_col)%*% V) %x% ones
    h_rT_h_r <- t(h_r)%*% h_r
    h_cT_h_c <- t(h_c) %*% h_c
    
    c = h_c %*% c(b_c)
    
    
    
    m_r = t(t(c(Z_d)-c)%*%h_r + prior_beta_mu %*% solve(prior_beta_var))
    
    V_r = solve((h_rT_h_r + solve(prior_beta_var)))
    
    b_r <-  rmvnorm(1, V_r%*%m_r, V_r) ## .94, 9.28, 98
    beta_rowINPUT[i, ] <- b_r
    print(b_r)
    ## we need to condition on rowXB to get colXB
    c2 = h_r %*% c(b_r)
    
    
    m_c = t(t(c(Z_d)-c2)%*%h_c + prior_beta_mu %*% solve(prior_beta_var))
    
    V_c = solve((h_cT_h_c + solve(prior_beta_var)))
    
    
    b_c <- rmvnorm(1, V_c%*%m_c, V_c)
    #  print(b_c)
    #   print(b_c)
    
    beta_row[i,] <- b_r
    beta_col[i,] <- b_c
    beta_colINPUT[i, ] <- b_c
    
  }
  
  
  return(list(beta_col = beta_col, beta_row = beta_row))
  
  
}

# Function to draw beta for dyadic covariates
# Main function to run MCMC, all covariates are assumed to be dep. on community
#' @param Z: current Z
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param U/V/ lambda:current values for U, V, lambda
#' @param beta_rowINPUTD: current dep beta for row dyad
#' @param beta_colINPUTD: current dep beta for row dyad
#' @param prior_beta_mu: k dimensional vector with prior mean for beta
#' @param prior_beta_var:k x k matrix for prior variance for beta
#' @param a/b: current a/b random effect values
#' @param to/td: current decorrelation 
#' @param numEigs: number of eigen values/vectors to consider for dyadic decomp.
#' @param  h_r_past/h_c_past: past values
rbeta_DECORR_DEPDYAD <-function(Z, X_d, U,V,lambda, beta_rowINPUTD, beta_colINPUTD,
                                prior_beta_mu=rep(1, dim(U)[2]), 
                                prior_beta_var = diag(3,dim(U)[2]), a, b,  to ,td,
                                numEigs,
                                h_r_past, h_c_past){
  
  n = dim(X_d)[1]
  p = dim(X_d)[3]
  Z <- Z - U%*%lambda%*%t(V) - outer(a, b, "+")
  beta_row_dyad <- beta_col_dyad <- matrix(0, nrow = p, ncol = dim(U)[2] ) # + alpha  + a+ b
  
  for(i in 1:p){
    print("I")
    print(i)
    
    XB_ind_get    <- XBETA_DEP_UPDATEDDYAD(X_d, beta_rowINPUTD, beta_colINPUTD,
                                           U,V, 
                                           FALSE, 
                                           h_r_past, 
                                           h_c_past)
    
    XBETA_DYAD    <- XB_ind_get$XBETADYAD
    
    if(p ==2){
      
      Ztemp <- Z - XBETA_DYAD[,,-i] 
      
    }else{
      Ztemp <- Z - rowSums(XBETA_DYAD[,,-i], dims = 2) 
      
    }
    
    Z_d <-td*Ztemp + to*t(Ztemp) 
    
    b_dr   <- beta_rowINPUTD[i,]
    
    X_dyad <- X_d[,,i]
    
    
    eigX <- svd(X_dyad)
    valsE <- eigX$d
    Uv    <- eigX$u
    Vv    <- eigX$v
    
    
    
    
    eigXT <- svd(t(X_dyad))
    valsET <- eigXT$d
    UvT    <- eigXT$u
    VvT    <- eigXT$v
    
    h_dc <- h_dr <- matrix(0, n^2, dim(beta_rowINPUTD)[2])
    
    
    
    
    for(j in 1:numEigs){
      currH_dcC <- (diag(c(Vv[,j]))%*%V) %x% (diag(c(U%*%b_dr))%*% (valsE[j]*Uv[,j]))
      #currH_dcD <- to*(diag(c(U%*%b_dr)) %*% Vv[,j]) %x% (diag(c(valsE[j]*Uv[,j]))%*% U)
      
      currH_dcD <- (diag(c(U%*%b_dr)) %*% VvT[,j]) %x% (diag(c(valsET[j]*UvT[,j]))%*% V)
      
      h_dc <- td*currH_dcC + to*currH_dcD + h_dc
      
    }
    print("HDC")
    print(summary(h_dc))
    h_cT_h_c <- t(h_dc)%*% h_dc
    
    m_c = t(t(c(Z_d))%*%h_dc + prior_beta_mu %*% solve(prior_beta_var))
    
    V_c = solve((h_cT_h_c + solve(prior_beta_var)))
    
    print("VAR C")
    print(V_c)
    
    
    b_dc <- rmvnorm(1, V_c%*%m_c, V_c)
    print("dyad BC")
    print(b_dc)
    b_dc <- c(b_dc)
    
    
    for(j in 1:numEigs){
      
      currH_drC <- (diag(c(V%*%b_dc)) %*% Vv[,j]) %x% (diag(c(valsE[j]*Uv[,j]))%*% U)
      currH_drD <- (diag(c(VvT[,j]))%*%U) %x% (diag(c(V%*%b_dc))%*% (valsET[j]*UvT[,j]))
      
      h_dr <- td*currH_drC + to*currH_drD + h_dr
      
    }
    
    
    print("HDR")
    print(summary(h_dr))
    
    h_rT_h_r <- t(h_dr)%*% h_dr
    
    
    
    m_r = t(t(c(Z_d))%*%h_dr + prior_beta_mu %*% solve(prior_beta_var))
    
    V_r = solve((h_rT_h_r + solve(prior_beta_var)))
    
    
    print("VAR R")
    print(V_r)
    
    b_dr <-  rmvnorm(1, V_r%*%m_r, V_r) ## .94, 9.28, 98
    b_dr <- c(b_dr)
    print("BETAR")
    print(b_dr)
    
    #   print(b_c)
    
    beta_row_dyad[i,]  <- b_dr
    beta_col_dyad[i,]  <- b_dc
    
    beta_rowINPUTD[i, ] <- b_dr
    beta_colINPUTD[i, ] <- b_dc
    
  }
  
  return(list(beta_col_dyad = beta_col_dyad, 
              beta_row_dyad = beta_row_dyad))
  
  
}






# CURRENTLY FOR INDEP. WE ONLY HAVE U (not U and V) cability and binary network
#' @param Z: current Z value
#' @param X_r: row covariates
#' @param X_c: col covariates
#' @param U/lambda: current values
#' @param p: number of indep covariates
#' @param beta_rowINPUT: current beta for row
#' @param beta_colINPUT: current beta for col
#' @param alpha: current intercept values
#' @param a/b: current row/col random effects
#' @param XB_dep: dependent XB
#' @param td/to: current decorrelation values
rbeta_INDEPDECORR <- function(Z, X_r, X_c, U,lambda, p, beta_rowINPUT, beta_colINPUT,
                              alpha, a , b, XB_dep,td,to){
  n = dim(X_r)[1]
  ## prior
  prior_beta <-0
  prior_beta_var <-1
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  ## get all xrow+xcol
  elseX<- td*(alpha  + XB_dep + outer(a, b, "+")) + to*(t(alpha  + XB_dep + outer(a, b, "+")))
  Z <- td*Z + to*t(Z)
  beta_row <- beta_col <- c()
  
  if(p == 1){
    b_c <- beta_colINPUT
    b_r <- beta_rowINPUT
    X_row <- td*X_r + to*t(X_r)
    X_col <- td*X_c   + to*t(X_c)
    
    dU <- c(td*U%*%lambda%*%t(U))+c(to*t(U%*%lambda%*%t(U)))
    
    
    c = c(X_col)* b_c+ dU + c(elseX)
    
    
    m_r = t(t(c(Z))%*%c(X_row) -t(c)%*%c(X_row) + prior_beta %*% solve(prior_beta_var))
    
    V_r = 1/(c(X_row) %*% c(X_row)+ 1/prior_beta_var)
    
    b_r <-  rnorm(1, V_r%*%m_r, V_r) ## .94, 9.28, 98
    
    print(b_r)
    ## we need to condition on rowXB to get colXB
    c = b_r*c(X_row)+ dU + c(elseX)
    
    
    m_c = t(t(c(Z))%*%c(X_col)-t(c)%*%c(X_col) + prior_beta %*% solve(prior_beta_var))
    
    V_c = 1/(c(X_col) %*% c(X_col)+ 1/prior_beta_var)
    
    
    b_c <- rnorm(1, V_c%*%m_c, V_c)
    
    beta_row <- b_r
    beta_col<- b_c
    
  }else{
    X_else <- Xbeta_indep_Individual(X_r, X_c,beta_rowINPUT, beta_colINPUT,p) 
    for(i in 1:p){
      b_c <- beta_colINPUT[i]
      b_r <- beta_rowINPUT[i]
      X_row <- td*X_r[,,i]+to*t(X_r[,,i])
      X_col <- td*X_c[,,i]+to*t(X_c[,,i])   
      
      X_elsesub <- td*(X_else[,,-i] + outer(a, b, "+") + alpha + XB_dep)+to*t((X_else[,,-i] + outer(a, b, "+") + alpha + XB_dep))
      
      if(is.na(dim(X_elsesub)[3])){
        elseX<-  X_elsesub
        
      }else{
        elseX<-  rowSums(X_elsesub, dims = 2)
        
      }
      dU <- c(td*U%*%lambda%*%t(U))+c(to*t(U%*%lambda%*%t(U)))
      
      c = c(X_col)* b_c+ dU + c(elseX)
      
      
      m_r = t(t(c(Z))%*%c(X_row) -t(c)%*%c(X_row) + prior_beta %*% solve(prior_beta_var))
      
      V_r = 1/(c(X_row) %*% c(X_row)+ 1/prior_beta_var)
      
      b_r <-  rnorm(1, V_r%*%m_r, V_r) ## .94, 9.28, 98
      
      print(b_r)
      ## we need to condition on rowXB to get colXB
      c = b_r*c(X_row)+ dU + c(elseX)
      
      
      m_c = t(t(c(Z))%*%c(X_col)-t(c)%*%c(X_col) + prior_beta %*% solve(prior_beta_var))
      
      V_c = 1/(c(X_col) %*% c(X_col)+ 1/prior_beta_var)
      
      
      b_c <- rnorm(1, V_c%*%m_c, V_c)
      
      
      beta_row[i] <- b_r
      beta_col[i] <- b_c
    }
    
  }
  
  
  
  
  
  return(list(beta_col = beta_col, beta_row = beta_row))
  
  
}


# Function using MH step to draw rho from # from 'ame' package version 1.3 # 
#' @param Z: current Z value
#' @param rho: current rho value
#' @param s2: variance set at 1
#' @param offset: Expected value of Z
#' @param asp: set to TRUE

rrho_mh_new <- function(Z,rho,s2=1,offset=0,asp=NULL)
{
  if(is.null(asp)){ asp<-TRUE } 
  
  E<-Z-offset
  EM<-cbind(E[upper.tri(E)],t(E)[upper.tri(E)] )/sqrt(s2)
  emcp<-sum(EM[,1]*EM[,2])
  emss<-sum(EM^2)
  
  m<- nrow(EM)
  sr<- 2*(1-cor(EM)[1,2]^2)/sqrt(m)
  
  rho1<-rho+sr*qnorm( runif(1,pnorm( (-1-rho)/sr),pnorm( (1-rho)/sr)))
  
  lhr<-(-.5*(m*log(1-rho1^2)+(emss-2*rho1*emcp)/(1-rho1^2)))-
    (-.5*(m*log(1-rho^2 )+(emss-2*rho*emcp )/(1-rho^2 )))   +
    asp*( (-.5*log(1-rho1^2)) - (-.5*log(1-rho^2)) )
  
  if(log(runif(1))<lhr) { rho<-rho1 }
  min(abs(rho),.995)*sign(rho)
}


# from ame package code: draws Sab
#' @param a/b: current row/col random effects
#' @param Sab0/eta0: prior parameters
rSab_fc<-function(a,b,Sab0=NULL,eta0=NULL) 
{
  
  if(is.null(Sab0)){ Sab0<-diag(2) } 
  if(is.null(eta0)){ eta0<-4 }
  
  solve(amen::rwish(solve(eta0*Sab0+crossprod(cbind(a, b))), eta0+length(a)))
}


#' Simulate Z based on a probit model # from 'ame' package version 1.3 # 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and a binary relational matrix Y
#' @usage rZ_bin_fc(Z, EZ, rho, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param Y square binary relational matrix
#' @return a square matrix , the new value of Z
#' @author Peter Hoff
#' @export rZ_bin_fc
rZ_bin_fc_new <-  function(Z,EZ,rho,Y)
{ 
  # simulates Z under the contraints
  # (1)  Y[i,j]=1   => Z[i,j]>0
  # (2)  Y[i,j]=0   => Z[i,j]<0
  diag(Y) <- NA
  sz<-sqrt(1-rho^2)
  ut<-upper.tri(EZ)
  lt<-lower.tri(EZ)
  
  Y[is.na(Y)]<- -1
  for(y in c((-1):1))
  { 
    lb<-c(-Inf,-Inf,0)[y+2] ; ub<-c(Inf,0,Inf)[y+2]
    
    for(tri in 1:2)
    { 
      if(tri==1){ up<-ut & Y==y }
      if(tri==2){ up<-lt & Y==y }
      
      ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
      zup<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
      zerr<-which(abs(zup)==Inf)
      if(length(zerr)>0){ zup[zerr]<-(Z[up])[zerr] }
      Z[up]<-zup
    }
  }
  
  ##
  c<-(sqrt(1+rho) + sqrt(1-rho))/2
  d<-(sqrt(1+rho) - sqrt(1-rho))/2
  E<-matrix(rnorm(nrow(Y)^2),nrow(Y),nrow(Y))
  ZP<-EZ + c*E + d*t(E) 
  A<-( (Y== -1) | ( sign(ZP) == sign(Y-.5)) ) ; diag(A)<-TRUE
  A<-A & t(A)
  Z[A]<-ZP[A]
  ##
  
  diag(Z)<-rnorm(nrow(Z),diag(EZ),sqrt(1+rho))
  ##
  
  Z
}


#' Simulate Z given fixed rank nomination data# from 'ame' package version 1.3 # 
#' 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and censored binary nomination data
#' 
#' simulates Z under the constraints (1) Y[i,j]=1, Y[i,k]=0 => Z[i,j]>Z[i,k] ,
#' (2) Y[i,j]=1 => Z[i,j]>0 , (3) Y[i,j]=0 & odobs[i]<odmax[i] => Z[i,j]<0
#' 
#' @usage rZ_cbin_fc(Z, EZ, rho, Y, odmax, odobs)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param Y square matrix of ranked nomination data
#' @param odmax a scalar or vector giving the maximum number of nominations for
#' each individual
#' @param odobs observed outdegree
#' @return a square matrix, the new value of Z
#' @author Peter Hoff
#' @export rZ_cbin_fc
rZ_cbin_fc_new <-function(Z,EZ,rho,Y,odmax,odobs){
  # simulates Z under the contraints 
  # (1)  Y[i,j] > Y[i,k]              => Z[i,j]>Z[i,k]  
  # (2)  Y[i,j]>0                     => Z[i,j]>0    
  # (3)  Y[i,j]=0 & odobs[i]<odmax[i] => Z[i,j]<0
  
  sz<-sqrt(1-rho^2)
  ut<-upper.tri(EZ)
  lt<-lower.tri(EZ)
  
  for(y in sample(0:1))
  {
    if(y==1)
    { 
      ub<- Inf
      lbm<-matrix(pmax(apply(Z-(Y!=0)*(Inf^(Y!=0)),1,max,na.rm=TRUE),0),
                  nrow(Z),nrow(Z))
    }
    
    if(y==0)
    { 
      lb<- -Inf
      ubm<-matrix(apply(Z+(Y!=1)*(Inf^(Y!=1)),1,min,na.rm=TRUE),nrow(Z), 
                  nrow(Z))
      ubm[ odobs<odmax ] <- 0 
    }
    
    up<- ut & Y==y 
    if(y==0)  { ub<-ubm[up] }  
    if(y==1)  { lb<-lbm[up] }
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
    
    up<- lt & Y==y
    if(y==0)  { ub<-ubm[up] }
    if(y==1)  { lb<-lbm[up] }
    ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
    Z[up]<-ez+sz*qnorm(runif(sum(up),pnorm((lb-ez)/sz),pnorm((ub-ez)/sz)))
  }
  
  diag(Z)<-rnorm(nrow(Z),diag(EZ),sqrt(1+rho))
  Z[is.na(Y)]<- rnorm(sum(is.na(Y)),EZ[is.na(Y)],1) # not right - doesn't use rho
  
  Z
}

# Function to simulate a network based off of current state # from 'ame' package version 1.3 # 
#' @param EZ: current expected Z value
#' @param rho: current rho value
simY_bin <-function (EZ, rho) 
{
  
  ZS <- simZ(EZ, rho)
  YS <- 1 * (ZS > 0)
  diag(YS) <- NA
  YS
}


# Function to simulate FRN Y # from 'ame' package version 1.3 # 
#' @param EZ: current expected Z value
#' @param rho: current rho value
#' @param odmax: max number of connections a node can have
#' @param Y0: NA values
simY_frn <-function (EZ, rho, odmax, YO = NULL) 
{
  if (length(odmax) == 1) {
    odmax <- rep(odmax, nrow(EZ))
  }
  
  ZS <- simZ(EZ, rho)
  diag(ZS) <- -Inf
  if (!is.null(YO)) {
    ZS[is.na(YO)] <- -Inf
  }
  YS <- ZS * 0
  for (i in 1:nrow(EZ)) {
    rs <- rank(ZS[i, ]) - (nrow(EZ) - odmax[i])
    YS[i, ] <- rs * (rs > 0) * (ZS[i, ] > 0)
    YS[i, YS[i, ] > 0] <- match(YS[i, YS[i, ] > 0], sort(unique(YS[i, 
                                                                   YS[i, ] > 0])))
  }
  diag(YS) <- NA
  if (!is.null(YO)) {
    YS[is.na(YO)] <- NA
  }
  YS
}


# Function to draw lambda
#' @param Z: current Z value 
#' @param XB: current XB 
#' @param U/V: current U/V value
#' @param  prior.lambda.prec: prior precision matrix for lambda (k x k)
#' @param  prior.lambda.mu: k^2 vector with prior mean for lambda
#' @param symmetricLambda: TRUE/FALSE indicates if lambda is symmetric
drawLambda <-function(Z, XB, U,V, prior.lambda.var,
                      prior.lambda.mu, td, to, symmetricLambda = TRUE){
  prodUV <-V %x% U
  
  prodVU <- U %x% V
  XB <- td*XB + to*t(XB)
  Z <- td*Z + to*t(Z)
  star_UV <- td*prodUV + to*prodVU
  
  
  Sigma <- solve( t(star_UV)%*%star_UV +  solve(prior.lambda.var))
  
  mu <- c(c(Z-XB)%*% star_UV) +  t(prior.lambda.mu%*%solve(prior.lambda.var))
  
  lambdaVect <- rmvnorm(1, as.matrix(Sigma%*%mu), Sigma)
  
  lambdaVect[is.na(lambdaVect)] <- 0
  
  lambda <- matrix(lambdaVect, nrow = dim(U)[2])
  
  makeSymm <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
  
  if(symmetricLambda ==TRUE){
    lambda <- makeSymm(lambda)
  }
  
  return(lambda)
}



# function for updating UV when we have censoringi

#' Function to draw U, U/V
#' Function to calculate dep X beta that is NOT all added together (returns array of XB)
#' @param lambda: current lambda (k x k matrix)
#' @param U: current U (n x k matrix)
#' @param V: current V (n x k matrix)
#' @param XB: current XB value
#' @param a,b: n vectors with current row/col random effects
#' @param Y: Observed network
#' @param UandV: TRUE/FALSE (are both estimated or just U)
#' @param counter: counts total iterations with burnin
#' @param MH: counts acceptance of MH step
#' @param X_r: n x n x p_dep array of row covar
#' @param X_c: n x n x p_dep array of col covar
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param beta_r: current beta_r, p_dep x k matrix
#' @param beta_c: current beta_c, p_dep x k matrix
#' @param alpha: current alpha value
#' @param beta_dr: current beta_dr, p_dep_dyad x k matrix
#' @param beta_dc: current beta_dc, p_dep_dyad x k matrix
#' @param h_r_past: Previous calculation of h_r
#' @param h_c_past: Previous calculation of h_c
#' @param Z: current Z value (n x n matrix)
#' @param td, to: current de-correlation values
#' @param Sab: covar from a,b
#' @param rho: current rho value
#' @param predY: Binary version of our network that is predicted from the Z values 
#' @param c: Estimated cutoff, non zero for censored individuals
#' @param numU: The number of nodes that were changed



# function to draw U/V values and to perform metropolis step under the censored binary likelihood
drawUV_probitLikelihood_ypred_cbin <- function (lambda, U, V,  XB, a,b,Y,c, UandV,counter, MH,
                                                X_r, X_c,X_d, beta_r, beta_c, alpha,
                                                beta_dr, beta_dc,
                                                h_r_past, h_c_past,Z,td,to,Sab,rho, prior_prop_in_groups=NULL) {
  
  
  
  z <- dim(U)[2]
  y.post <-Y
  Z2 <- Z
  u.post <- U
  v.post <- V
  
  if(is.null(prior_prop_in_groups)){
    prior_prop_in_groups <- rep(1/z, z)
  }
  
  
  ####################################################################### Metropolis
  
  
  # sampling BOTH U and V
  if(UandV == TRUE){
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    numV  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    oldV <-  apply(V[numV,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- prop_groupV <- c()
    currU       <- u.post[numU,]
    currV       <- v.post[numV,]
    curr_groupU <- curr_groupV <- c()
    currRowU <- currRowV    <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      currRowV   <- currV[k,]
      curr_groupV[k] <- which(currRowV == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1, prob =prior_prop_in_groups)
      temp_propV     <- sample(1:z,1, prob =prior_prop_in_groups)
      
      currRowU       <- currU[j,]
      currRowV       <- currV[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      curr_groupV[j] <- which(currRowV ==1)
      prop_groupV[j] <- temp_propV
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    temp_v.post                <- v.post
    
    temp_v.post[numV,]          <- rep(0,z)
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      temp_v.post[numV[h],prop_groupV[h]] <- 1
      
      
    }
    ########################################### Compare
    # recalculate new XB values for proposed U/V
    XB_ind_getP    <- XBETA_DEP_partition(X_r, X_c,X_d, beta_r, beta_c,
                                          beta_dr, beta_dc, temp_u.post,
                                          temp_v.post,
                                          TRUE, h_r_past, h_c_past)
    
    XBETA_ROWP     <- XB_ind_getP$XBETAROW
    XBETA_COLP     <- XB_ind_getP$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYADP    <- rowSums(XB_ind_getP$XBETADYAD,dims = 2)
    }else{
      XBETA_DYADP <- 0
    }
    
    XB_indP        <- XBETA_ROWP + XBETA_COLP 
    XBP            <- rowSums(XB_indP, dims =2 )+ alpha + 
      XBETA_DYADP # calculate XB
    
    
    # calculate likelihoods
    accept_prob      <-  probitLikelihood_cbin(temp_u.post,temp_v.post, lambda, 
                                               XBP+outer(a,b, "+")+c ,predY,Y,c(numU,numV),rho,td,to,c,Z)
    
    accept_prob_past <- probitLikelihood_cbin(u.post,v.post, lambda,
                                              XB+ outer(a,b,"+")+c, predY,Y,c(numU,numV),rho,td,to,c,Z)
    
    
    
    prob_groups <- prior_prop_in_groups
    prior_weight <- past_weight <- c()
    
    for(q in 1:length(prop_groupU)){
      prior_weight[q] <- prob_groups[prop_groupU[q]]
      past_weight[q] <-  prob_groups[curr_groupU[q]]
      
    }
    accept_prior <- sum(log(prior_weight))
    accept_past  <- sum(log(past_weight))
    ratio_pro <- accept_prob +accept_prior -accept_prob_past- accept_past
    runi <- runif(1)
    
    ################################################### accept reject
    if( ratio_pro > log(runi)) {
      
      ################################################### accept reject
      u.post <- temp_u.post
      v.post <- temp_v.post
      updateU <- TRUE
      swapPeep <- c(numU, numV)
      MH <- MH + 1
    }else{
      u.post <- u.post
      v.post <- v.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    counter <- counter + 1
  }else{
    # we only sample U
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- c()
    currU       <- u.post[numU,]
    curr_groupU <-  c()
    currRowU   <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      temp_propU     <- sample(1:z,1, prob = prior_prop_in_groups)
      currRowU       <- currU[j,]
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
    }
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      
      
    }
    # calculate new XB for proposal
    XB_ind_getP    <- XBETA_DEP_partition(X_r, X_c,X_d, beta_r, beta_c,
                                          beta_dr, beta_dc, temp_u.post,
                                          temp_u.post,
                                          TRUE, h_r_past, h_c_past)
    
    XBETA_ROWP     <- XB_ind_getP$XBETAROW
    XBETA_COLP     <- XB_ind_getP$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYADP    <- rowSums(XB_ind_getP$XBETADYAD,dims = 2)
    }else{
      XBETA_DYADP <- 0
    }
    
    XB_indP        <- XBETA_ROWP + XBETA_COLP 
    XBP            <- rowSums(XB_indP, dims =2 )+ alpha +XBETA_DYADP # calculate XB
    
    ########################################### Compare
    
    
    accept_prob      <-  probitLikelihood_cbin(temp_u.post,temp_u.post, lambda, 
                                               XBP+outer(a,b, "+")+c ,predY,Y,numU,rho,td,to,c,Z)
    
    accept_prob_past <- probitLikelihood_cbin(u.post,u.post, lambda,
                                              XB+ outer(a,b,"+")+c, predY,Y,numU,rho,td,to,c,Z)
    
    
    
    
    prob_groups <- prior_prop_in_groups
    prior_weight <- past_weight <- c()
    
    for(q in 1:length(prop_groupU)){
      prior_weight[q] <- prob_groups[prop_groupU[q]]
      past_weight[q] <-  prob_groups[curr_groupU[q]]
      
    }
    accept_prior <- sum(log(prior_weight))
    accept_past  <- sum(log(past_weight))
    ratio_pro <- accept_prob +accept_prior -accept_prob_past- accept_past
    runi <- runif(1)
    
    ###################################################3 accept reject
    if( ratio_pro > log(runi)) {
      u.post <- temp_u.post
      v.post <- temp_u.post
      updateU <- TRUE
      swapPeep <- c(numU)
      MH <- MH + 1
      print("RATIOS")
      
    }else{
      u.post <- u.post
      v.post <- u.post
      updateU <- FALSE
      swapPeep <- NULL 
    }
    
    
    U <- u.post
    V <- v.post
    
    
    counter <- counter + 1
  }
  print("MH")
  print(MH/counter)
  # swapPeep is who we changed, U,V are curr. U and V values, updateU indicates that U/V was updated, 
  # counter is overall iter number and MH is number of accepted runs
  return(list(U = U, V = V, updateU = updateU,
              swapPeep = swapPeep, 
              MH =MH, counter = counter))
}




## Drawing UV when we have NO community dep.
drawUV_normalLiknodep <- function (lambda, U, V,  XB, a,b,Y, UandV,counter, MH,
                                   h_r_past, h_c_past,Z,td,to,Sab,rho) {
  z <- dim(U)[2]
  
  #y.post    <- Y[lower.tri(Y, diag = FALSE)]
  
  y.post <-Y
  Z2 <- Z
  u.post <- U
  v.post <- V
  
  
  ####################################################################### Metropolis
  
  normalLikelihood <- function(U, V, lambda,XB, Z2,td,to,numU){
    Z2 <- Z2
    
    EZ2 <- (U%*%lambda %*%t(V) + XB)
    
    Zd <- c(td*Z2 + to*t(Z2))
    EZd <- c(td*EZ2 + to*t(EZ2))
    # EZdU <- c(EZd[numU,])#, EZd[,numU])
    # ZdU <- c(Zd[numU,])#, Zd[,numU])
    print(head(cbind(c(Zd), c(EZd))))
    print(tail(cbind(c(Zd), c(EZd))))
    print(numU)
    print(sum(c(Zd-EZd)))
    
    f <- dnorm(Zd-EZd,0, 1, log = TRUE)
    
    f2 <- sum(f)
    print(f2)
    
    return(f2)
  }
  
  
  probitLikelihood <- function(U, V, lambda, XB, Z2, Y,numU,rho,td,to){
    ULV <- U%*%lambda%*%t(V)
    n <- dim(Z2)[1]
    meanMat <- XB + ULV
    
    lh <-0
    # what we care about...
    
    for(i in numU){
      for(j in 1:n){
        if(i !=j){
          if(Y[i,j]==1 && Y[j,i] == 1){
            
            
            tmp3 <-  pbivnorm(meanMat[i,j], y = meanMat[j,i], rho) 
            if(tmp3 < .000000001){tmp3 <- .0000000001}
            if((1-tmp3) < abs(.00001)){
              tmp3 <- .9999999999
            }   
            lh <- lh + log(tmp3)
            #lh<-lh+log(tmp1) +log(tmp2)
          }else if(Y[i,j]==1 && Y[j,i] == 0){
            
            
            tmp3 <-  pbivnorm(meanMat[i,j], y = -meanMat[j,i], -rho) 
            if(tmp3 < .000000001){tmp3 <- .0000000001}
            if((1-tmp3) < abs(.00001)){
              tmp3 <- .9999999999
            }          
            lh <- lh + log(tmp3)
            
            # lh<-lh+log(tmp1) +log(1-tmp2)
            
          }else if(Y[i,j]==0 && Y[j,i] == 1){
            
            tmp3 <-  pbivnorm(-meanMat[i,j], y = meanMat[j,i], -rho) 
            if(tmp3 < .000000001){tmp3 <- .0000000001}
            if((1-tmp3) < abs(.00001)){
              tmp3 <- .9999999999
            }           
            
            lh <- lh + log(tmp3)
            
          }else{
            tmp3 <-  pbivnorm(-meanMat[i,j], y = -meanMat[j,i], rho) 
            if(tmp3 < .000000001){tmp3 <- .0000000001}
            if((1-tmp3) < abs(.00001)){
              tmp3 <- .9999999999
            }   
            lh <- lh + log(tmp3)
            
            # lh <- lh+log(1-tmp1) +log(1-tmp2)
            
          }
          # print(tmp3);print(i);print(j); print(Y[i,j])
        }
      }
      
    }
    
    tmp1 <- pnorm(diag(meanMat)[numU], sd =sqrt(1 + rho))
    
    tmp1 <- ifelse(tmp1 == 0, .000000000001, tmp1)
    tmp1 <- ifelse(tmp1 == 1, .999999999999, tmp1)
    lh <-lh+ sum(log(1-tmp1))
    
    return(lh)
  }
  
  if(UandV == TRUE){
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    numV  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    oldV <-  apply(V[numV,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- prop_groupV <- c()
    currU       <- u.post[numU,]
    currV       <- v.post[numV,]
    curr_groupU <- curr_groupV <- c()
    currRowU <- currRowV    <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      currRowV   <- currV[k,]
      curr_groupV[k] <- which(currRowV == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1)
      temp_propV     <- sample(1:z,1)
      
      currRowU       <- currU[j,]
      currRowV       <- currV[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      curr_groupV[j] <- which(currRowV ==1)
      prop_groupV[j] <- temp_propV
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    temp_v.post                <- v.post
    
    temp_v.post[numV,]          <- rep(0,z)
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      temp_v.post[numV[h],prop_groupV[h]] <- 1
      
      
    }
    
    
    accept_prob      <-  probitLikelihood(temp_u.post,temp_v.post, lambda, 
                                          XB+outer(a,b, "+") ,Z2,Y,c(numU,numV),rho,td,to)
    
    accept_prob_past <- probitLikelihood(u.post,v.post, lambda,
                                         XB+ outer(a,b,"+"), Z2,Y,c(numU,numV),rho,td,to)
    
    
    
    
    runi <- runif(1)
    ratio_pro <- accept_prob -accept_prob_past
    
    print(ratio_pro)
    ###################################################3 accept reject
    if( ratio_pro > log(runi)) {
      
      ###################################################3 accept reject
      u.post <- temp_u.post
      v.post <- temp_v.post
      updateU <- TRUE
      swapPeep <- c(numU, numV)
      MH <- MH + 1
    }else{
      u.post <- u.post
      v.post <- v.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    counter <- counter + 1
  }else{
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- c()
    currU       <- u.post[numU,]
    curr_groupU <-  c()
    currRowU   <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1)
      
      currRowU       <- currU[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      
      
    }
    
    
    
    ########################################### Compare
    
    
    accept_prob      <-  probitLikelihood(temp_u.post,temp_u.post, lambda, 
                                          XB+outer(a,b, "+") ,Z2,Y,numU,rho,td,to)
    
    accept_prob_past <- probitLikelihood(u.post,u.post, lambda,
                                         XB+ outer(a,b,"+"), Z2,Y,numU,rho,td,to)
    
    
    
    
    runi <- runif(1)
    ratio_pro <- accept_prob -accept_prob_past
    
    print(ratio_pro)
    ###################################################3 accept reject
    if( ratio_pro > log(runi)) {
      u.post <- temp_u.post
      v.post <- temp_u.post
      updateU <- TRUE
      swapPeep <- c(numU)
      MH <- MH + 1
      print("RATIOS")
      
    }else{
      u.post <- u.post
      v.post <- u.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    
    counter <- counter + 1
  }
  ## pick random people to switch
  
  
  return(list(U = U, V = V, updateU = updateU,
              swapPeep = swapPeep, 
              MH =MH, counter = counter))
}





# Function to calculate the censored probit likelihood
probitLikelihood_cbin <- function(U, V, lambda, XB, predY,Y, numU,rho,td,to, c, Z){
  ULV <- U%*%lambda%*%t(V)
  n <- dim(Y)[1]
  meanMat <- XB + ULV
  
  
  # initalize likelihood
  lh <-0
  
  for(i in numU){
    for(j in 1:n){
      if(i != j){
        
        if(Y[i,j]==1 && Y[j,i] == 1){
          
          
          tmp3 <-  pbivnorm(meanMat[i,j]-c[i], y = meanMat[j,i]-c[j], rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          } 
          
          lh <- lh + log(tmp3)
          #lh<-lh+log(tmp1) +log(tmp2)
        }else if(Y[i,j]==1 && Y[j,i] == 0){
          
          tmp3 <-  pbivnorm(meanMat[i,j]- c[i], y = -meanMat[j,i] + c[j], -rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }
          lh <- lh + log(tmp3)
          
        }else if(Y[i,j]==0 && Y[j,i] == 1){
          
          
          tmp3 <-  pbivnorm(-meanMat[i,j] + c[i], y = meanMat[j,i]-c[j], -rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }
          lh <- lh + log(tmp3)
          
        }else{
          tmp3 <-  pbivnorm(-meanMat[i,j] + c[i], y = -meanMat[j,i] + c[j], rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }   
          lh <- lh + log(tmp3)
        }
        
      }
    }
  }
  
  
  # get diag entries
  tmp1 <- pnorm(diag(meanMat)[numU], sd =sqrt(1 + rho))
  
  tmp1 <- ifelse(tmp1 == 0, .000000000001, tmp1)
  tmp1 <- ifelse(tmp1 == 1, .999999999999, tmp1)
  lh <-lh+ sum(log(1-tmp1))
  
  return(lh)
}

# Function to draw U, U/V
# Function to calculate dep X beta that is NOT all added together (returns array of XB)
#' @param lambda: current lambda (k x k matrix)
#' @param U: current U (n x k matrix)
#' @param V: current V (n x k matrix)
#' @param XB: current XB value
#' @param a,b: n vectors with current row/col random effects
#' @param Y: Observed network
#' @param UandV: TRUE/FALSE (are both estimated or just U)
#' @param counter: counts total iterations with burnin
#' @param MH: counts acceptance of MH step
#' @param X_r: n x n x p_dep array of row covar
#' @param X_c: n x n x p_dep array of col covar
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param beta_r: current beta_r, p_dep x k matrix
#' @param beta_c: current beta_c, p_dep x k matrix
#' @param alpha: current alpha value
#' @param beta_dr: current beta_dr, p_dep_dyad x k matrix
#' @param beta_dc: current beta_dc, p_dep_dyad x k matrix
#' @param h_r_past: Previous calculation of h_r
#' @param h_c_past: Previous calculation of h_c
#' @param Z: current Z value (n x n matrix)
#' @param td, to: current de-correlation values
#' @param Sab: covar from a,b
#' @param rho: current rho value
#' @param indepBeta: indicates indep covar
#' @param XB_indep_get: XB just for indep covariates

drawUV_probitLikelihood_indepDep <- function (lambda, U, V,  XB, a,b,Y, UandV,counter, MH,
                                              X_r, X_c,X_d, beta_r, beta_c, alpha,
                                              beta_dr, beta_dc,
                                              h_r_past, h_c_past,Z,td,to,Sab,rho,indepBeta, XB_indep_get) {
  z <- dim(U)[2]
  
  #y.post    <- Y[lower.tri(Y, diag = FALSE)]
  
  y.post <-Y
  Z2 <- Z
  u.post <- U
  v.post <- V
  
  
  ####################################################################### Metropolis
  
  
  
  if(UandV == TRUE){
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    numV  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    oldV <-  apply(V[numV,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- prop_groupV <- c()
    currU       <- u.post[numU,]
    currV       <- v.post[numV,]
    curr_groupU <- curr_groupV <- c()
    currRowU <- currRowV    <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      currRowV   <- currV[k,]
      curr_groupV[k] <- which(currRowV == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1)
      temp_propV     <- sample(1:z,1)
      
      currRowU       <- currU[j,]
      currRowV       <- currV[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      curr_groupV[j] <- which(currRowV ==1)
      prop_groupV[j] <- temp_propV
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    temp_v.post                <- v.post
    
    temp_v.post[numV,]          <- rep(0,z)
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      temp_v.post[numV[h],prop_groupV[h]] <- 1
      
      
    }
    ########################################### Compare
    
    XB_ind_getP <- XBETA_DEP_allowIndep(X_r[,,-indepBeta], X_c[,,-indepBeta], beta_r,
                                        beta_c, temp_u.post, dim(beta_r)[1],
                                        TRUE, h_r_past, h_c_past,indepBeta) 
    
    
    XBETA_ROWP     <- XB_ind_getP$XBETAROW
    XBETA_COLP     <- XB_ind_getP$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYADP    <- rowSums(XB_ind_getP$XBETADYAD,dims = 2)
    }else{
      XBETA_DYADP <- 0
    }
    
    XB_indP        <- XBETA_ROWP + XBETA_COLP 
    
    if(dim(beta_r)[1]==1){
      
      XBP            <- XB_indP + alpha + XBETA_DYADP + XB_indep_get # calculate XB
    }else{
      
      XBP            <- rowSums(XB_indP, dims =2 )+ alpha +XB_indep_get +
        XBETA_DYADP # calculate XB
    }
    
    accept_prob      <-  probitLikelihood(temp_u.post,temp_v.post, lambda, 
                                          XBP+outer(a,b, "+") ,Z2,Y,c(numU,numV),rho,td,to)
    
    accept_prob_past <- probitLikelihood(u.post,v.post, lambda,
                                         XB+ outer(a,b,"+"), Z2,Y,c(numU,numV),rho,td,to)
    
    
    runi <- runif(1)
    ratio_pro <- accept_prob -accept_prob_past
    
    
    ################################################### accept reject
    if( ratio_pro > log(runi)) {
      u.post <- temp_u.post
      v.post <- temp_v.post
      updateU <- TRUE
      swapPeep <- c(numU, numV)
      MH <- MH + 1
    }else{
      u.post <- u.post
      v.post <- v.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    counter <- counter + 1
  }else{
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- c()
    currU       <- u.post[numU,]
    curr_groupU <-  c()
    currRowU   <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1)
      
      currRowU       <- currU[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      
      
    }
    
    p_dep <-ifelse(is.null(dim(beta_r)[1]), 1, dim(beta_r)[1])
    
    XB_ind_getP    <- XBETA_DEP_partition(X_r[,,-indepBeta], X_c[,,-indepBeta],
                                          X_d, beta_r, beta_c,
                                          beta_dr, beta_dc, temp_u.post,temp_u.post,
                                          updatedUV=TRUE, h_r_past, h_c_past)
    XBETA_ROWP     <- XB_ind_getP$XBETAROW
    XBETA_COLP     <- XB_ind_getP$XBETACOL
    
    if(!is.null(X_d)){
      
      XBETA_DYADP    <- rowSums(XB_ind_getP$XBETADYAD,dims = 2)
    }else{
      XBETA_DYADP <- 0
    }
    
    XB_indP        <-  XBETA_ROWP+XBETA_COLP
    
    if(p_dep==1){
      
      XBP            <- XB_indP + alpha + XBETA_DYADP + XB_indep_get # calculate XB
    }else{
      
      XBP            <- rowSums(XB_indP, dims =2 )+ alpha +XB_indep_get +
        XBETA_DYADP # calculate XB
    }
    
    
    
    
    
    accept_prob      <-  probitLikelihood(temp_u.post,temp_u.post, lambda, 
                                          XBP+outer(a,b, "+") ,Z2,Y,numU,rho,td,to)
    
    accept_prob_past <- probitLikelihood(u.post,u.post, lambda,
                                         XB+ outer(a,b,"+"), Z2,Y,numU,rho,td,to)
    
    
    runi <- runif(1)
    ratio_pro <- accept_prob -accept_prob_past
    
    print(ratio_pro)
    ###################################################3 accept reject
    if( ratio_pro > log(runi)) {
      u.post <- temp_u.post
      v.post <- temp_u.post
      updateU <- TRUE
      swapPeep <- c(numU)
      MH <- MH + 1
      print("RATIOS")
      
    }else{
      u.post <- u.post
      v.post <- u.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    
    counter <- counter + 1
  }
  ## pick random people to switch
  
  
  return(list(U = U, V = V, updateU = updateU,
              swapPeep = swapPeep, 
              MH =MH, counter = counter))
}

# drawing UV for our comm dep model no censoring


# Function to calc. probit likelihood
probitLikelihood <- function(U, V, lambda, XB, Z2, Y,numU,rho,td,to){
  ULV <- U%*%lambda%*%t(V)
  n <- dim(Z2)[1]
  meanMat <- XB + ULV
  
  lh <-0
  # what we care about...
  
  for(i in numU){
    for(j in 1:n){
      if(i !=j){
        if(Y[i,j]==1 && Y[j,i] == 1){
          
          
          tmp3 <-  pbivnorm(meanMat[i,j], y = meanMat[j,i], rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }   
          lh <- lh + log(tmp3)
          #lh<-lh+log(tmp1) +log(tmp2)
        }else if(Y[i,j]==1 && Y[j,i] == 0){
          
          
          tmp3 <-  pbivnorm(meanMat[i,j], y = -meanMat[j,i], -rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }          
          lh <- lh + log(tmp3)
          
          # lh<-lh+log(tmp1) +log(1-tmp2)
          
        }else if(Y[i,j]==0 && Y[j,i] == 1){
          
          tmp3 <-  pbivnorm(-meanMat[i,j], y = meanMat[j,i], -rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }           
          
          lh <- lh + log(tmp3)
          
        }else{
          tmp3 <-  pbivnorm(-meanMat[i,j], y = -meanMat[j,i], rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }   
          lh <- lh + log(tmp3)
          
          # lh <- lh+log(1-tmp1) +log(1-tmp2)
          
        }
        # print(tmp3);print(i);print(j); print(Y[i,j])
      }
    }
    
  }
  
  tmp1 <- pnorm(diag(meanMat)[numU], sd =sqrt(1 + rho))
  
  tmp1 <- ifelse(tmp1 == 0, .000000000001, tmp1)
  tmp1 <- ifelse(tmp1 == 1, .999999999999, tmp1)
  lh <-lh+ sum(log(1-tmp1))
  
  return(lh)
}


# Function to draw U, U/V
# Function to calculate dep X beta that is NOT all added together (returns array of XB)
#' @param lambda: current lambda (k x k matrix)
#' @param U: current U (n x k matrix)
#' @param V: current V (n x k matrix)
#' @param XB: current XB value
#' @param a,b: n vectors with current row/col random effects
#' @param Y: Observed network
#' @param UandV: TRUE/FALSE (are both estimated or just U)
#' @param counter: counts total iterations with burnin
#' @param MH: counts acceptance of MH step
#' @param X_r: n x n x p_dep array of row covar
#' @param X_c: n x n x p_dep array of col covar
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param beta_r: current beta_r, p_dep x k matrix
#' @param beta_c: current beta_c, p_dep x k matrix
#' @param alpha: current alpha value
#' @param beta_dr: current beta_dr, p_dep_dyad x k matrix
#' @param beta_dc: current beta_dc, p_dep_dyad x k matrix
#' @param h_r_past: Previous calculation of h_r
#' @param h_c_past: Previous calculation of h_c
#' @param Z: current Z value (n x n matrix)
#' @param td, to: current de-correlation values
#' @param Sab: covar from a,b
#' @param rho: current rho value

drawUV_probitlLikelihood <- function (lambda, U, V,  XB, a,b,Y, UandV,counter, MH,
                                      X_r, X_c,X_d, beta_r, beta_c, alpha,
                                      beta_dr, beta_dc,
                                      h_r_past, h_c_past,Z,td,to,Sab,rho) {
  z <- dim(U)[2]
  
  #y.post    <- Y[lower.tri(Y, diag = FALSE)]
  
  y.post <-Y
  Z2 <- Z
  u.post <- U
  v.post <- V
  
  
  ####################################################################### Metropolis
  
  
  
  
  if(UandV == TRUE){
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    numV  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    oldV <-  apply(V[numV,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- prop_groupV <- c()
    currU       <- u.post[numU,]
    currV       <- v.post[numV,]
    curr_groupU <- curr_groupV <- c()
    currRowU <- currRowV    <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      currRowV   <- currV[k,]
      curr_groupV[k] <- which(currRowV == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1)
      temp_propV     <- sample(1:z,1)
      
      currRowU       <- currU[j,]
      currRowV       <- currV[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      curr_groupV[j] <- which(currRowV ==1)
      prop_groupV[j] <- temp_propV
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    temp_v.post                <- v.post
    
    temp_v.post[numV,]          <- rep(0,z)
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      temp_v.post[numV[h],prop_groupV[h]] <- 1
      
      
    }
    ########################################### Compare
    XB_ind_getP    <- XBETA_DEP_partition(X_r, X_c,X_d, beta_r, beta_c,
                                          beta_dr, beta_dc, temp_u.post,
                                          temp_v.post,
                                          TRUE, h_r_past, h_c_past)
    
    XBETA_ROWP     <- XB_ind_getP$XBETAROW
    XBETA_COLP     <- XB_ind_getP$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYADP    <- rowSums(XB_ind_getP$XBETADYAD,dims = 2)
    }else{
      XBETA_DYADP <- 0
    }
    
    XB_indP        <- XBETA_ROWP + XBETA_COLP 
    XBP            <- rowSums(XB_indP, dims =2 )+ alpha +
      XBETA_DYADP # calculate XB
    
    
    accept_prob      <-  probitLikelihood(temp_u.post,temp_v.post, lambda, 
                                          XBP+outer(a,b, "+") ,Z2,Y,c(numU,numV),rho,td,to)
    
    accept_prob_past <- probitLikelihood(u.post,v.post, lambda,
                                         XB+ outer(a,b,"+"), Z2,Y,c(numU,numV),rho,td,to)
    
    
    
    
    runi <- runif(1)
    ratio_pro <- accept_prob -accept_prob_past
    
    print(ratio_pro)
    ###################################################3 accept reject
    if( ratio_pro > log(runi)) {
      
      ###################################################3 accept reject
      u.post <- temp_u.post
      v.post <- temp_v.post
      updateU <- TRUE
      swapPeep <- c(numU, numV)
      MH <- MH + 1
    }else{
      u.post <- u.post
      v.post <- v.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    counter <- counter + 1
  }else{
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- c()
    currU       <- u.post[numU,]
    curr_groupU <-  c()
    currRowU   <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1)
      
      currRowU       <- currU[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      
      
    }
    
    XB_ind_getP    <- XBETA_DEP_partition(X_r, X_c,X_d, beta_r, beta_c,
                                          beta_dr, beta_dc, temp_u.post,
                                          temp_u.post,
                                          TRUE, h_r_past, h_c_past)
    
    XBETA_ROWP     <- XB_ind_getP$XBETAROW
    XBETA_COLP     <- XB_ind_getP$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYADP    <- rowSums(XB_ind_getP$XBETADYAD,dims = 2)
    }else{
      XBETA_DYADP <- 0
    }
    
    XB_indP        <- XBETA_ROWP + XBETA_COLP 
    XBP            <- rowSums(XB_indP, dims =2 )+ alpha +XBETA_DYADP # calculate XB
    
    ########################################### Compare
    
    
    accept_prob      <-  probitLikelihood(temp_u.post,temp_u.post, lambda, 
                                          XBP+outer(a,b, "+") ,Z2,Y,numU,rho,td,to)
    
    accept_prob_past <- probitLikelihood(u.post,u.post, lambda,
                                         XB+ outer(a,b,"+"), Z2,Y,numU,rho,td,to)
    
    
    
    
    runi <- runif(1)
    ratio_pro <- accept_prob -accept_prob_past
    
    print(ratio_pro)
    ###################################################3 accept reject
    if( ratio_pro > log(runi)) {
      u.post <- temp_u.post
      v.post <- temp_u.post
      updateU <- TRUE
      swapPeep <- c(numU)
      MH <- MH + 1
      print("RATIOS")
      
    }else{
      u.post <- u.post
      v.post <- u.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    
    counter <- counter + 1
  }
  ## pick random people to switch
  
  
  return(list(U = U, V = V, updateU = updateU,
              swapPeep = swapPeep, 
              MH =MH, counter = counter))
}

### MAIN FUNCTIONS FOR COMM DEP MODEL




# Methods for estimating U to initialize

# Main function to run MCMC, all covariates are assumed to be dep. on community
# Main function to run MCMC, all covariates are assumed to be dep. on community
#' @param X_r: n x n x p_dep array of row covar
#' @param X_c: n x n x p_dep array of col covar
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param Y: Observed network (can be binary or censored binary)
#' @parm iter: number of iterations to run mcmc for 
#' @param numGroup: number of communities (also denote k)
#' @param prior_beta_mu: k dimensional vector with prior mean for beta
#' @param prior_beta_var:k x k matrix for prior variance for beta
#' @param prior_alpha_mu: constant for alpha(intercept) mean
#' @param prior_alpha_var:constant for alpha(intercept) var
#' @param start_beta_c: p_dep x k matrix for initializing beta_c
#' @param start_beta_r:  p_dep x k matrix for initializing beta_c
#' @param start_beta_dyad:  p_dep_dyad x k matrix for initializing beta_c
#' @param start_alpha: constant to initialize alpha
#' @param UandV: Default is FALSE, indicates if U and V should be calculated or just U
#' @param symmetricLambda: Default is TRUE, indicates if lambda should be symmetric
#' @param keep.UV: Default is TRUE, store U/V?
#' @param dcor: Default is TRUE, indicates if correlation in error terms should be estimated
#' @param symmetric: Default is FALSE, indicates that Y is not symmetric
#' @param model:  "cbin" or "bin", indicates what type of network Y is 
#' @param odmax: if model = "cbin", this is the max out degree a node can have
#' @param odens = 10: how often to save output
#' @param burnin: amount of burnin iterations
#' @param badInit: start with completely random initialization for U/V
#' @param prior_group_probs: For future use (if we want to change prior prob on group membership)
library(MCMCpack)

subAmen_allDepBeta_censored <-function(X_r, X_c, X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                       prior_alpha_mu, prior_alpha_var, start_beta_c,
                                       start_beta_r,start_beta_dyad, start_alpha, 
                                       UandV = FALSE, symmetricLambda = TRUE,  keep.UV = TRUE,
                                       dcor = TRUE, symmetric =FALSE, model = "cbin",odmax = NULL,odens = 10, 
                                       burnin = 1000,badInit = FALSE,prior_group_probs = NULL){
  
  print("CRNS")
  
  p <- dim(X_r)[3]
  n <- dim(X_r)[1]
  ###################### sampling
  beta_r <- start_beta_r
  beta_c <- start_beta_c
  
  if(!is.null(X_d)){
    beta_dc <- beta_dr <- start_beta_dyad
    
  }else{
    beta_dc <- beta_dr <- 0
    
  }
  
  
  if(symmetric == FALSE){
    rvar <- cvar <-  TRUE
    nvar <- FALSE
  }else{
    rvar <-cvar <- TRUE
    nvar <- TRUE
  }
  counter <- 0
  MH <-0
  
  
  Y_binary <- 1 * (Y > 0)
  
  
  if (is.element(model, c("cbin", "frn", "rrl"))) {
    odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
    if (length(odmax) == 1) {
      odmax <- rep(odmax, nrow(Y))
    }
  }
  
  sigma_c <- 1
  sigma_c_alpha <- 10; sigma_c_beta <- 100
  
  diag(Y) <- NA
  
  if (is.element(model, c("cbin", "frn"))) {
    Z <- Y
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      zi <- zscores(yi)
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
        zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 0.001
      }
      if (length(rnkd) < odmax[i]) {
        urnkd <- which(!is.na(yi) & yi == 0)
        if (max(zi[urnkd]) > 0) {
          zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 0.001
        }
      }
      Z[i, ] <- zi
    }
  }
  
  diag(Y) <- 0
  
  ######
  ## initialize U using regularizer from Rohe
  diag(Y_binary) <- 0
  UV_start <- specClust(Y_binary, tauAvg = TRUE,  numGroup,UandV)
  U <- V <- matrix(UV_start$U_hat, ncol = numGroup)
  image(U %*%t(U))
  #V <-  matrix(UV_start$V_hat, ncol = numGroup)
  
  if(badInit == TRUE){
    set.seed(13)
    k<- dim(U)[2]
    for(i in 1:n){
      membership <- sample(1:k, 1)
      U[i,membership] <- 1
      U[i, seq(1,k)[-membership]] <- 0
    }
    V <-U
    image(U%*%t(U))
  }
  
  ## initialize lambda
  prior.lambda.mu  <- rep(1, numGroup^2)
  prior.lambda.var <- prior.lambda.prec <- diag(1, numGroup^2) ## have had .2
  lambda <- matrix(prior.lambda.mu, nrow = numGroup)
  
  
  ## alpha
  prior_alpha     <-  prior_alpha_mu
  prior_var_alpha <- prior_alpha_var
  alpha           <- start_alpha
  
  # get initial XBeta
  h_r_past <- NULL
  h_c_past <- NULL
  XB_ind_get <- XBETA_DEP(X_r, X_c, X_d, beta_r, 
                          beta_c, beta_dr, beta_dc, 
                          U, V,TRUE, h_r_past, h_c_past) 
  
  XB_ind     <- XB_ind_get$XBETA
  
  if(!is.null(X_d)){
    XB_dyad    <- XB_ind_get$XBETADYAD
    XB         <- rowSums(XB_dyad, dims = 2) + rowSums(XB_ind, dims =2 )+ alpha # calculate XB
    
  }else{
    XB         <-  rowSums(XB_ind, dims =2 )+ alpha # calculate XB
    
  }
  
  # initialize censoring
  censored_individuals <- which(rowSums(Y_binary) == odmax)
  

  pred_Y <- Y_binary
  c <- rep(0,n)
  c[censored_individuals] <- -.0001
  E <- Z- XB - U %*% lambda %*%t(V) - c
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  # from AMEN code (Hoff et al)
  prior <- list()
  if(model!="nrm")
  {  
    ydist<-table(Y)
    ymode<-as.numeric(names(ydist)[ ydist==max(ydist) ])[1] 
 
    YB<-1*(Y!=ymode) 
    ybar<-mean(YB,na.rm=TRUE) ; mu<-qnorm(ybar)
    E<- (YB - ybar)/dnorm(qnorm(ybar)) ; diag(E)<-0
    a<-rowMeans(E,na.rm=TRUE)  ; a[is.na(a)]<-0 
    b<-colMeans(E,na.rm=TRUE)  ; b[is.na(b)]<-0
    vscale<-mean(diag(cov(cbind(a,b))))
    PHAT<-pnorm(mu+outer(a,b,"+"))
    vdfmlt<-.25/mean(PHAT*(1-PHAT))
    if(is.null(prior$Sab0)){ prior$Sab0<-diag(2)*vscale }
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$kappa0)){ prior$kappa0<-round((2*numGroup+2)*vdfmlt) }  
    if(is.null(prior$g)){ prior$g<-sum(!is.na(Y)) }
  }
  
  Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar)) 
  E <- E- outer(a,b, "+")
  s2<-1
  if(model=="nrm"){s2<-mean(E^2)}
  rho<-cor( c(E[upper.tri(E)]), c(t(E)[upper.tri(E)]) )*dcor  
  
  ## initialize things to store
  for(i in 1: dim(X_r)[3]){
    nam <- paste("beta_rMat", i, sep = "")
    assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
    
    namc <- paste("beta_cMat", i, sep = "")
    assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
    
    
    
  }
  
  if(!is.null(X_d)){
    for(i in 1:dim(X_d)[3]){
      
      nam <- paste("beta_rMatDYAD", i, sep = "")
      assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
      
      namc <- paste("beta_cMatDYAD", i, sep = "")
      assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
      
    }
  }
  
  ulv <- U %*% lambda %*% t(V)
  lambda_mat <- matrix(nrow = iter/odens, ncol = numGroup^2)
  if(keep.UV == TRUE){
    UV_mat      <-  matrix(nrow = iter/odens, ncol = length(c(U)))
    
    
  }
  
  
  
  if(is.null(prior_group_probs)){
    prior_group_probs <- rep(1/numGroup, numGroup)
  }
  
  Z_mat <- matrix(nrow = iter/odens, ncol = n)
  
  
  zmean <- matrix(nrow = iter/odens, ncol = 1)
  Ys <-matrix(nrow = iter/odens, ncol = 1)
  
  updatedUV <- TRUE
  XB_ind_get    <- XBETA_DEP_partition(X_r, X_c, X_d, beta_r, beta_c,beta_dr, beta_dc, U,V,updatedUV, h_r_past, h_c_past)
  XBETA_ROW     <- XB_ind_get$XBETAROW
  XBETA_COL     <- XB_ind_get$XBETACOL
  
  if(!is.null(X_d)){
    
    XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
  }else{
    XBETA_DYAD <- 0
  }
  XB_ind        <- XBETA_ROW + XBETA_COL 
  h_r_past      <- XB_ind_get$h_r_past
  h_c_past      <- XB_ind_get$h_c_past
  
  
  XB            <- rowSums(XB_ind, dims =2 )+ XBETA_DYAD + alpha # calculate XB
  
  EZ <- XB + ulv  + outer(a, b, "+") +c
  updatedUV <- TRUE
  diag(Y_binary) <-0

  ## start sampling
  diag(Z) <- mean(na.omit(c(Z)))
  Se   <-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2 <-mhalf(solve(Se))
  td   <-iSe2[1,1] ; to<-iSe2[1,2]
  
  
  # j = 0
  for(i in 1:(iter+burnin)){
    
    
    # Draw latent lambda and U
    lambda <- drawLambda(Z, XB + outer(a,b, "+") + c, U,V, prior.lambda.var,
                         prior.lambda.mu,td, to,  symmetricLambda = symmetricLambda)
    
    
    
    UV_info      <- drawUV_probitLikelihood_ypred_cbin(lambda, U, V,  XB ,a,b,Y,
                                                       c, UandV,counter, MH,
                                                       X_r, X_c,X_d, beta_r, beta_c,alpha,
                                                       beta_dr, beta_dc,
                                                       h_r_past, h_c_past,Z,td,to,Sab, rho)#, prior_group_probs)
    
    
    U            <- UV_info$U
    V            <- UV_info$V
    
    
    print(colSums(U))
    
    updatedUV     <- UV_info$updateU
    swapPeep     <- UV_info$swapPeep
    counter      <- UV_info$counter
    MH           <- UV_info$MH
    
    
    ULV <- U %*%lambda%*%t(V)
    
    
    
    
    XB_ind_get    <- XBETA_DEP_partition(X_r, X_c,X_d, beta_r, beta_c,
                                         beta_dr, beta_dc, U,V,
                                         updatedUV, h_r_past, h_c_past)
    XBETA_ROW     <- XB_ind_get$XBETAROW
    XBETA_COL     <- XB_ind_get$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
    }else{
      XBETA_DYAD <- 0
    }
    
    XB_ind        <- XBETA_ROW + XBETA_COL 
    h_r_past      <- XB_ind_get$h_r_past
    h_c_past      <- XB_ind_get$h_c_past
    XB            <- rowSums(XB_ind, dims =2 )+ alpha +XBETA_DYAD # calculate XB
    #
    EZ <- XB + ULV + outer(a, b, "+") + c 
    
    
    if(updatedUV == TRUE){
      for(h in 1:length(swapPeep)){
        Z <-    update_Z_bin(swapPeep[h], Y, EZ, rho,  n,Z)
      }
      
    }
    
    
    # Draw rho
    
    if(symmetric){rho <- min(0.9999, 1 - 1/sqrt(i))}
    if(dcor){ rho <- rrho_mh_new(Z-EZ,rho,s2,offset=0,asp = NULL)} 
    Se   <-matrix(c(1,rho,rho,1),2,2)*s2
    iSe2 <-mhalf(solve(Se))
    td   <-iSe2[1,1] ; to<-iSe2[1,2]
    
    
    ## draw intercept
    XB_tmp <- XB -alpha
    
    alpha_tmp <-  rALPHA(Z-c, ULV,
                         alpha, a , b, XB,td, to,prior_alpha_mu,prior_alpha_var)
    alpha <- alpha_tmp$alpha 
    
    XB_updatedALPHA <- XB_tmp + alpha
    
    if(!is.null(X_d)){
      
      
      
      numEigs <- 1 ### CHANGE
      
      tmpDYAD <- rbeta_DECORR_DEPDYAD(Z-XB_updatedALPHA+XBETA_DYAD-c, X_d, U,V, lambda, 
                                      beta_dr, beta_dc, 
                                      prior_beta_mu, 
                                      prior_beta_var , a, b,   to ,td,numEigs, h_r_past,
                                      h_c_past)
      
      
      beta_dr <- tmpDYAD$beta_row_dyad
      beta_dc <- tmpDYAD$beta_col_dyad
      print("DYADF")
      print(tmpDYAD)
      
      XBtmp <- XB_updatedALPHA -XBETA_DYAD
      XB_updateDyad_get    <- XBETA_DEP_UPDATEDDYAD(X_d, 
                                                    beta_dr, beta_dc,
                                                    U,V)
      XBETA_DYAD    <- rowSums(XB_updateDyad_get$XBETADYAD, dims =2 )
      XB            <-  XBtmp +XBETA_DYAD # calculate XB
      EZ <- XB + ULV + outer(a, b, "+") + c
      
      
    }else{
      XB <- XB_updatedALPHA
      EZ <- XB + ULV + outer(a, b, "+") + c
      
    }
    #
    
    sigma_c <- draw_covar_c(sigma_c_alpha, sigma_c_beta, 0, c, n)
    # draw C
    c <- rC(Z, ULV, a, b, XB, td,to,c, censored_individuals, prior_c = rep(0,n), 
            sigma_c)
    

    # Draw Beta
    
    
    tmp <- rbeta_DECORR_DEP(Z, X_r, X_c,XBETA_DYAD,U,V, lambda, dim(X_r)[3], beta_r, beta_c,
                            alpha + c, prior_beta_mu, 
                            prior_beta_var, a, b,  h_r_past, h_c_past, to ,td,updatedUV )
    
    
    
    beta_r <- tmp$beta_row
    beta_c <- tmp$beta_col
    
    
    XBtmp <- XB - rowSums(XB_ind, dims =2 )
    
    # recalculate XB
    
    
    XB_updateRC <- XBETA_DEP_UPDATEROWCOLUMN(X_r, X_c, beta_r, beta_c, 
                                             U,V,
                                             FALSE, 
                                             h_r_past, 
                                             h_c_past)
    
    XBETA_ROW     <-  XB_updateRC$XBETAROW
    XBETA_COL     <-  XB_updateRC$XBETACOL
    XB_ind        <-  XBETA_ROW + XBETA_COL 
    XB            <-  XBtmp + rowSums(XB_ind, dims =2)
    #
    EZ <- XB + ULV + outer(a, b, "+") + c
    
    
    Z  <- rZ_bin_fc_new(Z, EZ, rho, Y)
    
    
    # draw random effects
    ab <- r_ab_DECORR_DEP(Z, XB + c, Sab,rho, s2, ULV)
    print(summary(c(Z))); print(summary(c(EZ)))
    a  <- ab$a * rvar
    b  <- ab$b * cvar
    
    if(symmetric == TRUE){
      a <- b <- (a + b)/2
    }
    EZ <- XB + ULV + outer(a, b, "+") + c
    
    # draw Sab and a -> z
    if(symmetric == FALSE){
      Sab <-rSab_fc(a,b,Sab0=prior$Sab0,eta0=prior$eta0)
      
      if(model=="bin"){
        tmp<-raSab_bin_fc_new(Z,Y,a,b,Sab,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      if(model=="cbin")
      {
        # tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs)
        tmp<-raSab_bin_fc_new(Z,Y,a,b,Sab,Sab0=prior$Sab0,eta0=prior$eta0)
        
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      
    } else{
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]
    }
    
    EZ <- XB + ULV + outer(a, b, "+") + c
    
    
    # posterior pred
    if(symmetric){ EZ<-(EZ+t(EZ))/2 } 
    
    
    
    if(model=="cbin"){ Ys_draw<-1*(simY_frn(EZ,rho,odmax,YO=Y)>0) }
    
    zsim <- Z #simZ(EZ, rho)
    Ysim <- 1*(zsim > 0)
    replacement_Y <- Ysim[censored_individuals,]
    
    pred_Y <- Y_binary
    pred_Y[censored_individuals, ] <- replacement_Y
    print("PRED Y")
    print(all(Ysim == pred_Y))
    
    degreeCounter <- rowSums(1*(zsim > 0))
    # Store results
    if(i == burnin + 1){ j = 0}
    if(i > burnin){
      
      if(i %% odens==0){
        j = j + 1
        
        for(g in 1:p){
          nam <-sprintf("beta_rMat%s[%s,]", g, j)
          nam2 <-sprintf("beta_cMat%s[%s,]", g, j)
          
          eval(parse(text= paste(nam, "<- beta_r[g,]")))        
          eval(parse(text= paste(nam2, "<- beta_c[g,]")))        
          
          
          
          
          
          
        }
        if(!is.null(X_d)){
          
          for(k in 1:dim(X_d)[3]){
            nam3 <-sprintf("beta_rMatDYAD%s[%s,]", k, j)
            nam4 <-sprintf("beta_cMatDYAD%s[%s,]", k, j)
            
            eval(parse(text= paste(nam3, "<- beta_dr[k,]")))        
            eval(parse(text= paste(nam4, "<- beta_dc[k,]")))        
            
          }
        }
        lambda_mat[j,] <- c(lambda)
        zmean[j,] <- mean(c(Z))
        Ys[j,] <-sum(abs(Ys_draw-Y),na.rm = TRUE)
        
        if(keep.UV == TRUE){
          UV_mat[j, ] <- c(U)
        }
        Z_mat[j, ] <- degreeCounter
      }
      
    }
    
    if(i %% 10000 ==0){
      results <- list()
      
      h = 1
      k=1
      while(h < 2*dim(X_r)[3]){
        
        nam <- paste("beta_rMat", k, sep = "")
        nam2 <- paste("beta_cMat", k, sep = "")
        
        results[[h]]   <- eval(parse(text= nam))
        names(results)[h] <- nam
        
        results[[h+1]] <- eval(parse(text= nam2))
        names(results)[h+1] <- nam2
        
        h <- h + 2
        k <-k+1
      }
      
      
      if(!is.null(X_d)){
        h = 1
        k=1
        while(h < 2*dim(X_d)[3]){
          
          nam3 <- paste("beta_rMatDYAD", k, sep = "")
          nam4 <- paste("beta_cMatDYAD", k, sep = "")
          
          results[[(2*dim(X_r)[3] + h)]]   <- eval(parse(text= nam3))
          names(results)[(2*dim(X_r)[3] + h)] <- nam3
          
          results[[(2*dim(X_r)[3] + h+1)]] <- eval(parse(text= nam4))
          names(results)[(2*dim(X_r)[3] + h +1)] <- nam4
          
          h <- h + 2
          k <-k+1
        }
        
        
      }
      if(keep.UV == TRUE){
        results[[length(results)+1]] <- UV_mat
        names(results)[[length(results)]] <- "U"
      }
      results[[length(results)+1]] <-   Z_mat 
      names(results)[[length(results)]] <- "degcount"
      
      results[[length(results)+1]] <-   lambda_mat
      names(results)[[length(results)]] <- "lambda"
      # saveRDS(results, sprintf("/hpc/group/volfovskylab/hmm40/CBIN/results_cbin_disprop_comm_censored_AMEN_333_10_%s_OCT%s%s.rds","givenName", n,i))
      saveRDS(results, sprintf("/hpc/group/volfovskylab/hmm40/CBIN/results_cbin_disprop_comm_censored_AMEN_333_10_%s_JLY21%s%s.rds","givenName", n,i))
      
    }
    
  }
  
  
  ## make list to return
  
  results <- list()
  
  h = 1
  k=1
  while(h < 2*dim(X_r)[3]){
    
    nam <- paste("beta_rMat", k, sep = "")
    nam2 <- paste("beta_cMat", k, sep = "")
    
    results[[h]]   <- eval(parse(text= nam))
    names(results)[h] <- nam
    
    results[[h+1]] <- eval(parse(text= nam2))
    names(results)[h+1] <- nam2
    
    h <- h + 2
    k <-k+1
  }
  
  
  if(!is.null(X_d)){
    h = 1
    k=1
    while(h < 2*dim(X_d)[3]){
      
      nam3 <- paste("beta_rMatDYAD", k, sep = "")
      nam4 <- paste("beta_cMatDYAD", k, sep = "")
      
      results[[(2*dim(X_r)[3] + h)]]   <- eval(parse(text= nam3))
      names(results)[(2*dim(X_r)[3] + h)] <- nam3
      
      results[[(2*dim(X_r)[3] + h+1)]] <- eval(parse(text= nam4))
      names(results)[(2*dim(X_r)[3] + h +1)] <- nam4
      
      h <- h + 2
      k <-k+1
    }
    
    
  }
  if(keep.UV == TRUE){
    results[[length(results)+1]] <- UV_mat
    names(results)[[length(results)]] <- "U"
  }
  results[[length(results)+1]] <-   Z_mat 
  names(results)[[length(results)]] <- "degcount"
  
  results[[length(results)+1]] <-   lambda_mat
  names(results)[[length(results)]] <- "lambda"
  
  
  print("COUNTER")
  print(counter/MH)
  print(counter)
  return(list(results=results,counter=counter, MH=MH))
  
  
}










# all given covariates are allowed to have a dependency structure
# accepts binary and fixed rank likelihoods ("frn", "bin")

## UandV (if true, then we draw membership for soc. and pop.)
## symmetricLambda (make Lambda symmetric or not)

# all given covariates are allowed to have a dependency structure
# accepts binary and fixed rank likelihoods ("frn", "bin")

## UandV (if true, then we draw membership for soc. and pop.)
## symmetricLambda (make Lambda symmetric or not)

subAmen_allDepBeta_stay_at_initialU <-function(X_r, X_c, X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                               prior_alpha_mu, prior_alpha_var, start_beta_c,
                                               start_beta_r,start_beta_dyad, start_alpha, 
                                               UandV = FALSE, symmetricLambda = TRUE,  keep.UV = FALSE,
                                               dcor = FALSE, symmetric = TRUE, model = "frn",odmax = NULL,odens = 10, 
                                               burnin = 1000,badInit = FALSE,prior_group_probs = NULL){
  
  
  
  
  
  p <- dim(X_r)[3]
  n <- dim(X_r)[1]
  ###################### sampling
  beta_r <- start_beta_r
  beta_c <- start_beta_c
  
  if(!is.null(X_d)){
    beta_dc <- beta_dr <- start_beta_dyad
    
  }else{
    beta_dc <- beta_dr <- 0
    
  }
  
  
  if(symmetric == FALSE){
    rvar <- cvar <-  TRUE
    nvar <- FALSE
  }else{
    rvar <-cvar <- TRUE
    nvar <- TRUE
  }
  counter <- 0
  MH <-0
  
  
  Y_binary <- 1 * (Y > 0)
  
  
  if (is.element(model, c("cbin", "frn", "rrl"))) {
    odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
    if (length(odmax) == 1) {
      odmax <- rep(odmax, nrow(Y))
    }
  }
  
  
  
  diag(Y) <- NA
  
  
  
  if (is.element(model, c("frn", "rrl"))) {
    
    ymx <- max(apply(1 * (Y > 0), 1, sum, na.rm = TRUE))
    YL <- NULL
    warn <- FALSE
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(yi[rnkd]) > length(unique(yi[rnkd]))) {
        warn <- TRUE
      }
      yi[rnkd] <- rank(yi[rnkd], ties.method = "random")
      Y[i, ] <- yi
      YL <- rbind(YL, match(1:ymx, yi))
    }
    if (warn) {
      cat("WARNING: Random reordering used to break ties in ranks\n")
    }
  }
  
  if (model == "bin") {
    Y <- 1 * (Y > 0)
    
    Z <- matrix(zscores(Y), nrow(Y), nrow(Y))
    z01 <- 0.5 * (max(Z[Y == 0], na.rm = TRUE) + min(Z[Y == 
                                                         1], na.rm = TRUE))
    Z <- Z - z01
  }
  if (is.element(model, c("cbin", "frn"))) {
    Z <- Y
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      zi <- zscores(yi)
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
        zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 0.001
      }
      if (length(rnkd) < odmax[i]) {
        urnkd <- which(!is.na(yi) & yi == 0)
        if (max(zi[urnkd]) > 0) {
          zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 0.001
        }
      }
      Z[i, ] <- zi
    }
  }
  
  diag(Y) <- 0
  
  ######
  ## initialize U using regularizer from Rohe
  diag(Y_binary) <- 0
  UV_start <- specClust(Y_binary, tauAvg = TRUE,  numGroup,UandV)
  U <- V <- matrix(UV_start$U_hat, ncol = numGroup)
  #image(U %*%t(U))
  #V <-  matrix(UV_start$V_hat, ncol = numGroup)
 
  
  if(badInit == TRUE){
    set.seed(13)
    k<- dim(U)[2]
    for(i in 1:n){
      membership <- sample(1:k, 1)
      U[i,membership] <- 1
      U[i, seq(1,k)[-membership]] <- 0
    }
    V <-U
    image(U%*%t(U))
  }
  
  ## initialize lambda
  prior.lambda.mu  <- rep(1, numGroup^2)
  prior.lambda.var <- prior.lambda.prec <- diag(1, numGroup^2) ## have had .2
  lambda <- matrix(prior.lambda.mu, nrow = numGroup)
  
  ## alpha
  prior_alpha     <-  prior_alpha_mu
  prior_var_alpha <- prior_alpha_var
  alpha           <- start_alpha
  
  # get initial XBeta
  h_r_past <- NULL
  h_c_past <- NULL
  XB_ind_get <- XBETA_DEP(X_r, X_c, X_d, beta_r, 
                          beta_c, beta_dr, beta_dc, 
                          U, V,TRUE, h_r_past, h_c_past) 
  
  XB_ind     <- XB_ind_get$XBETA
  
  if(!is.null(X_d)){
    XB_dyad    <- XB_ind_get$XBETADYAD
    XB         <- rowSums(XB_dyad, dims = 2) + rowSums(XB_ind, dims =2 )+ alpha # calculate XB
    
  }else{
    XB         <-  rowSums(XB_ind, dims =2 )+ alpha # calculate XB
    
  }
  
  E <- Z- XB
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  prior <- list()
  if(model!="nrm")
  {  
    ydist<-table(Y)
    ymode<-as.numeric(names(ydist)[ ydist==max(ydist) ])[1] 
    ## eg, in a sparse binary network, ymode will be zero 
    YB<-1*(Y!=ymode) 
    ybar<-mean(YB,na.rm=TRUE) ; mu<-qnorm(ybar)
    E<- (YB - ybar)/dnorm(qnorm(ybar)) ; diag(E)<-0
    a<-rowMeans(E,na.rm=TRUE)  ; a[is.na(a)]<-0 
    b<-colMeans(E,na.rm=TRUE)  ; b[is.na(b)]<-0
    vscale<-mean(diag(cov(cbind(a,b))))
    PHAT<-pnorm(mu+outer(a,b,"+"))
    vdfmlt<-.25/mean(PHAT*(1-PHAT))
    if(is.null(prior$Sab0)){ prior$Sab0<-diag(2)*vscale }
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$kappa0)){ prior$kappa0<-round((2*numGroup+2)*vdfmlt) }  
    if(is.null(prior$g)){ prior$g<-sum(!is.na(Y)) }
  }
  
  Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar)) 
  E <- E- outer(a,b, "+")
  s2<-1
  if(model=="nrm"){s2<-mean(E^2)}
  rho<-cor( c(E[upper.tri(E)]), c(t(E)[upper.tri(E)]) )*dcor  
  
  ## initialize things to store
  for(i in 1: dim(X_r)[3]){
    nam <- paste("beta_rMat", i, sep = "")
    assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
    
    namc <- paste("beta_cMat", i, sep = "")
    assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
    
    
    
  }
  
  if(!is.null(X_d)){
    for(i in 1:dim(X_d)[3]){
      
      nam <- paste("beta_rMatDYAD", i, sep = "")
      assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
      
      namc <- paste("beta_cMatDYAD", i, sep = "")
      assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
      
    }
  }
  
  ulv <- U %*% lambda %*% t(V)
  lambda_mat <- matrix(nrow = iter/odens, ncol = numGroup^2)
  if(keep.UV == TRUE){
    UV_mat      <-  matrix(nrow = iter/odens, ncol = length(c(U)))
    
    
  }
  
  
  
  if(is.null(prior_group_probs)){
    prior_group_probs <- rep(1/numGroup, numGroup)
  }
  
  Z_mat <- matrix(nrow = iter/odens, ncol = n)
  
  
  zmean <- matrix(nrow = iter/odens, ncol = 1)
  Ys <-matrix(nrow = iter/odens, ncol = 1)
  
  updatedUV <- TRUE
  XB_ind_get    <- XBETA_DEP_partition(X_r, X_c, X_d, beta_r, beta_c,beta_dr, beta_dc, U,V,updatedUV, h_r_past, h_c_past)
  XBETA_ROW     <- XB_ind_get$XBETAROW
  XBETA_COL     <- XB_ind_get$XBETACOL
  
  if(!is.null(X_d)){
    
    XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
  }else{
    XBETA_DYAD <- 0
  }
  XB_ind        <- XBETA_ROW + XBETA_COL 
  h_r_past      <- XB_ind_get$h_r_past
  h_c_past      <- XB_ind_get$h_c_past
  
  
  XB            <- rowSums(XB_ind, dims =2 )+ XBETA_DYAD + alpha # calculate XB
  
  EZ <- XB + ulv  + outer(a, b, "+")
  updatedUV <- TRUE
  diag(Y_binary) <-0
  #diag(Z) <- -1
  ## start sampling
  diag(Z) <- mean(na.omit(c(Z)))
  Se   <-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2 <-mhalf(solve(Se))
  td   <-iSe2[1,1] ; to<-iSe2[1,2]
  #if(model == "bin"){Z  <- rZ_bin_fc_new(Z, EZ, rho, Y) }
  
  # j = 0
  for(i in 1:(iter+burnin)){
    
    
    # Draw latent lambda and U
    lambda <- drawLambda(Z, XB + outer(a,b, "+"), U,V, prior.lambda.var,
                         prior.lambda.mu, td, to, symmetricLambda = symmetricLambda)
    
    # 
    # UV_info      <- drawUV_normalLik_OLD(lambda, U, V,  XB ,a,b,
    #                                      Y_binary, UandV,counter, MH,
    #                                      X_r, X_c,X_d, beta_r, beta_c,alpha,
    #                                      beta_dr, beta_dc,
    #                                      h_r_past, h_c_past,Z,td,to,Sab, rho)#, prior_group_probs)
    # 
    # UV_info <- drawUV_normalLik_alpha(lambda, U, V,  XB ,a,b,
    #                             Y_binary, UandV,counter, MH,
    #                             X_r, X_c,X_d, beta_r, beta_c,alpha,
    #                             beta_dr, beta_dc,
    #                             h_r_past, h_c_past,Z,td,to,Sab, rho)
    # 
    #  U            <- UV_info$U
    #  V            <- UV_info$V
    
    # print(UandV)
    #  print(colSums(U))
    #  print(colSums(V))
    
    
    # updatedUV     <- UV_info$updateU
    #swapPeep     <- UV_info$swapPeep
    #counter      <- UV_info$counter
    #MH           <- UV_info$MH
    
    ULV <- U %*%lambda%*%t(V)
    
    print(rho)
    
    XB_ind_get    <- XBETA_DEP_partition(X_r, X_c,X_d, beta_r, beta_c,
                                         beta_dr, beta_dc, U,V,
                                         updatedUV, h_r_past, h_c_past)
    XBETA_ROW     <- XB_ind_get$XBETAROW
    XBETA_COL     <- XB_ind_get$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
    }else{
      XBETA_DYAD <- 0
    }
    
    XB_ind        <- XBETA_ROW + XBETA_COL 
    h_r_past      <- XB_ind_get$h_r_past
    h_c_past      <- XB_ind_get$h_c_past
    XB            <- rowSums(XB_ind, dims =2 )+ alpha +XBETA_DYAD # calculate XB
    #
    EZ <- XB + ULV + outer(a, b, "+")
    
    if(model=="nrm"){ s2<-rs2_fc(Z-EZ,rho,nu0=prior$nu0,s20=prior$s20)  } 
    
    # Draw rho
    
    if(symmetric){rho <- min(0.9999, 1 - 1/sqrt(i))}
    if(dcor){ rho <- rrho_mh_new(Z-EZ,rho,s2,offset=0,asp = NULL)} 
    ## for decorrelation
    Se   <-matrix(c(1,rho,rho,1),2,2)*s2
    iSe2 <-mhalf(solve(Se))
    td   <-iSe2[1,1] ; to<-iSe2[1,2]
    
    
    ## draw intercept
    XB_tmp <- XB -alpha
    
    alpha_tmp <-  rALPHA(Z, ULV,
                         alpha, a , b, XB,td, to,prior_alpha_mu,prior_alpha_var)
    alpha <- alpha_tmp$alpha 
    
    XB_updatedALPHA <- XB_tmp + alpha
    
    if(!is.null(X_d)){
      
      
      
      numEigs <- 1 ### CHANGE
      
      tmpDYAD <- rbeta_DECORR_DEPDYAD(Z-XB_updatedALPHA+XBETA_DYAD, X_d, U,V, lambda, 
                                      beta_dr, beta_dc, 
                                      prior_beta_mu, 
                                      prior_beta_var , a, b,   to ,td,numEigs, h_r_past,
                                      h_c_past)
      
      
      beta_dr <- tmpDYAD$beta_row_dyad
      beta_dc <- tmpDYAD$beta_col_dyad
      print("DYADF")
      print(tmpDYAD)
      
      XBtmp <- XB_updatedALPHA -XBETA_DYAD
      XB_updateDyad_get    <- XBETA_DEP_UPDATEDDYAD(X_d, 
                                                    beta_dr, beta_dc,
                                                    U,V, 
                                                    FALSE, 
                                                    h_r_past, 
                                                    h_c_past)
      XBETA_DYAD    <- rowSums(XB_updateDyad_get$XBETADYAD, dims =2 )
      XB            <-  XBtmp +XBETA_DYAD # calculate XB
      EZ <- XB + ULV + outer(a, b, "+")
      
      
    }else{
      XB <- XB_updatedALPHA
      EZ <- XB + ULV + outer(a, b, "+")
      
    }
    #
    
    
    
    # Draw Beta
    
    
    tmp <- rbeta_DECORR_DEP(Z, X_r, X_c,XBETA_DYAD,U,V, lambda, dim(X_r)[3], beta_r, beta_c,
                            alpha, prior_beta_mu, 
                            prior_beta_var, a, b,  h_r_past, h_c_past, to ,td,updatedUV )
    
    print("ITER")
    print(i)
    
    beta_r <- tmp$beta_row
    beta_c <- tmp$beta_col
    
    
    XBtmp <- XB - rowSums(XB_ind, dims =2 )
    
    # recalculate XB
    
    
    XB_updateRC <- XBETA_DEP_UPDATEROWCOLUMN(X_r, X_c, beta_r, beta_c, 
                                             U,V,
                                             FALSE, 
                                             h_r_past, 
                                             h_c_past)
    
    XBETA_ROW     <-  XB_updateRC$XBETAROW
    XBETA_COL     <-  XB_updateRC$XBETACOL
    XB_ind        <-  XBETA_ROW + XBETA_COL 
    XB            <-  XBtmp + rowSums(XB_ind, dims =2)
    #
    EZ <- XB + ULV + outer(a, b, "+")
    
    
    
    
    
    
    
    # draw Z
    if(model == "frn"){Z <- rZ_frn_fc_new(Z, EZ, rho, Y, YL, odmax, odobs)}
    if(model == "bin"){Z  <- rZ_bin_fc_new(Z, EZ, rho, Y) }
    
    if(model=="nrm"){ Z<-rZ_nrm_fc(Z,EZ,rho,s2,Y) } 
    if(model=="ord"){ Z<-rZ_ord_fc(Z,EZ,rho,Y) }
    if(model=="cbin"){Z<-rZ_cbin_fc_reject(Z,EZ,rho,Y,odmax,odobs)}
    if(model=="rrl"){ Z<-rZ_rrl_fc(Z,EZ,rho,Y,YL)} 
    
    
    # draw random effects
    ab <- r_ab_DECORR_DEP(Z, XB, Sab,rho, s2, ULV)
    print(summary(c(Z))); print(summary(c(EZ)))
    a  <- ab$a * rvar
    b  <- ab$b * cvar
    
    if(symmetric == TRUE){
      a <- b <- (a + b)/2
    }
    EZ <- XB + ULV + outer(a, b, "+")
    
    # draw Sab and a -> z
    if(symmetric == FALSE){
      Sab <-rSab_fc(a,b,Sab0=prior$Sab0,eta0=prior$eta0)
      if(model=="frn"){
        tmp<-raSab_frn_fc_new(Z,Y,YL,a,b,Sab,odmax, odobs,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
        
      }
      if(model=="bin"){
        tmp<-raSab_bin_fc_new(Z,Y,a,b,Sab,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      if(model=="cbin")
      {
        tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs,
                           Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      
    } else{
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]
    }
    
    
    
    # posterior pred
    if(symmetric){ EZ<-(EZ+t(EZ))/2 } 
    
    
    
    if(model=="frn") { Ys_draw<-simY_frn(EZ, rho, odmax, YO = Y) }
    if(model=="bin") { Ys_draw<-simY_bin(EZ,rho) }
    if(model=="cbin"){ Ys_draw<-1*(simY_frn(EZ,rho,odmax,YO=Y)>0) }
    if(model=="rrl") { Ys<-simY_rrl(EZ,rho,odobs,YO=Y ) }
    if(model=="nrm") { Ys<-simY_nrm(EZ,rho,s2) }
    if(model=="ord") { Ys<-simY_ord(EZ,rho,Y) } 
    zsim <- simZ(EZ, rho)
    Ysim <- 1*(zsim > 0)
    degreeCounter <- rowSums(1*(zsim > 0))
    
    # Store results
    if(i == burnin + 1){ j = 0}
    if(i > burnin){
      
      if(i %% odens==0){
        j = j + 1
        
        for(g in 1:p){
          nam <-sprintf("beta_rMat%s[%s,]", g, j)
          nam2 <-sprintf("beta_cMat%s[%s,]", g, j)
          
          eval(parse(text= paste(nam, "<- beta_r[g,]")))        
          eval(parse(text= paste(nam2, "<- beta_c[g,]")))        
          
          
          
          
          
          
        }
        if(!is.null(X_d)){
          
          for(k in 1:dim(X_d)[3]){
            nam3 <-sprintf("beta_rMatDYAD%s[%s,]", k, j)
            nam4 <-sprintf("beta_cMatDYAD%s[%s,]", k, j)
            
            eval(parse(text= paste(nam3, "<- beta_dr[k,]")))        
            eval(parse(text= paste(nam4, "<- beta_dc[k,]")))        
            
          }
        }
        lambda_mat[j,] <- c(lambda)
        zmean[j,] <- mean(c(Z))
        Ys[j,] <-sum(abs(Ys_draw-Y),na.rm = TRUE)
        
        if(keep.UV == TRUE){
          UV_mat[j, ] <- c(U)
        }
        Z_mat[j, ] <- degreeCounter
      }
      
    }
    
  }
  
  ## make list to return
  
  results <- list()
  
  h = 1
  k=1
  while(h < 2*dim(X_r)[3]){
    
    nam <- paste("beta_rMat", k, sep = "")
    nam2 <- paste("beta_cMat", k, sep = "")
    
    results[[h]]   <- eval(parse(text= nam))
    names(results)[h] <- nam
    
    results[[h+1]] <- eval(parse(text= nam2))
    names(results)[h+1] <- nam2
    
    h <- h + 2
    k <-k+1
  }
  
  
  if(!is.null(X_d)){
    h = 1
    k=1
    while(h < 2*dim(X_d)[3]){
      
      nam3 <- paste("beta_rMatDYAD", k, sep = "")
      nam4 <- paste("beta_cMatDYAD", k, sep = "")
      
      results[[(2*dim(X_r)[3] + h)]]   <- eval(parse(text= nam3))
      names(results)[(2*dim(X_r)[3] + h)] <- nam3
      
      results[[(2*dim(X_r)[3] + h+1)]] <- eval(parse(text= nam4))
      names(results)[(2*dim(X_r)[3] + h +1)] <- nam4
      
      h <- h + 2
      k <-k+1
    }
    image(U%*%t(U))
    
    
  }
  if(keep.UV == TRUE){
    results[[length(results)+1]] <- UV_mat
    names(results)[[length(results)]] <- "U"
  }
  results[[length(results)+1]] <-   Z_mat 
  names(results)[[length(results)]] <- "degcount"
  
  
  print("COUNTER")
  print(counter/MH)
  print(counter)
  return(list(results=results,counter=counter, MH=MH))
  
  
}


# all given covariates are allowed to have a dependency structure
# accepts binary and fixed rank likelihoods ("frn", "bin")

## UandV (if true, then we draw membership for soc. and pop.)
## symmetricLambda (make Lambda symmetric or not)

# all given covariates are allowed to have a dependency structure
# accepts binary and fixed rank likelihoods ("frn", "bin")

## UandV (if true, then we draw membership for soc. and pop.)
## symmetricLambda (make Lambda symmetric or not)

subAmen_allDepBeta_giveTRUTH <-function(X_r, X_c, X_d, Y,iter,numGroup, prior_beta_mu, prior_beta_var,
                                        prior_alpha_mu, prior_alpha_var, start_beta_c,
                                        start_beta_r,start_beta_dyad, start_alpha, 
                                        UandV = FALSE, symmetricLambda = TRUE,  keep.UV = FALSE,
                                        dcor = FALSE, symmetric = TRUE, model = "frn",odmax = NULL,odens = 10, 
                                        burnin = 1000,badInit = FALSE,prior_group_probs = NULL){
  
  
  
  
  
  p <- dim(X_r)[3]
  n <- dim(X_r)[1]
  ###################### sampling
  beta_r <- start_beta_r
  beta_c <- start_beta_c
  
  if(!is.null(X_d)){
    beta_dc <- beta_dr <- start_beta_dyad
    
  }else{
    beta_dc <- beta_dr <- 0
    
  }
  
  
  if(symmetric == FALSE){
    rvar <- cvar <-  TRUE
    nvar <- FALSE
  }else{
    rvar <-cvar <- TRUE
    nvar <- TRUE
  }
  counter <- 0
  MH <-0
  
  
  Y_binary <- 1 * (Y > 0)
  
  
  if (is.element(model, c("cbin", "frn", "rrl"))) {
    odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
    if (length(odmax) == 1) {
      odmax <- rep(odmax, nrow(Y))
    }
  }
  
  
  
  diag(Y) <- NA
  
  
  
  if (is.element(model, c("frn", "rrl"))) {
    
    ymx <- max(apply(1 * (Y > 0), 1, sum, na.rm = TRUE))
    YL <- NULL
    warn <- FALSE
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(yi[rnkd]) > length(unique(yi[rnkd]))) {
        warn <- TRUE
      }
      yi[rnkd] <- rank(yi[rnkd], ties.method = "random")
      Y[i, ] <- yi
      YL <- rbind(YL, match(1:ymx, yi))
    }
    if (warn) {
      cat("WARNING: Random reordering used to break ties in ranks\n")
    }
  }
  
  if (model == "bin") {
    Y <- 1 * (Y > 0)
    
    Z <- matrix(zscores(Y), nrow(Y), nrow(Y))
    z01 <- 0.5 * (max(Z[Y == 0], na.rm = TRUE) + min(Z[Y == 
                                                         1], na.rm = TRUE))
    Z <- Z - z01
  }
  if (is.element(model, c("cbin", "frn"))) {
    Z <- Y
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      zi <- zscores(yi)
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
        zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 0.001
      }
      if (length(rnkd) < odmax[i]) {
        urnkd <- which(!is.na(yi) & yi == 0)
        if (max(zi[urnkd]) > 0) {
          zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 0.001
        }
      }
      Z[i, ] <- zi
    }
  }
  
  diag(Y) <- 0
  
  ######
  ## initialize U using regularizer from Rohe
  diag(Y_binary) <- 0
  UV_start <- specClust(Y_binary, tauAvg = TRUE,  numGroup,UandV)
  # U <- V <- matrix(UV_start$U_hat, ncol = numGroup)
  U <- V <- matrix(c(rep(c(1,0,0), n/numGroup), rep(c(0,1,0), n/numGroup), 
                     rep(c(0,0,1), n/numGroup)),
                   ncol = numGroup,nrow = n, byrow = TRUE)
  #image(U %*%t(U))
  #V <-  matrix(UV_start$V_hat, ncol = numGroup)
  #  starting point
  
  if(badInit == TRUE){
    set.seed(13)
    k<- dim(U)[2]
    for(i in 1:n){
      membership <- sample(1:k, 1)
      U[i,membership] <- 1
      U[i, seq(1,k)[-membership]] <- 0
    }
    V <-U
    image(U%*%t(U))
  }
  
  ## initialize lambda
  prior.lambda.mu  <- rep(1, numGroup^2)
  prior.lambda.var <- prior.lambda.prec <- diag(1, numGroup^2) ## have had .2
  lambda <- matrix(prior.lambda.mu, nrow = numGroup)
  
  ## alpha
  prior_alpha     <-  prior_alpha_mu
  prior_var_alpha <- prior_alpha_var
  alpha           <- start_alpha
  
  # get initial XBeta
  h_r_past <- NULL
  h_c_past <- NULL
  XB_ind_get <- XBETA_DEP(X_r, X_c, X_d, beta_r, 
                          beta_c, beta_dr, beta_dc, 
                          U, V,TRUE, h_r_past, h_c_past) 
  
  XB_ind     <- XB_ind_get$XBETA
  
  if(!is.null(X_d)){
    XB_dyad    <- XB_ind_get$XBETADYAD
    XB         <- rowSums(XB_dyad, dims = 2) + rowSums(XB_ind, dims =2 )+ alpha # calculate XB
    
  }else{
    XB         <-  rowSums(XB_ind, dims =2 )+ alpha # calculate XB
    
  }
  
  E <- Z- XB
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  prior <- list()
  if(model!="nrm")
  {  
    ydist<-table(Y)
    ymode<-as.numeric(names(ydist)[ ydist==max(ydist) ])[1] 
    ## eg, in a sparse binary network, ymode will be zero 
    YB<-1*(Y!=ymode) 
    ybar<-mean(YB,na.rm=TRUE) ; mu<-qnorm(ybar)
    E<- (YB - ybar)/dnorm(qnorm(ybar)) ; diag(E)<-0
    a<-rowMeans(E,na.rm=TRUE)  ; a[is.na(a)]<-0 
    b<-colMeans(E,na.rm=TRUE)  ; b[is.na(b)]<-0
    vscale<-mean(diag(cov(cbind(a,b))))
    PHAT<-pnorm(mu+outer(a,b,"+"))
    vdfmlt<-.25/mean(PHAT*(1-PHAT))
    if(is.null(prior$Sab0)){ prior$Sab0<-diag(2)*vscale }
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$kappa0)){ prior$kappa0<-round((2*numGroup+2)*vdfmlt) }  
    if(is.null(prior$g)){ prior$g<-sum(!is.na(Y)) }
  }
  
  Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar)) 
  E <- E- outer(a,b, "+")
  s2<-1
  if(model=="nrm"){s2<-mean(E^2)}
  rho<-cor( c(E[upper.tri(E)]), c(t(E)[upper.tri(E)]) )*dcor  
  
  ## initialize things to store
  for(i in 1: dim(X_r)[3]){
    nam <- paste("beta_rMat", i, sep = "")
    assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
    
    namc <- paste("beta_cMat", i, sep = "")
    assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
    
    
    
  }
  
  if(!is.null(X_d)){
    for(i in 1:dim(X_d)[3]){
      
      nam <- paste("beta_rMatDYAD", i, sep = "")
      assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
      
      namc <- paste("beta_cMatDYAD", i, sep = "")
      assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
      
    }
  }
  
  ulv <- U %*% lambda %*% t(V)
  lambda_mat <- matrix(nrow = iter/odens, ncol = numGroup^2)
  if(keep.UV == TRUE){
    UV_mat      <-  matrix(nrow = iter/odens, ncol = length(c(U)))
    
    
  }
  
  
  
  if(is.null(prior_group_probs)){
    prior_group_probs <- rep(1/numGroup, numGroup)
  }
  
  Z_mat <- matrix(nrow = iter/odens, ncol = n)
  
  
  zmean <- matrix(nrow = iter/odens, ncol = 1)
  Ys <-matrix(nrow = iter/odens, ncol = 1)
  
  updatedUV <- TRUE
  XB_ind_get    <- XBETA_DEP_partition(X_r, X_c, X_d, beta_r, beta_c,beta_dr, beta_dc, U,V,updatedUV, h_r_past, h_c_past)
  XBETA_ROW     <- XB_ind_get$XBETAROW
  XBETA_COL     <- XB_ind_get$XBETACOL
  
  if(!is.null(X_d)){
    
    XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
  }else{
    XBETA_DYAD <- 0
  }
  XB_ind        <- XBETA_ROW + XBETA_COL 
  h_r_past      <- XB_ind_get$h_r_past
  h_c_past      <- XB_ind_get$h_c_past
  
  
  XB            <- rowSums(XB_ind, dims =2 )+ XBETA_DYAD + alpha # calculate XB
  
  EZ <- XB + ulv  + outer(a, b, "+")
  updatedUV <- TRUE
  diag(Y_binary) <-0
  #diag(Z) <- -1
  ## start sampling
  diag(Z) <- mean(na.omit(c(Z)))
  Se   <-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2 <-mhalf(solve(Se))
  td   <-iSe2[1,1] ; to<-iSe2[1,2]
  #if(model == "bin"){Z  <- rZ_bin_fc_new(Z, EZ, rho, Y) }
  
  # j = 0
  for(i in 1:(iter+burnin)){
    
    
    # Draw latent lambda and U
    lambda <- drawLambda(Z, XB + outer(a,b, "+"), U,V, prior.lambda.var,
                         prior.lambda.mu, td, to,symmetricLambda = symmetricLambda)
    
    # 
    # UV_info      <- drawUV_normalLik_OLD(lambda, U, V,  XB ,a,b,
    #                                      Y_binary, UandV,counter, MH,
    #                                      X_r, X_c,X_d, beta_r, beta_c,alpha,
    #                                      beta_dr, beta_dc,
    #                                      h_r_past, h_c_past,Z,td,to,Sab, rho)#, prior_group_probs)
    # 
    # UV_info <- drawUV_normalLik_alpha(lambda, U, V,  XB ,a,b,
    #                             Y_binary, UandV,counter, MH,
    #                             X_r, X_c,X_d, beta_r, beta_c,alpha,
    #                             beta_dr, beta_dc,
    #                             h_r_past, h_c_past,Z,td,to,Sab, rho)
    # 
    #  U            <- UV_info$U
    #  V            <- UV_info$V
    
    # print(UandV)
    #  print(colSums(U))
    #  print(colSums(V))
    
    
    # updatedUV     <- UV_info$updateU
    #swapPeep     <- UV_info$swapPeep
    #counter      <- UV_info$counter
    #MH           <- UV_info$MH
    
    ULV <- U %*%lambda%*%t(V)
    
    print(rho)
    
    XB_ind_get    <- XBETA_DEP_partition(X_r, X_c,X_d, beta_r, beta_c,
                                         beta_dr, beta_dc, U,V,
                                         updatedUV, h_r_past, h_c_past)
    XBETA_ROW     <- XB_ind_get$XBETAROW
    XBETA_COL     <- XB_ind_get$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
    }else{
      XBETA_DYAD <- 0
    }
    
    XB_ind        <- XBETA_ROW + XBETA_COL 
    h_r_past      <- XB_ind_get$h_r_past
    h_c_past      <- XB_ind_get$h_c_past
    XB            <- rowSums(XB_ind, dims =2 )+ alpha +XBETA_DYAD # calculate XB
    #
    EZ <- XB + ULV + outer(a, b, "+")
    
    if(model=="nrm"){ s2<-rs2_fc(Z-EZ,rho,nu0=prior$nu0,s20=prior$s20)  } 
    
    # Draw rho
    
    if(symmetric){rho <- min(0.9999, 1 - 1/sqrt(i))}
    if(dcor){ rho <- rrho_mh_new(Z-EZ,rho,s2,offset=0,asp = NULL)} 
    ## for decorrelation
    Se   <-matrix(c(1,rho,rho,1),2,2)*s2
    iSe2 <-mhalf(solve(Se))
    td   <-iSe2[1,1] ; to<-iSe2[1,2]
    
    
    ## draw intercept
    XB_tmp <- XB -alpha
    
    alpha_tmp <-  rALPHA(Z, ULV,
                         alpha, a , b, XB,td, to,prior_alpha_mu,prior_alpha_var)
    alpha <- alpha_tmp$alpha 
    
    XB_updatedALPHA <- XB_tmp + alpha
    
    if(!is.null(X_d)){
      
      
      
      numEigs <- 1 ### CHANGE
      
      tmpDYAD <- rbeta_DECORR_DEPDYAD(Z-XB_updatedALPHA+XBETA_DYAD, X_d, U,V, lambda, 
                                      beta_dr, beta_dc, 
                                      prior_beta_mu, 
                                      prior_beta_var , a, b,   to ,td,numEigs, h_r_past,
                                      h_c_past)
      
      
      beta_dr <- tmpDYAD$beta_row_dyad
      beta_dc <- tmpDYAD$beta_col_dyad
      print("DYADF")
      print(tmpDYAD)
      
      XBtmp <- XB_updatedALPHA -XBETA_DYAD
      XB_updateDyad_get    <- XBETA_DEP_UPDATEDDYAD(X_d, 
                                                    beta_dr, beta_dc,
                                                    U,V, 
                                                    FALSE, 
                                                    h_r_past, 
                                                    h_c_past)
      XBETA_DYAD    <- rowSums(XB_updateDyad_get$XBETADYAD, dims =2 )
      XB            <-  XBtmp +XBETA_DYAD # calculate XB
      EZ <- XB + ULV + outer(a, b, "+")
      
      
    }else{
      XB <- XB_updatedALPHA
      EZ <- XB + ULV + outer(a, b, "+")
      
    }
    #
    
    
    
    # Draw Beta
    
    
    tmp <- rbeta_DECORR_DEP(Z, X_r, X_c,XBETA_DYAD,U,V, lambda, dim(X_r)[3], beta_r, beta_c,
                            alpha, prior_beta_mu, 
                            prior_beta_var, a, b,  h_r_past, h_c_past, to ,td,updatedUV )
    
    print("ITER")
    print(i)
    
    beta_r <- tmp$beta_row
    beta_c <- tmp$beta_col
    
    
    XBtmp <- XB - rowSums(XB_ind, dims =2 )
    
    # recalculate XB
    
    
    XB_updateRC <- XBETA_DEP_UPDATEROWCOLUMN(X_r, X_c, beta_r, beta_c, 
                                             U,V,
                                             FALSE, 
                                             h_r_past, 
                                             h_c_past)
    
    XBETA_ROW     <-  XB_updateRC$XBETAROW
    XBETA_COL     <-  XB_updateRC$XBETACOL
    XB_ind        <-  XBETA_ROW + XBETA_COL 
    XB            <-  XBtmp + rowSums(XB_ind, dims =2)
    #
    EZ <- XB + ULV + outer(a, b, "+")
    
    
    
    
    
    
    
    # draw Z
    if(model == "frn"){Z <- rZ_frn_fc_new(Z, EZ, rho, Y, YL, odmax, odobs)}
    if(model == "bin"){Z  <- rZ_bin_fc_new(Z, EZ, rho, Y) }
    
    if(model=="nrm"){ Z<-rZ_nrm_fc(Z,EZ,rho,s2,Y) } 
    if(model=="ord"){ Z<-rZ_ord_fc(Z,EZ,rho,Y) }
    if(model=="cbin"){Z<-rZ_cbin_fc_reject(Z,EZ,rho,Y,odmax,odobs)}
    if(model=="rrl"){ Z<-rZ_rrl_fc(Z,EZ,rho,Y,YL)} 
    
    
    # draw random effects
    ab <- r_ab_DECORR_DEP(Z, XB, Sab,rho, s2, ULV)
    print(summary(c(Z))); print(summary(c(EZ)))
    a  <- ab$a * rvar
    b  <- ab$b * cvar
    
    if(symmetric == TRUE){
      a <- b <- (a + b)/2
    }
    EZ <- XB + ULV + outer(a, b, "+")
    
    # draw Sab and a -> z
    if(symmetric == FALSE){
      Sab <-rSab_fc(a,b,Sab0=prior$Sab0,eta0=prior$eta0)
      if(model=="frn"){
        tmp<-raSab_frn_fc_new(Z,Y,YL,a,b,Sab,odmax, odobs,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
        
      }
      if(model=="bin"){
        tmp<-raSab_bin_fc_new(Z,Y,a,b,Sab,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      if(model=="cbin")
      {
        tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs,
                           Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      
    } else{
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]
    }
    
    
    
    # posterior pred
    if(symmetric){ EZ<-(EZ+t(EZ))/2 } 
    
    
    
    if(model=="frn") { Ys_draw<-simY_frn(EZ, rho, odmax, YO = Y) }
    if(model=="bin") { Ys_draw<-simY_bin(EZ,rho) }
    if(model=="cbin"){ Ys_draw<-1*(simY_frn(EZ,rho,odmax,YO=Y)>0) }
    if(model=="rrl") { Ys<-simY_rrl(EZ,rho,odobs,YO=Y ) }
    if(model=="nrm") { Ys<-simY_nrm(EZ,rho,s2) }
    if(model=="ord") { Ys<-simY_ord(EZ,rho,Y) } 
    zsim <- simZ(EZ, rho)
    Ysim <- 1*(zsim > 0)
    degreeCounter <- rowSums(1*(zsim > 0))
    
    # Store results
    if(i == burnin + 1){ j = 0}
    if(i > burnin){
      
      if(i %% odens==0){
        j = j + 1
        
        for(g in 1:p){
          nam <-sprintf("beta_rMat%s[%s,]", g, j)
          nam2 <-sprintf("beta_cMat%s[%s,]", g, j)
          
          eval(parse(text= paste(nam, "<- beta_r[g,]")))        
          eval(parse(text= paste(nam2, "<- beta_c[g,]")))        
          
          
          
          
          
          
        }
        if(!is.null(X_d)){
          
          for(k in 1:dim(X_d)[3]){
            nam3 <-sprintf("beta_rMatDYAD%s[%s,]", k, j)
            nam4 <-sprintf("beta_cMatDYAD%s[%s,]", k, j)
            
            eval(parse(text= paste(nam3, "<- beta_dr[k,]")))        
            eval(parse(text= paste(nam4, "<- beta_dc[k,]")))        
            
          }
        }
        lambda_mat[j,] <- c(lambda)
        zmean[j,] <- mean(c(Z))
        Ys[j,] <-sum(abs(Ys_draw-Y),na.rm = TRUE)
        
        if(keep.UV == TRUE){
          UV_mat[j, ] <- c(U)
        }
        Z_mat[j, ] <- degreeCounter
      }
      
    }
    
  }
  
  ## make list to return
  
  results <- list()
  
  h = 1
  k=1
  while(h < 2*dim(X_r)[3]){
    
    nam <- paste("beta_rMat", k, sep = "")
    nam2 <- paste("beta_cMat", k, sep = "")
    
    results[[h]]   <- eval(parse(text= nam))
    names(results)[h] <- nam
    
    results[[h+1]] <- eval(parse(text= nam2))
    names(results)[h+1] <- nam2
    
    h <- h + 2
    k <-k+1
  }
  
  
  if(!is.null(X_d)){
    h = 1
    k=1
    while(h < 2*dim(X_d)[3]){
      
      nam3 <- paste("beta_rMatDYAD", k, sep = "")
      nam4 <- paste("beta_cMatDYAD", k, sep = "")
      
      results[[(2*dim(X_r)[3] + h)]]   <- eval(parse(text= nam3))
      names(results)[(2*dim(X_r)[3] + h)] <- nam3
      
      results[[(2*dim(X_r)[3] + h+1)]] <- eval(parse(text= nam4))
      names(results)[(2*dim(X_r)[3] + h +1)] <- nam4
      
      h <- h + 2
      k <-k+1
    }
    image(U%*%t(U))
    
    
  }
  if(keep.UV == TRUE){
    results[[length(results)+1]] <- UV_mat
    names(results)[[length(results)]] <- "U"
  }
  results[[length(results)+1]] <-   Z_mat 
  names(results)[[length(results)]] <- "degcount"
  
  
  print("COUNTER")
  print(counter/MH)
  print(counter)
  return(list(results=results,counter=counter, MH=MH))
  
  
}


# Main function to run MCMC, all covariates are assumed to be dep. on community
#' @param X_r: n x n x p_dep array of row covar
#' @param X_c: n x n x p_dep array of col covar
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param Y: Observed network (can be binary or censored binary)
#' @parm iter: number of iterations to run mcmc for 
#' @param numGroup: number of communities (also denote k)
#' @param prior_beta_mu: k dimensional vector with prior mean for beta
#' @param prior_beta_var:k x k matrix for prior variance for beta
#' @param prior_alpha_mu: constant for alpha(intercept) mean
#' @param prior_alpha_var:constant for alpha(intercept) var
#' @param start_beta_c: p_dep x k matrix for initializing beta_c
#' @param start_beta_r:  p_dep x k matrix for initializing beta_c
#' @param start_beta_dyad:  p_dep_dyad x k matrix for initializing beta_c
#' @param start_alpha: constant to initialize alpha
#' @param UandV: Default is FALSE, indicates if U and V should be calculated or just U
#' @param symmetricLambda: Default is TRUE, indicates if lambda should be symmetric
#' @param keep.UV: Default is TRUE, store U/V?
#' @param dcor: Default is TRUE, indicates if correlation in error terms should be estimated
#' @param symmetric: Default is FALSE, indicates that Y is not symmetric
#' @param model:  "cbin" or "bin", indicates what type of network Y is 
#' @param odmax: if model = "cbin", this is the max out degree a node can have
#' @param odens = 10: how often to save output
#' @param burnin: amount of burnin iterations
#' @param badInit: start with completely random initialization for U/V
#' @param prior_group_probs: For future use (if we want to change prior prob on group membership)


subAmen_allDepBeta <-function(X_r, X_c, X_d, Y,iter,numGroup, prior_beta_mu,
                              prior_beta_var,
                              prior_alpha_mu, prior_alpha_var, start_beta_c,
                              start_beta_r,start_beta_dyad, start_alpha, 
                              UandV = FALSE, symmetricLambda = TRUE,  keep.UV = TRUE,
                              dcor = TRUE, symmetric = FALSE, model = "bin",
                              odmax = NULL,odens = 10, 
                              burnin = 1000,badInit = FALSE,prior_group_probs = NULL){
  
  
  
  
  
  p <- dim(X_r)[3]
  n <- dim(X_r)[1]
  ###################### sampling
  beta_r <- start_beta_r
  beta_c <- start_beta_c
  
  if(!is.null(X_d)){
    beta_dc <- beta_dr <- start_beta_dyad
    
  }else{
    beta_dc <- beta_dr <- 0
    
  }
  
  
  if(symmetric == FALSE){
    rvar <- cvar <-  TRUE
    nvar <- FALSE
  }else{
    rvar <-cvar <- TRUE
    nvar <- TRUE
  }
  counter <- 0
  MH <-0
  
  
  Y_binary <- 1 * (Y > 0)
  
  
  if (is.element(model, c("cbin"))) {
    odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
    if (length(odmax) == 1) {
      odmax <- rep(odmax, nrow(Y))
    }
  }
  
  
  diag(Y) <- NA
  
  if (model == "bin") {
    Y <- 1 * (Y > 0)
    
    Z <- matrix(zscores(Y), nrow(Y), nrow(Y))
    z01 <- 0.5 * (max(Z[Y == 0], na.rm = TRUE) + min(Z[Y == 
                                                         1], na.rm = TRUE))
    Z <- Z - z01
  }
  if (is.element(model, c("cbin"))) {
    Z <- Y
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      zi <- zscores(yi)
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
        zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 0.001
      }
      if (length(rnkd) < odmax[i]) {
        urnkd <- which(!is.na(yi) & yi == 0)
        if (max(zi[urnkd]) > 0) {
          zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 0.001
        }
      }
      Z[i, ] <- zi
    }
  }
  
  diag(Y) <- 0
  
  ######
  ## initialize U using regularizer from Rohe
  diag(Y_binary) <- 0
  UV_start <- specClust(Y_binary, tauAvg = TRUE,  numGroup,UandV)
  U <- V <- matrix(UV_start$U_hat, ncol = numGroup)
  
  
  if(badInit == TRUE){
    set.seed(13)
    k<- dim(U)[2]
    for(i in 1:n){
      membership <- sample(1:k, 1)
      U[i,membership] <- 1
      U[i, seq(1,k)[-membership]] <- 0
    }
    V <-U
    image(U%*%t(U))
  }
  
  ## initialize lambda
  prior.lambda.mu  <- rep(1, numGroup^2)
  prior.lambda.var <- prior.lambda.prec <- diag(1, numGroup^2) ## have had .2
  lambda <- matrix(prior.lambda.mu, nrow = numGroup)
  
  ## alpha
  prior_alpha     <-  prior_alpha_mu
  prior_var_alpha <- prior_alpha_var
  alpha           <- start_alpha
  
  # get initial XBeta
  h_r_past <- NULL
  h_c_past <- NULL
  XB_ind_get <- XBETA_DEP(X_r, X_c, X_d, beta_r, 
                          beta_c, beta_dr, beta_dc, 
                          U, V,TRUE, h_r_past, h_c_past) 
  
  XB_ind     <- XB_ind_get$XBETA
  
  if(!is.null(X_d)){
    XB_dyad    <- XB_ind_get$XBETADYAD
    XB         <- rowSums(XB_dyad, dims = 2) + rowSums(XB_ind, dims =2 )+ alpha # calculate XB
    
  }else{
    XB         <-  rowSums(XB_ind, dims =2 )+ alpha # calculate XB
    
  }
  
  E <- Z- XB - U%*%lambda%*%t(V)
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  prior <- list()
  if(model!="nrm")
  {  
    ydist<-table(Y)
    ymode<-as.numeric(names(ydist)[ ydist==max(ydist) ])[1] 
    ## eg, in a sparse binary network, ymode will be zero 
    YB<-1*(Y!=ymode) 
    ybar<-mean(YB,na.rm=TRUE) ; mu<-qnorm(ybar)
    E<- (YB - ybar)/dnorm(qnorm(ybar)) ; diag(E)<-0
    a<-rowMeans(E,na.rm=TRUE)  ; a[is.na(a)]<-0 
    b<-colMeans(E,na.rm=TRUE)  ; b[is.na(b)]<-0
    vscale<-mean(diag(cov(cbind(a,b))))
    PHAT<-pnorm(mu+outer(a,b,"+"))
    vdfmlt<-.25/mean(PHAT*(1-PHAT))
    if(is.null(prior$Sab0)){ prior$Sab0<-diag(2)*vscale }
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$kappa0)){ prior$kappa0<-round((2*numGroup+2)*vdfmlt) }  
    if(is.null(prior$g)){ prior$g<-sum(!is.na(Y)) }
  }
  
  Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar)) 
  E <- E- outer(a,b, "+")
  s2<-1
  rho<-cor( c(E[upper.tri(E)]), c(t(E)[upper.tri(E)]) )*dcor  
  
  ## initialize things to store
  for(i in 1: dim(X_r)[3]){
    nam <- paste("beta_rMat", i, sep = "")
    assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
    
    namc <- paste("beta_cMat", i, sep = "")
    assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
    
    
    
  }
  
  if(!is.null(X_d)){
    for(i in 1:dim(X_d)[3]){
      
      nam <- paste("beta_rMatDYAD", i, sep = "")
      assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
      
      namc <- paste("beta_cMatDYAD", i, sep = "")
      assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
      
    }
  }
  
  ulv<- U %*% lambda %*% t(V)
  lambda_mat <- matrix(nrow = iter/odens, ncol = numGroup^2)
  if(keep.UV == TRUE){
    UV_mat      <-  matrix(nrow = iter/odens, ncol = length(c(U)))
    
    
  }
  
  
  
  if(is.null(prior_group_probs)){
    prior_group_probs <- rep(1/numGroup, numGroup)
  }
  
  Z_mat <- matrix(nrow = iter/odens, ncol = n)
  
  
  zmean <- matrix(nrow = iter/odens, ncol = 1)
  Ys <-matrix(nrow = iter/odens, ncol = 1)
  
  updatedUV <- TRUE
  XB_ind_get    <- XBETA_DEP_partition(X_r, X_c, X_d, beta_r, beta_c,beta_dr, beta_dc, U,V,updatedUV, h_r_past, h_c_past)
  XBETA_ROW     <- XB_ind_get$XBETAROW
  XBETA_COL     <- XB_ind_get$XBETACOL
  
  if(!is.null(X_d)){
    
    XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
  }else{
    XBETA_DYAD <- 0
  }
  XB_ind        <- XBETA_ROW + XBETA_COL 
  h_r_past      <- XB_ind_get$h_r_past
  h_c_past      <- XB_ind_get$h_c_past
  
  
  XB            <- rowSums(XB_ind, dims =2 )+ XBETA_DYAD + alpha # calculate XB
  
  EZ <- XB + ulv  + outer(a, b, "+")
  updatedUV <- TRUE
  diag(Y_binary) <-0
  ## start sampling
  diag(Z) <- mean(na.omit(c(Z)))
  Se   <-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2 <-mhalf(solve(Se))
  td   <-iSe2[1,1] ; to<-iSe2[1,2]
  
  
  for(i in 1:(iter+burnin)){
    print(i)
    
    # Draw latent lambda and U
    lambda <- drawLambda(Z, XB + outer(a,b, "+"), U,V, prior.lambda.prec,
                         prior.lambda.mu,td,to, symmetricLambda = symmetricLambda)
    
    
    UV_info      <- drawUV_probitlLikelihood(lambda, U, V,  XB ,a,b,
                                             Y_binary, UandV,counter, MH,
                                             X_r, X_c,X_d, beta_r, beta_c,alpha,
                                             beta_dr, beta_dc,
                                             h_r_past, h_c_past,Z,td,to,Sab, rho)#, prior_group_probs)
    
    
    U            <- UV_info$U
    V            <- UV_info$V
    
    updatedUV     <- UV_info$updateU
    swapPeep     <- UV_info$swapPeep
    counter      <- UV_info$counter
    MH           <- UV_info$MH
    
    
    
    
    
    
    
    ULV <- U %*%lambda%*%t(V)
    
    
    XB_ind_get    <- XBETA_DEP_partition(X_r, X_c,X_d, beta_r, beta_c,
                                         beta_dr, beta_dc, U,V,
                                         updatedUV, h_r_past, h_c_past)
    XBETA_ROW     <- XB_ind_get$XBETAROW
    XBETA_COL     <- XB_ind_get$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
    }else{
      XBETA_DYAD <- 0
    }
    
    XB_ind        <- XBETA_ROW + XBETA_COL 
    h_r_past      <- XB_ind_get$h_r_past
    h_c_past      <- XB_ind_get$h_c_past
    XB            <- rowSums(XB_ind, dims =2 )+ alpha +XBETA_DYAD # calculate XB
    #
    EZ <- XB + ULV + outer(a, b, "+")
    
    if(updatedUV == TRUE){
      for(h in 1:length(swapPeep)){
        Z <-     update_Z_bin(swapPeep[h], Y, EZ, rho,  n,Z)
      }
      
    }
    
    
    
    
    # Draw rho
    
    if(symmetric){rho <- min(0.9999, 1 - 1/sqrt(i))}
    if(dcor){ rho <- rrho_mh_new(Z-EZ,rho,s2,offset=0,asp = NULL)} 
    ## for decorrelation
    Se   <-matrix(c(1,rho,rho,1),2,2)*s2
    iSe2 <-mhalf(solve(Se))
    td   <-iSe2[1,1] ; to<-iSe2[1,2]
    
    
    ## draw intercept
    XB_tmp <- XB -alpha
    
    alpha_tmp <-  rALPHA(Z, ULV,
                         alpha, a , b, XB,td, to,prior_alpha_mu,prior_alpha_var)
    alpha <- alpha_tmp$alpha 
    
    XB_updatedALPHA <- XB_tmp + alpha
    
    if(!is.null(X_d)){
      
      
      
      numEigs <- 1 ## might want to update this, for comp. eff.
      
      tmpDYAD <- rbeta_DECORR_DEPDYAD(Z-XB_updatedALPHA+XBETA_DYAD, X_d, U,V, lambda, 
                                      beta_dr, beta_dc, 
                                      prior_beta_mu, 
                                      prior_beta_var , a, b,   to ,td,numEigs, h_r_past,
                                      h_c_past)
      
      
      beta_dr <- tmpDYAD$beta_row_dyad
      beta_dc <- tmpDYAD$beta_col_dyad
      
      
      XBtmp <- XB_updatedALPHA -XBETA_DYAD
      XB_updateDyad_get    <- XBETA_DEP_UPDATEDDYAD(X_d, 
                                                    beta_dr, beta_dc,
                                                    U,V)
      XBETA_DYAD    <- rowSums(XB_updateDyad_get$XBETADYAD, dims =2 )
      XB            <-  XBtmp +XBETA_DYAD # calculate XB
      EZ <- XB + ULV + outer(a, b, "+")
      
      
    }else{
      XB <- XB_updatedALPHA
      EZ <- XB + ULV + outer(a, b, "+")
      
    }
    #
    
    
    
    # Draw Beta
    
    
    tmp <- rbeta_DECORR_DEP(Z, X_r, X_c,XBETA_DYAD,U,V, lambda, dim(X_r)[3], beta_r, beta_c,
                            alpha, prior_beta_mu, 
                            prior_beta_var, a, b,  h_r_past, h_c_past, to ,td,updatedUV )
    
    
    
    beta_r <- tmp$beta_row
    beta_c <- tmp$beta_col
    
    
    XBtmp <- XB - rowSums(XB_ind, dims =2 )
    
    # recalculate XB
    
    
    XB_updateRC <- XBETA_DEP_UPDATEROWCOLUMN(X_r, X_c, beta_r, beta_c, 
                                             U,V,
                                             FALSE, 
                                             h_r_past, 
                                             h_c_past)
    
    XBETA_ROW     <-  XB_updateRC$XBETAROW
    XBETA_COL     <-  XB_updateRC$XBETACOL
    XB_ind        <-  XBETA_ROW + XBETA_COL 
    XB            <-  XBtmp + rowSums(XB_ind, dims =2)
    #
    EZ <- XB + ULV + outer(a, b, "+")
    
    
    
    
    
    
    
    # draw Z
    if(model == "bin"){Z  <- rZ_bin_fc_new(Z, EZ, rho, Y) }
    if(model=="cbin"){Z<-rZ_cbin_fc_reject(Z,EZ,rho,Y,odmax,odobs)}
    
    
    # draw random effects
    ab <- r_ab_DECORR_DEP(Z, XB, Sab,rho, s2, ULV)
    a  <- ab$a * rvar
    b  <- ab$b * cvar
    
    if(symmetric == TRUE){
      a <- b <- (a + b)/2
    }
    EZ <- XB + ULV + outer(a, b, "+")
    
    # draw Sab and a -> z
    if(symmetric == FALSE){
      Sab <-rSab_fc(a,b,Sab0=prior$Sab0,eta0=prior$eta0)
      if(model=="frn"){
        tmp<-raSab_frn_fc_new(Z,Y,YL,a,b,Sab,odmax, odobs,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
        
      }
      if(model=="bin"){
        tmp<-raSab_bin_fc_new(Z,Y,a,b,Sab,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      if(model=="cbin")
      {
        tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs,
                           Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      
    } else{
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]
    }
    
    
    
    # posterior pred
    if(symmetric){ EZ<-(EZ+t(EZ))/2 } 
    
    
    
    if(model=="bin") { Ys_draw<-simY_bin(EZ,rho) }
    if(model=="cbin"){ Ys_draw<-1*(simY_frn(EZ,rho,odmax,YO=Y)>0) }
    
    zsim <- simZ(EZ, rho)
    Ysim <- 1*(zsim > 0)
    degreeCounter <- rowSums(1*(zsim > 0))
    
    # Store results
    if(i == burnin + 1){ j = 0}
    if(i > burnin){
      
      if(i %% odens==0){
        j = j + 1
        
        for(g in 1:p){
          nam <-sprintf("beta_rMat%s[%s,]", g, j)
          nam2 <-sprintf("beta_cMat%s[%s,]", g, j)
          
          eval(parse(text= paste(nam, "<- beta_r[g,]")))        
          eval(parse(text= paste(nam2, "<- beta_c[g,]")))        
          
          
          
          
          
          
        }
        if(!is.null(X_d)){
          
          for(k in 1:dim(X_d)[3]){
            nam3 <-sprintf("beta_rMatDYAD%s[%s,]", k, j)
            nam4 <-sprintf("beta_cMatDYAD%s[%s,]", k, j)
            
            eval(parse(text= paste(nam3, "<- beta_dr[k,]")))        
            eval(parse(text= paste(nam4, "<- beta_dc[k,]")))        
            
          }
        }
        lambda_mat[j,] <- c(lambda)
        zmean[j,] <- mean(c(Z))
        Ys[j,] <-sum(abs(Ys_draw-Y),na.rm = TRUE)
        
        if(keep.UV == TRUE){
          UV_mat[j, ] <- c(U)
        }
        Z_mat[j, ] <- degreeCounter
      }
      
    }
    
    
    
    if(i %% 10000 ==0){
      results <- list()
      
      h = 1
      k=1
      while(h < 2*dim(X_r)[3]){
        
        nam <- paste("beta_rMat", k, sep = "")
        nam2 <- paste("beta_cMat", k, sep = "")
        
        results[[h]]   <- eval(parse(text= nam))
        names(results)[h] <- nam
        
        results[[h+1]] <- eval(parse(text= nam2))
        names(results)[h+1] <- nam2
        
        h <- h + 2
        k <-k+1
      }
      
      
      if(!is.null(X_d)){
        h = 1
        k=1
        while(h < 2*dim(X_d)[3]){
          
          nam3 <- paste("beta_rMatDYAD", k, sep = "")
          nam4 <- paste("beta_cMatDYAD", k, sep = "")
          
          results[[(2*dim(X_r)[3] + h)]]   <- eval(parse(text= nam3))
          names(results)[(2*dim(X_r)[3] + h)] <- nam3
          
          results[[(2*dim(X_r)[3] + h+1)]] <- eval(parse(text= nam4))
          names(results)[(2*dim(X_r)[3] + h +1)] <- nam4
          
          h <- h + 2
          k <-k+1
        }
        
        
      }
      if(keep.UV == TRUE){
        results[[length(results)+1]] <- UV_mat
        names(results)[[length(results)]] <- "U"
      }
      results[[length(results)+1]] <-   Z_mat 
      names(results)[[length(results)]] <- "degcount"
      
      results[[length(results)+1]] <-   lambda_mat
      names(results)[[length(results)]] <- "lambda"
      saveRDS(results, sprintf("/hpchome/volfovskylab/hmm40/CBIN/results_cbin_disprop_comm_censored_AMEN_333_10_%s_OCT%s%s.rds","givenName", n,i))
    }
    
  }
  
  
  
  
  ## make list to return
  
  results <- list()
  
  h = 1
  k=1
  while(h < 2*dim(X_r)[3]){
    
    nam <- paste("beta_rMat", k, sep = "")
    nam2 <- paste("beta_cMat", k, sep = "")
    
    results[[h]]   <- eval(parse(text= nam))
    names(results)[h] <- nam
    
    results[[h+1]] <- eval(parse(text= nam2))
    names(results)[h+1] <- nam2
    
    h <- h + 2
    k <-k+1
  }
  
  
  if(!is.null(X_d)){
    h = 1
    k=1
    while(h < 2*dim(X_d)[3]){
      
      nam3 <- paste("beta_rMatDYAD", k, sep = "")
      nam4 <- paste("beta_cMatDYAD", k, sep = "")
      
      results[[(2*dim(X_r)[3] + h)]]   <- eval(parse(text= nam3))
      names(results)[(2*dim(X_r)[3] + h)] <- nam3
      
      results[[(2*dim(X_r)[3] + h+1)]] <- eval(parse(text= nam4))
      names(results)[(2*dim(X_r)[3] + h +1)] <- nam4
      
      h <- h + 2
      k <-k+1
    }
    image(U%*%t(U))
    
    
  }
  if(keep.UV == TRUE){
    results[[length(results)+1]] <- UV_mat
    names(results)[[length(results)]] <- "U"
  }
  results[[length(results)+1]] <-   Z_mat 
  names(results)[[length(results)]] <- "degcount"
  
  
  print("COUNTER")
  print(counter/MH)
  print(counter)
  return(list(results=results,counter=counter, MH=MH))
  
  
}

# Main function to run MCMC, all covariates are assumed to be dep. on community
#' @param X_r: n x n x p array of row covar
#' @param X_c: n x n x p array of col covar
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param Y: Observed network (can be binary or censored binary)
#' @param iter: number of iterations to run mcmc for 
#' @param indepBeta: which covariates are indep of community
#' @param numGroup: number of communities (also denote k)
#' @param prior_beta_mu: k dimensional vector with prior mean for beta
#' @param prior_beta_var:k x k matrix for prior variance for beta
#' @param prior_alpha_mu: constant for alpha(intercept) mean
#' @param prior_alpha_var:constant for alpha(intercept) var
#' @param start_beta_c: p_dep x k matrix for initializing beta_c
#' @param start_beta_r:  p_dep x k matrix for initializing beta_c
#' @param start_beta_dyad:  p_dep_dyad x k matrix for initializing beta_c
#' @param start_alpha: constant to initialize alpha
#' @param UandV: Default is FALSE, indicates if U and V should be calculated or just U
#' @param symmetricLambda: Default is TRUE, indicates if lambda should be symmetric
#' @param keep.UV: Default is TRUE, store U/V?
#' @param dcor: Default is TRUE, indicates if correlation in error terms should be estimated
#' @param symmetric: Default is FALSE, indicates that Y is not symmetric
#' @param model:  "cbin" or "bin", indicates what type of network Y is 
#' @param odmax: if model = "cbin", this is the max out degree a node can have
#' @param odens = 10: how often to save output
#' @param burnin: amount of burnin iterations

subAmen_Dep_And_IndepBeta<-function(X_r, X_c, X_d, Y,iter,indepBeta, numGroup, prior_beta_mu, prior_beta_var,
                                    prior_alpha_mu, prior_alpha_var, start_beta_c,
                                    start_beta_r,start_beta_dyad, start_alpha, 
                                    UandV = FALSE, symmetricLambda = TRUE,  keep.UV = TRUE,
                                    dcor =TRUE, symmetric = FALSE, model = "bin",odmax = NULL,odens = 10, 
                                    burnin = 1000){
  
  
  p_indep <- length(indepBeta)
  p_dep   <- dim(X_r)[3]- p_indep
  n <- dim(Y)[1]
  counter <- 0
  MH <-0
  GOF <- c()
  ######################  initialize
  beta_indep_r <- rep(0, p_indep)
  beta_indep_c <- rep(0, p_indep)
  beta_r <- start_beta_r
  beta_c <- start_beta_c
  
  if(!is.null(X_d)){
    beta_dc <- beta_dr <- start_beta_dyad
    
  }else{
    beta_dc <- beta_dr <- 0
    
  }
  
  
  if(symmetric == FALSE){
    rvar <- cvar <-  TRUE
    nvar <- FALSE
  }else{
    rvar <-cvar <- TRUE
    nvar <- TRUE
  }
  counter <- 0
  MH <-0
  
  
  Y_binary <- 1 * (Y > 0)
  
  
  if (is.element(model, c("cbin"))) {
    odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
    if (length(odmax) == 1) {
      odmax <- rep(odmax, nrow(Y))
    }
  }
  
  
  
  diag(Y) <- NA
  
  
  
  
  if (model == "bin") {
    Y <- 1 * (Y > 0)
    
    Z <- matrix(zscores(Y), nrow(Y), nrow(Y))
    z01 <- 0.5 * (max(Z[Y == 0], na.rm = TRUE) + min(Z[Y == 
                                                         1], na.rm = TRUE))
    Z <- Z - z01
  }
  if (is.element(model, c("cbin"))) {
    Z <- Y
    for (i in 1:nrow(Y)) {
      yi <- Y[i, ]
      zi <- zscores(yi)
      rnkd <- which(!is.na(yi) & yi > 0)
      if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
        zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 0.001
      }
      if (length(rnkd) < odmax[i]) {
        urnkd <- which(!is.na(yi) & yi == 0)
        if (max(zi[urnkd]) > 0) {
          zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 0.001
        }
      }
      Z[i, ] <- zi
    }
  }
  
  diag(Y) <- 0
  
  ######
  ## initialize U using regularizer from Rohe
  diag(Y_binary) <- 0
  UV_start <- specClust(Y_binary, tauAvg = TRUE,  numGroup,TRUE)
  U <- V <- matrix(UV_start$U_hat, ncol = numGroup)
  #V <-  matrix(UV_start$V_hat, ncol = numGroup)
  
  ## initialize lambda
  prior.lambda.mu  <- rep(1, numGroup^2)
  prior.lambda.var <- prior.lambda.prec <- diag(1, numGroup^2) ## have had .2
  lambda <- matrix(prior.lambda.mu, nrow = numGroup)
  
  ## alpha
  prior_alpha     <-  prior_alpha_mu
  prior_var_alpha <- prior_alpha_var
  alpha           <- start_alpha
  
  # get initial XBeta
  h_r_past <- NULL
  h_c_past <- NULL
  XB_indep_get <- XBETA_INDEP(X_r[,,indepBeta], X_c[,,indepBeta], beta_indep_r, beta_indep_c, p_indep)
  
  
  
  XB_dep_get   <- XBETA_DEP_partition(X_r[,,-indepBeta], X_c[,,-indepBeta],X_d, 
                                      beta_r, beta_c,
                                      beta_dr, beta_dc, U,V,
                                      updatedUV=TRUE, h_r_past, h_c_past)
  XB_ind     <-XB_dep_get$XBETA
  
  if(!is.null(X_d)){
    XB_dyad    <-XB_dep_get$XBETADYAD
    
    if(p_dep==1){
      XB         <- rowSums(XB_dyad, dims = 2) + 
        XB_ind + alpha + XB_indep_get# calculate XB
      
    }else{
      XB         <- rowSums(XB_dyad, dims = 2) + 
        rowSums(XB_ind, dims =2 )+ alpha + XB_indep_get# calculate XB
      
    }
    
  }else{
    if(p_dep==1){
      XB         <-  XB_ind+ alpha  + XB_indep_get# calculate XB
      
    }else{
      XB         <-  rowSums(XB_ind, dims =2 )+ alpha  + XB_indep_get# calculate XB
      
    }
    
    
  }
  
  E <- Z- XB - U%*%lambda%*%t(V)
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  prior <- list()
  if(model!="nrm")
  {  
    ydist<-table(Y)
    ymode<-as.numeric(names(ydist)[ ydist==max(ydist) ])[1] 
    ## eg, in a sparse binary network, ymode will be zero 
    YB<-1*(Y!=ymode) 
    ybar<-mean(YB,na.rm=TRUE) ; mu<-qnorm(ybar)
    E<- (YB - ybar)/dnorm(qnorm(ybar)) ; diag(E)<-0
    a<-rowMeans(E,na.rm=TRUE)  ; a[is.na(a)]<-0 
    b<-colMeans(E,na.rm=TRUE)  ; b[is.na(b)]<-0
    vscale<-mean(diag(cov(cbind(a,b))))
    PHAT<-pnorm(mu+outer(a,b,"+"))
    vdfmlt<-.25/mean(PHAT*(1-PHAT))
    if(is.null(prior$Sab0)){ prior$Sab0<-diag(2)*vscale }
    if(is.null(prior$eta0)){ prior$eta0<-round(4*vdfmlt) } 
    if(is.null(prior$kappa0)){ prior$kappa0<-round((2*numGroup+2)*vdfmlt) }  
    if(is.null(prior$g)){ prior$g<-sum(!is.na(Y)) }
  }
  
  Sab <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar)) 
  E <- E- outer(a,b, "+")
  s2<-1
  if(model=="nrm"){s2<-mean(E^2)}
  rho<-cor( c(E[upper.tri(E)]), c(t(E)[upper.tri(E)]) )*dcor  
  
  ## initialize things to store
  for(i in 1: dim(X_r)[3]){
    nam <- paste("beta_rMat", i, sep = "")
    assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
    
    namc <- paste("beta_cMat", i, sep = "")
    assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
    
    
    
  }
  
  if(!is.null(X_d)){
    for(i in 1:dim(X_d)[3]){
      
      nam <- paste("beta_rMatDYAD", i, sep = "")
      assign(nam, matrix(nrow = iter/odens, ncol = numGroup))
      
      namc <- paste("beta_cMatDYAD", i, sep = "")
      assign(namc, matrix(nrow = iter/odens, ncol = numGroup))
      
    }
  }
  
  for(i in 1: p_indep){
    nam5 <- paste("beta_rMatI", i, sep = "")
    assign(nam5, matrix(nrow = iter/odens, ncol = 1))
    
    nam6 <- paste("beta_cMatI", i, sep = "")
    assign(nam6, matrix(nrow = iter/odens, ncol =1))
  }
  
  # initialize ulv
  ulv <- U %*% lambda %*% t(V)
  lambda_mat <- matrix(nrow = iter/odens, ncol = numGroup^2)
  if(keep.UV == TRUE){
    UV_mat      <-  matrix(nrow = iter/odens, ncol = length(c(U)))
    
    
  }
  
  zmean <- matrix(nrow = iter/odens, ncol = 1)
  Ys <-matrix(nrow = iter/odens, ncol = 1)
  
  updatedUV <- TRUE
  XB_dep_get   <- XBETA_DEP_partition(X_r[,,-indepBeta], X_c[,,-indepBeta], X_d, beta_r, beta_c,beta_dr, beta_dc, U,V,
                                      updatedUV, h_r_past, h_c_past)
  XBETA_ROW     <-XB_dep_get$XBETAROW
  XBETA_COL     <-XB_dep_get$XBETACOL
  XB_indep_get <- XBETA_INDEP(X_r[,,indepBeta], X_c[,,indepBeta], beta_indep_r, beta_indep_c, p_indep)
  
  if(!is.null(X_d)){
    
    XBETA_DYAD    <- rowSums(XB_ind_get$XBETADYAD,dims = 2)
  }else{
    XBETA_DYAD <- 0
  }
  XB_ind        <- XBETA_ROW + XBETA_COL + XB_indep_get
  h_r_past      <-XB_dep_get$h_r_past
  h_c_past      <-XB_dep_get$h_c_past
  
  if(p_dep==1){
    XB            <- XB_ind+ XBETA_DYAD + alpha # calculate XB
    
  }else{
    XB            <- rowSums(XB_ind, dims =2 )+ XBETA_DYAD + alpha # calculate XB
    
  }
  
  diag(Y_binary) <-0
  ## start sampling
  diag(Z) <- mean(na.omit(c(Z)))
  Se   <-matrix(c(1,rho,rho,1),2,2)*s2
  iSe2 <-mhalf(solve(Se))
  td   <-iSe2[1,1] ; to<-iSe2[1,2]
  
  # j = 0
  for(i in 1:(iter+burnin)){
    
    # Draw latent lambda and U
    lambda <- drawLambda(Z, XB + outer(a,b, "+"), U,V, prior.lambda.prec,
                         prior.lambda.mu,td,to, symmetricLambda = symmetricLambda)
    
    UV_info      <- drawUV_probitLikelihood_indepDep(lambda, U, V,  XB ,a,b,
                                                     Y_binary, UandV,counter, MH,
                                                     X_r, X_c,X_d, beta_r, beta_c,alpha,
                                                     beta_dr, beta_dc,
                                                     h_r_past, h_c_past,Z,td,to,Sab, rho,indepBeta, XB_indep_get)
    
    
    U            <- UV_info$U
    V            <- UV_info$V
    
    
    updatedUV     <- UV_info$updateU
    swapPeep      <- UV_info$swapPeep
    counter       <- UV_info$counter
    MH            <- UV_info$MH
    
    ULV <- U %*%lambda%*%t(V)
    
    # get indep XB
    XB_indep_get <- XBETA_INDEP(X_r[,,indepBeta], X_c[,,indepBeta], beta_indep_r, beta_indep_c, p_indep)
    
    # get dep XB
    XB_dep_get   <- XBETA_DEP_partition(X_r[,,-indepBeta], X_c[,,-indepBeta],
                                        X_d, beta_r, beta_c,
                                        beta_dr, beta_dc, U,V,
                                        updatedUV, h_r_past, h_c_past)
    XBETA_ROW     <-XB_dep_get$XBETAROW
    XBETA_COL     <-XB_dep_get$XBETACOL
    if(!is.null(X_d)){
      
      XBETA_DYAD    <- rowSums(XB_dep_get$XBETADYAD,dims = 2)
    }else{
      XBETA_DYAD <- 0
    }
    
    XB_dep        <- XBETA_ROW + XBETA_COL 
    h_r_past      <-XB_dep_get$h_r_past
    h_c_past      <-XB_dep_get$h_c_past
    
    # recalculate XB
    if(p_dep==1){
      XB            <- XB_dep + alpha +XBETA_DYAD + XB_indep_get # calculate XB
      
    }else{
      XB            <- rowSums(XB_dep, dims =2 )+ alpha +XBETA_DYAD + XB_indep_get # calculate XB
      
    }
    
    EZ <- XB + ULV + outer(a, b, "+")
    
    
    if(updatedUV == TRUE){
      for(h in 1:length(swapPeep)){
        Z <-     update_Z_bin(swapPeep[h], Y, EZ, rho,  n,Z)
      }
      
    }
    
    
    # Draw rho
    
    if(symmetric){rho <- min(0.9999, 1 - 1/sqrt(i))}
    if(dcor){ rho <- rrho_mh_new(Z-EZ,rho,s2,offset=0,asp = NULL)} 
    
    ## for decorrelation
    Se   <-matrix(c(1,rho,rho,1),2,2)*s2
    iSe2 <-mhalf(solve(Se))
    td   <-iSe2[1,1] ; to<-iSe2[1,2]
    
    
    ## draw intercept
    XB_tmp <- XB -alpha
    
    alpha_tmp <-  rALPHA(Z, ULV,
                         alpha, a , b, XB,td, to,prior_alpha_mu,prior_alpha_var)
    print("alpha")
    
    print(alpha)
    alpha <- alpha_tmp$alpha 
    
    XB_updatedALPHA <- XB_tmp + alpha
    
    # get dyadic coef.
    if(!is.null(X_d)){
      
      numEigs <- 1 
      
      tmpDYAD <- rbeta_DECORR_DEPDYAD(Z-XB_updatedALPHA+XBETA_DYAD, X_d, U,V, lambda, 
                                      beta_dr, beta_dc, 
                                      prior_beta_mu, 
                                      prior_beta_var , a, b,   to ,td,numEigs, h_r_past,
                                      h_c_past)
      
      
      beta_dr <- tmpDYAD$beta_row_dyad
      beta_dc <- tmpDYAD$beta_col_dyad
      
      
      XBtmp <- XB_updatedALPHA -XBETA_DYAD
      XB_updateDyad_get    <- XBETA_DEP_UPDATEDDYAD(X_d, 
                                                    beta_dr, beta_dc,
                                                    U,V)
      XBETA_DYAD    <- rowSums(XB_updateDyad_get$XBETADYAD, dims =2 )
      XB            <-  XBtmp +XBETA_DYAD # calculate XB
      EZ <- XB + ULV + outer(a, b, "+")
      
      
    }else{
      XB <- XB_updatedALPHA
      EZ <- XB + ULV + outer(a, b, "+")
      
    }
    
    
    
    
    # Draw Beta
    ## DRAW INDEP BETA
    if(!is.null(indepBeta)){
      
      XB_dep_get   <- XBETA_DEP_partition(X_r[,,-indepBeta], X_c[,,-indepBeta],
                                          X_d, beta_r, beta_c,
                                          beta_dr, beta_dc, U,V,
                                          updatedUV, h_r_past, h_c_past)
      
      indepBeta_samp <- rbeta_INDEPDECORR(Z, X_r[,,indepBeta], X_c[,,indepBeta], U,lambda,
                                          p_indep, beta_indep_r, beta_indep_c,
                                          alpha, a, b,
                                          XB_dep_get$XBETAROW +XB_dep_get$XBETACOL,td,to)
      
    }
    
    beta_indep_c <- indepBeta_samp$beta_col
    beta_indep_r <- indepBeta_samp$beta_row
    
    XB_indep_get <- XBETA_INDEP(X_r[,,indepBeta], X_c[,,indepBeta], beta_indep_r, beta_indep_c, p_indep)
    
    XB            <- XB_dep + alpha +XBETA_DYAD + XB_indep_get # calculate XB
    
    
    
    tmp <- rbeta_DECORR_DEP_ANDINDEP(Z, X_r[,,-indepBeta], X_c[,,-indepBeta], 
                                     U,lambda, p_dep, beta_r, beta_c,
                                     alpha, prior_beta_mu, 
                                     prior_beta_var , a, b,
                                     h_r_past, h_c_past, to ,td,updatedU, 
                                     XB_indep_get )
    
    
    beta_r <- tmp$beta_row
    beta_c <- tmp$beta_col
    
    
    
    if(p_dep==1){
      XBtmp <- XB - XB_dep
      
    }else{
      XBtmp <- XB - rowSums(XB_dep, dims =2 )
      
      
    }
    
    # recalculate XB
    
    
    XB_updateRC <- XBETA_DEP_UPDATEROWCOLUMN(X_r[,,-indepBeta], X_c[,,-indepBeta],
                                             beta_r, beta_c, 
                                             U,V,
                                             FALSE, 
                                             h_r_past, 
                                             h_c_past)
    
    XBETA_ROW     <-  XB_updateRC$XBETAROW
    XBETA_COL     <-  XB_updateRC$XBETACOL
    XB_dep        <-  XBETA_ROW + XBETA_COL# +XB_indep_get
    
    if(p_dep ==1){
      XB            <-  XBtmp + XB_dep# -XB_indep_get+alpha
      
    }else{
      XB            <-  XBtmp + rowSums(XB_dep, dims =2)#-XB_indep_get+alpha
      
      
    }
    #
    EZ <- XB + ULV + outer(a, b, "+")
    
    
    
    
    
    
    
    # draw Z
    if(model == "bin"){Z  <- rZ_bin_fc_new(Z, EZ, rho, Y) }
    if(model=="cbin"){Z<-rZ_cbin_fc(Z,EZ,rho,Y,odmax,odobs)}
    
    
    # draw random effects
    ab <- r_ab_DECORR_DEP(Z, XB, Sab,rho, s2, ULV)
    
    a  <- ab$a * rvar
    b  <- ab$b * cvar
    
    if(symmetric == TRUE){
      a <- b <- (a + b)/2
    }
    EZ <- XB + ULV + outer(a, b, "+")
    
    # draw Sab and a -> z
    if(symmetric == FALSE){
      Sab <-rSab_fc(a,b,Sab0=prior$Sab0,eta0=prior$eta0)
      if(model=="frn"){
        tmp<-raSab_frn_fc_new(Z,Y,YL,a,b,Sab,odmax, odobs,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
        
      }
      if(model=="bin"){
        tmp<-raSab_bin_fc_new(Z,Y,a,b,Sab,Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      if(model=="cbin")
      {
        tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs,
                           Sab0=prior$Sab0,eta0=prior$eta0)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      }
      
    } else{
      Sab[1,1]<-Sab[2,2]<-1/rgamma(1,(1+nrow(Y))/2,(1+sum(a^2))/2)
      Sab[1,2]<-Sab[2,1]<-.999*Sab[1,1]
    }
    
    
    
    # posterior pred
    if(symmetric){ EZ<-(EZ+t(EZ))/2 } 
    
    
    
    if(model=="bin") { Ys_draw<-simY_bin(EZ,rho) }
    if(model=="cbin"){ Ys<-1*(simY_frn(EZ,rho,odmax,YO=Y)>0) }
    
    # Store results
    if(i == burnin + 1){ j = 0}
    if(i > burnin){
      
      if(i %% odens==0){
        j = j + 1
        
        for(g in 1:dim(beta_r)[1]){
          nam <-sprintf("beta_rMat%s[%s,]", g, j)
          nam2 <-sprintf("beta_cMat%s[%s,]", g, j)
          
          eval(parse(text= paste(nam, "<- beta_r[g,]")))        
          eval(parse(text= paste(nam2, "<- beta_c[g,]")))        
          
          
          
        }
        
        
        for(g in 1:p_indep){
          nam <-sprintf("beta_rMatI%s[%s,]", g, j)
          nam2 <-sprintf("beta_cMatI%s[%s,]", g, j)
          
          eval(parse(text= paste(nam, "<- beta_indep_r[g]")))        
          eval(parse(text= paste(nam2, "<- beta_indep_c[g]")))        
          
          
          
          
        }
        if(!is.null(X_d)){
          
          for(k in 1:dim(X_d)[3]){
            nam3 <-sprintf("beta_rMatDYAD%s[%s,]", k, j)
            nam4 <-sprintf("beta_cMatDYAD%s[%s,]", k, j)
            
            eval(parse(text= paste(nam3, "<- beta_dr[k,]")))        
            eval(parse(text= paste(nam4, "<- beta_dc[k,]")))        
            
          }
        }
        
        
        lambda_mat[j,] <- c(lambda)
        zmean[j,] <- mean(c(Z))
        Ys[j,] <-sum(abs(Ys_draw-Y),na.rm = TRUE)
        
        if(keep.UV == TRUE){
          UV_mat[j, ] <- c(U)
        }
        
      }
      
    }
    
  }
  
  ## make list to return
  
  results <- list()
  
  h = 1
  k=1
  while(h < 2*dim(beta_r)[1]){
    
    nam <- paste("beta_rMat", k, sep = "")
    nam2 <- paste("beta_cMat", k, sep = "")
    
    results[[h]]   <- eval(parse(text= nam))
    names(results)[h] <- nam
    
    results[[h+1]] <- eval(parse(text= nam2))
    names(results)[h+1] <- nam2
    
    h <- h + 2
    k <-k+1
  }
  
  
  
  if(!is.null(X_d)){
    h = 1
    k=1
    while(h < 2*dim(X_d)[3]){
      
      nam3 <- paste("beta_rMatDYAD", k, sep = "")
      nam4 <- paste("beta_cMatDYAD", k, sep = "")
      
      results[[(2*dim(beta_r)[1] + h)]]   <- eval(parse(text= nam3))
      names(results)[(2*dim(X_r)[3] + h)] <- nam3
      
      results[[(2*dim(beta_r)[1] + h+1)]] <- eval(parse(text= nam4))
      names(results)[(2*dim(beta_r)[1] + h +1)] <- nam4
      
      h <- h + 2
      k <-k+1
    }
  }
  
  h = 1
  k=1
  
  if(is.null(X_d)){
    currDim <- (2*dim(beta_r)[1]) 
    
  }else{
    currDim <- (2*dim(beta_r)[1]+ 2*dim(X_d)[3]) 
    
  }
  while(h < 2*p_indep){
    
    nam5 <- paste("beta_rMatI", k, sep = "")
    nam6 <- paste("beta_cMatI", k, sep = "")
    
    results[[(currDim + h)]]   <- eval(parse(text= nam5))
    names(results)[(currDim + h)] <- nam5
    
    results[[(currDim + h+1)]] <- eval(parse(text= nam6))
    names(results)[(currDim + h+1)] <- nam6
    
    h <- h + 2
    k <-k+1
  }
  
  if(keep.UV == TRUE){
    results[[length(results)+1]] <- UV_mat
    names(results)[[length(results)]] <- "U"
  }
  
  print("COUNTER")
  print(counter/MH)
  print(counter)
  return(list(results=results,counter=counter, MH=MH))
  
  
}




# binary version so we are drawing from the BINARY likelihood
#' @param idx: index of node to update
#' @param Y: Observed network (can be binary or censored binary)
#' @param EZ: current expected value of Z
#' @param rho: correlation parameter
#' @param n: size of network
#' @param Z: current Z matrix 
update_Z_bin <- function(idx, Y, EZ, rho,  n,Z){
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2, byrow = T)
  
  # we need to change rows that had a U switch
  # first, get values for rows where Y equals 1 for both
  for(j in 1:n){
    if(Y[idx,j]==1 & Y[j, idx]==1){
      lb <- rep(0,2); ub <- rep(Inf, 2)  
      tmpZ_pair <-  rtmvnorm(1, c(EZ[idx,j], EZ[j,idx]), Sigma, lb, ub)
      Z[idx, j] <- tmpZ_pair[1]
      Z[j, idx] <- tmpZ_pair[2]
      
    }else if(Y[idx,j]==1 & Y[j, idx]==0){
      
      
      lb <- c(0, -Inf); ub <- c(Inf, 0)
      
      
      tmpZ_pair <-  rtmvnorm(1, c(EZ[idx,j], EZ[j,idx]), Sigma, lb, ub)
      Z[idx, j] <- tmpZ_pair[1]
      Z[j, idx] <- tmpZ_pair[2]
      
      
    }else if(Y[idx,j]==0 & Y[j, idx]==1){
      # check if censoring occurs
      
      lb <- c(-Inf, 0 ); ub <- c( 0, Inf)
      
      
      tmpZ_pair <-  rtmvnorm(1, c(EZ[idx,j], EZ[j,idx]), Sigma, lb, ub)
      Z[idx, j] <- tmpZ_pair[1]
      Z[j, idx] <- tmpZ_pair[2]
    }else{
      
      
      lb <- rep(-Inf, 2); ub <- rep(0, 2)
      
      tmpZ_pair <-  rtmvnorm(1, c(EZ[idx,j], EZ[j,idx]), Sigma, lb, ub)
      Z[idx, j] <- tmpZ_pair[1]
      Z[j, idx] <- tmpZ_pair[2]
    }
    
  }
  return(Z)
  
  
}

# Function to calculate dep X beta that is NOT all added together (returns array of XB)
#' @param X_row: n x n x p_dep array of row covar
#' @param X_col: n x n x p_dep array of col covar
#' @param X_dyad: n x n x p_dep_dyad array of dyadic covar
#' @param beta_r: current beta_r, p_dep x k matrix
#' @param beta_c: current beta_c, p_dep x k matrix
#' @param beta_dr: current beta_dr, p_dep_dyad x k matrix
#' @param beta_dc: current beta_dc, p_dep_dyad x k matrix
#' @param U: current U (n x k matrix)
#' @param V: curren V (n x k matrix)
#' @param updatedUV: indicates if U, V have been updated. If so, H_r/H_c need to be re calculated
#' @param h_r_past: Previous calculation of h_r
#' @param h_c_past: Previous calculation of h_c
#' @param indepBeta: indicates if there is an independent covariate 



XBETA_DEP_allowIndep<- function(X_row, X_col, beta_r, beta_c, U,num_covar, updatedU, h_r_past, h_c_past,indepBeta){
  n = dim(X_row)[1]
  XBETA <- array(dim = c(n,n, num_covar))
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  if(updatedU == TRUE){
    H_R <- list()
    H_C <- list()
  }else{
    H_R <- h_r_past
    H_C <- h_c_past
  }
  
  
  if(num_covar ==1){
    if(updatedU == TRUE){
      h_r <- (ones %x% ((diag(1,n) *X_row)%*%U)) ## for row
      h_c <- ((diag(1,n)*X_col)%*% U) %x% ones
      
      ind <- ifelse(indepBeta == 1, 2, 1)
      H_R[[ind]] <- h_r
      H_C[[ind]] <- h_c
    }else{
      h_r <- H_R[[ind]]
      h_c <- H_C[[ind]]
    }
    
    
    
    
    rowXB <- h_r %*% c(beta_r)
    colXB <- h_c %*% c(beta_c)
    XBETA <- matrix(rowXB,nrow = n) + matrix(colXB,nrow = n)
    
  }else{
    for(p in 1:num_covar){
      ## calculate XB
      if(updatedU == TRUE){
        h_r <- (ones %x% ((diag(1,n) *X_row[,,p])%*%U)) ## for row
        h_c <- ((diag(1,n)*X_col[,,p])%*% U) %x% ones
        
        H_R[[p]] <- h_r
        H_C[[p]] <- h_c
      }else{
        h_r <- H_R[[p]]
        h_c <- H_C[[p]]
      }
      
      
      rowXB <- h_r %*% beta_r[p,]
      colXB <- h_c %*% beta_c[p,]
      XBETA[,,p] <- matrix(rowXB,nrow = n) + matrix(colXB,nrow = n)
      
    }
  }
  
  return(list(XBETA = XBETA, h_r_past = H_R, h_c_past = H_C))
  
}


# Function to calculate dep X beta that is NOT all added together 
# (returns array of XBETAROW, COL, DYAD )
#' @param X_row: n x n x p_dep array of row covar
#' @param X_col: n x n x p_dep array of col covar
#' @param X_dyad: n x n x p_dep_dyad array of dyadic covar
#' @param beta_r: current beta_r, p_dep x k matrix
#' @param beta_c: current beta_c, p_dep x k matrix
#' @param beta_dr: current beta_dr, p_dep_dyad x k matrix
#' @param beta_dc: current beta_dc, p_dep_dyad x k matrix
#' @param U: current U (n x k matrix)
#' @param V: curren V (n x k matrix)
#' @param updatedUV: indicates if U, V have been updated. If so, H_r/H_c need to be re calculated
#' @param h_r_past: Previous calculation of h_r
#' @param h_c_past: Previous calculation of h_c


XBETA_DEP_partition <-function(X_row, X_col, X_dyad, 
                               beta_r, beta_c, 
                               beta_dr, beta_dc,
                               U,V,
                               updatedUV, 
                               h_r_past, 
                               h_c_past){
  
  n = dim(X_row)[1]
  num_covar = dim(beta_r)[1]
  
  if(is.null(num_covar)){
    num_covar <-1
  }
  
  if(!is.null(X_dyad)){
    num_covarDYAD = dim(beta_dr)[1]
    XBETADYAD <- array(dim = c(n,n, num_covarDYAD))
    
  }
  
  XBETAROW <- array(dim = c(n,n, num_covar))
  XBETACOL <- array(dim = c(n,n, num_covar))
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  if(updatedUV == TRUE){
    H_R <- list()
    H_C <- list()
    
  }else{
    H_R <- h_r_past
    H_C <- h_c_past
    
  }
  
  
  ## For row and column covariates
  if(num_covar ==1){
    X_r <- X_row
    X_c <- X_col
    if(updatedUV == TRUE){
      h_r <- (ones %x% ((diag(1,n)*(X_r))%*%U)) ## for row
      h_c <- ((diag(1,n)*(X_c))%*% V) %x% ones
      
      H_R[[1]] <- h_r
      H_C[[1]] <- h_c
    }else{
      h_r <- H_R[[1]]
      h_c <- H_C[[1]]
    }
    
    
    
    
    rowXB <- h_r %*% matrix(beta_r,nrow = dim(h_r)[2])#[1,]
    colXB <- h_c %*% matrix(beta_c,nrow = dim(h_r)[2])#[1,]
    XBETAROW <- matrix(rowXB,nrow = n)
    XBETACOL <- matrix(colXB,nrow = n)
    
    
  }else{
    for(p in 1:num_covar){
      X_r <- X_row[,,p]
      X_c <- X_col[,,p]
      
      ## calculate XB
      if(updatedUV == TRUE){
        h_r <- (ones %x% ((diag(1,n) *(X_r))%*%U)) ## for row
        h_c <- ((diag(1,n)*(X_c))%*% V) %x% ones
        
        H_R[[p]] <- h_r
        H_C[[p]] <- h_c
      }else{
        h_r <- H_R[[p]]
        h_c <- H_C[[p]]
      }
      
      
      
      
      rowXB <- h_r %*% beta_r[p,]
      colXB <- h_c %*% beta_c[p,]
      XBETAROW[,,p] <- matrix(rowXB,nrow = n)
      XBETACOL[,,p] <- matrix(colXB,nrow = n)
      
      
    }
    
  }
  
  
  if(!is.null(X_dyad)){
    for(p in 1:num_covarDYAD){
      X_d <- X_dyad[,,p]
      dyadXB <- diag(c(U%*%(beta_dr[p,]))) %*% X_d %*%  diag(c(V%*%(beta_dc[p,])))
      
      XBETADYAD[,,p] <- matrix(dyadXB,nrow = n)
      
      
    }
  }else{
    XBETADYAD <- NULL
  }
  
  return(list(XBETAROW = XBETAROW, 
              XBETACOL = XBETACOL, 
              XBETADYAD = XBETADYAD,
              h_r_past = H_R, h_c_past = H_C))
  
}

# Function to calculate dep X beta just for dyadic covar (returns array of XB)
#' @param X_dyad: n x n x p_dep_dyad array of dyadic covar
#' @param beta_dr: current beta_dr, p_dep_dyad x k matrix
#' @param beta_dc: current beta_dc, p_dep_dyad x k matrix


XBETA_DEP_UPDATEDDYAD <-function(X_dyad, 
                                 beta_dr, beta_dc,
                                 U,V){
  
  n = dim(X_dyad)[1]
  num_covarDYAD = dim(beta_dr)[1]
  
  
  XBETADYAD <- array(dim = c(n,n, num_covarDYAD))
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  
  
  ## for dyadic covariates
  for(p in 1:num_covarDYAD){
    X_d    <- X_dyad[,,p]
    dyadXB <- diag(c(U%*%(beta_dr[p,]))) %*% X_d %*%  diag(c(V%*%(beta_dc[p,])))
    
    XBETADYAD[,,p] <- matrix(dyadXB,nrow = n)
    
    
  }
  return(list(XBETADYAD = XBETADYAD))
  
}

# Function to calculate dep X beta JUST for row and column and returns XBETAROW, XBETACOL separate
#' @param X_row: n x n x p_dep array of row covar
#' @param X_col: n x n x p_dep array of col covar
#' @param beta_r: current beta_r, p_dep x k matrix
#' @param beta_c: current beta_c, p_dep x k matrix
#' @param U: current U (n x k matrix)
#' @param V: curren V (n x k matrix)
#' @param updatedUV: indicates if U, V have been updated. If so, H_r/H_c need to be re calculated
#' @param h_r_past: Previous calculation of h_r
#' @param h_c_past: Previous calculation of h_c
XBETA_DEP_UPDATEROWCOLUMN <-function(X_row, X_col, 
                                     beta_r, beta_c, 
                                     U,V,
                                     updatedUV, 
                                     h_r_past, 
                                     h_c_past){
  
  n = dim(X_row)[1]
  num_covar = dim(beta_r)[1]
  
  XBETAROW <- array(dim = c(n,n, num_covar))
  XBETACOL <- array(dim = c(n,n, num_covar))
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  if(updatedUV == TRUE){
    H_R <- list()
    H_C <- list()
    
  }else{
    H_R <- h_r_past
    H_C <- h_c_past
    
  }
  
  ## For row and column covariates
  
  for(p in 1:num_covar){
    
    if(num_covar==1){
      X_r <- X_row
      X_c <- X_col
      
    }else{
      X_r <- X_row[,,p]
      X_c <- X_col[,,p]
      
    }
    
    
    ## calculate XB
    if(updatedUV == TRUE){
      h_r <- (ones %x% ((diag(1,n) *(X_r))%*%U)) ## for row
      h_c <- ((diag(1,n)*(X_c))%*% V) %x% ones
      
      H_R[[p]] <- h_r
      H_C[[p]] <- h_c
    }else{
      h_r <- H_R[[p]]
      h_c <- H_C[[p]]
    }
    
    
    
    
    rowXB <- h_r %*% beta_r[p,]
    colXB <- h_c %*% beta_c[p,]
    
    if(num_covar==1){
      XBETAROW <- matrix(rowXB,nrow = n)
      XBETACOL <- matrix(colXB,nrow = n)
      
    }else{
      XBETAROW[,,p] <- matrix(rowXB,nrow = n)
      XBETACOL[,,p] <- matrix(colXB,nrow = n)
      
    }
    
    
  }
  
  
  
  return(list(XBETAROW = XBETAROW, 
              XBETACOL = XBETACOL, 
              h_r_past = H_R, h_c_past = H_C))
  
}


# Function to calculate dep X beta that is NOT all added together (returns array of XB)
#' @param X_row: n x n x p_dep array of row covar
#' @param X_col: n x n x p_dep array of col covar
#' @param X_dyad: n x n x p_dep_dyad array of dyadic covar
#' @param beta_r: current beta_r, p_dep x k matrix
#' @param beta_c: current beta_c, p_dep x k matrix
#' @param beta_dr: current beta_dr, p_dep_dyad x k matrix
#' @param beta_dc: current beta_dc, p_dep_dyad x k matrix
#' @param U: current U (n x k matrix)
#' @param V: curren V (n x k matrix)
#' @param updatedUV: indicates if U, V have been updated. If so, H_r/H_c need to be re calculated
#' @param h_r_past: Previous calculation of h_r
#' @param h_c_past: Previous calculation of h_c
#' @param indepBeta: indicates if there is an independent covariate 
XBETA_DEP <- function(X_row, X_col,
                      X_dyad,
                      beta_r, beta_c,
                      beta_dr, beta_dc, 
                      U, V,  updatedUV,
                      h_r_past, h_c_past,indepBeta){
  
  
  n = dim(X_row)[1]
  num_covar <- dim(beta_r)[1]
  XBETA <- array(dim = c(n,n, num_covar))
  
  
  if(!is.null(X_dyad)){
    num_covarDYAD <- dim(beta_dr)[1]
    
    XBETADYAD <- array(dim = c(n,n, num_covarDYAD))      
  }else{
    XBETADYAD <- 0 
    
  }
  
  
  ones <- as.matrix(rep(1,n), ncol =1 )
  
  if(updatedUV == TRUE){
    H_R <- list()
    H_C <- list()
  }else{
    H_R <- h_r_past
    H_C <- h_c_past
  }
  
  
  if(num_covar ==1){
    if(updatedUV == TRUE){
      h_r <- (ones %x% ((diag(1,n) *X_row)%*%U)) ## for row
      h_c <- ((diag(1,n)*X_col)%*% V) %x% ones  ## col
      
      ind <- ifelse(indepBeta == 1, 2, 1)
      H_R[[ind]] <- h_r
      H_C[[ind]] <- h_c
    }else{
      h_r <- H_R[[ind]]
      h_c <- H_C[[ind]]
    }
    
    
    
    
    rowXB <- h_r %*% c(beta_r)
    colXB <- h_c %*% c(beta_c)
    XBETA <- matrix(rowXB,nrow = n) + matrix(colXB,nrow = n)
    
  }else{
    for(p in 1:num_covar){
      ## calculate XB
      if(updatedUV == TRUE){
        h_r <- (ones %x% ((diag(1,n) *X_row[,,p])%*%U)) ## for row
        h_c <- ((diag(1,n)*X_col[,,p])%*% V) %x% ones ## col
        
        H_R[[p]] <- h_r
        H_C[[p]] <- h_c
      }else{
        h_r <- H_R[[p]]
        h_c <- H_C[[p]]
      }
      
      
      
      
      rowXB <- h_r %*% beta_r[p,]
      colXB <- h_c %*% beta_c[p,]
      XBETA[,,p] <- matrix(rowXB,nrow = n) + matrix(colXB,nrow = n)
      
    }
  }
  
  
  if(!is.null(X_dyad)){
    ## for dyadic covariates
    for(p in 1:num_covarDYAD){
      X_d <- X_dyad[,,p]
      dyadXB <- diag(c(U%*%(beta_dr[p,]))) %*% X_d %*%  diag(c(V%*%(beta_dc[p,])))
      
      XBETADYAD[,,p] <- matrix(dyadXB,nrow = n)
      
      
    }
  }
  
  
  
  
  return(list(XBETA = XBETA,XBETADYAD = XBETADYAD,
              h_r_past = H_R, h_c_past = H_C))
  
  
}


# Function to calculate X beta that is NOT all added together (returns array of XB)
#' @param X_row: n x n x p_indep array of row covar
#' @param X_col: n x n x p_indep array of col covar
#' @param beta_r: current beta_r, p_indep vector 
#' @param beta_c: current beta_c, p_indep vector 
#' @param num_covar: number of indep covariates

Xbeta_indep_Individual <-function(X_row, X_col, beta_r, beta_c, num_covar){
  n = dim(X_row)[1]
  XBETA <- array(dim = c(n,n, num_covar))
  
  if(num_covar == 1){
    
    rowXB <- beta_r *X_row 
    colXB <- X_col  *beta_c
    XBETA <- rowXB + colXB
  }else{
    for(p in 1:num_covar){
      ## calculate XB
      rowXB <- beta_r[p] *X_row[,,p] 
      colXB <- X_col[,,p]  *beta_c[p]
      XBETA[,,p] <- rowXB + colXB
      
    }
    
  }
  
  return(XBETA)
  
}


# Function to calculate X beta for indep specified covariate
#' @param X_row: n x n x p_indep array of row covar
#' @param X_col: n x n x p_indep array of col covar
#' @param beta_r: current beta_r, p_indep vector 
#' @param beta_c: current beta_c, p_indep vector 
#' @param num_covar: number of indep covariates

XBETA_INDEP <- function(X_row, X_col, beta_r, beta_c, num_covar){
  
  n <- dim(X_row)[1]
  XBETA <- matrix(0,n,n)
  
  if(num_covar == 1){
    rowXB <- beta_r *X_row
    colXB <- X_col *beta_c
    XB <- rowXB + colXB
    XBETA <- XB
    
  }else{
    for(p in 1:num_covar){
      ## calculate XB
      
      rowXB <- beta_r[p] *X_row[,,p] 
      colXB <- X_col[,,p] *beta_c[p]
      XB <- rowXB + colXB
      XBETA <- XBETA + XB
      
    }
  }
  
  
  return(XBETA)
  
}



# covariate assisted spectral clustering from rho to give initialization for
# community membership
#' @param Y: observed network
#' @param X_r: array of row covariates  (n x  n x num_row covar)
#' @param X_c: array of row covariates  (n x  n x num_col covar)
#' @param tauAvg: adjust for average row sums (out degree) set to TRUE
#' @param covarAssist: should this be covariate assisted (always TRUE)
#' @param groupNum: number of communities to look for 

specClust_xx         <- function(Y,X_r,X_c, tauAvg = TRUE, covarAssist = TRUE, groupNum  ){
  
  
  X<-c()
  for(i in 1:dim(X_r)[3]){
    x_cur <- X_r[,1,i]
    X <- cbind(X,x_cur)
    
  }
  
  for(i in 1:dim(X_c)[3]){
    x_cur2 <- X_c[1,,i]
    X <- cbind(X,x_cur2)
    
  }
  D <- rowSums(Y)
  D2 <- diag(D)
  n <- dim(Y)[1]
  tau2 <- ifelse(tauAvg == TRUE, sum(D2)/n, 0)
  
  newD <- D + tau2
  D_sqrt <- diag(1/sqrt(newD))
  L <- D_sqrt %*% Y %*% D_sqrt
  
  if(covarAssist == TRUE){
    XXt <- X %*% t(X)
    
    newL <- L + XXt
  } else{
    newL <- L
  }
  
  
  
  ## Find eigen values and vectors of L
  eigenDecomp        <- irlba(newL, nu = groupNum, nv = groupNum)
  
  
  
  communityClust <- kmeans(eigenDecomp$u, centers = groupNum)
  
  
  
  memberShip <- as.factor(communityClust$cluster)
  memberShip <- as.data.frame(memberShip)
  
  
  
  memberModel <- model.matrix(~memberShip-1, memberShip)
  U_hat <- data.matrix(memberModel)
  
  U_hat 
}


########## estimate U/V using Rohe spectral clustering (not covar assisted)
#' @param A: observed network
#' @param tauAvg: account for average row/col sums (out degree) (set to TRUE)
#' @param groupNum: number of communities
#' @param UandV: Should U and V be estimated or just U?
specClust           <- function(A, tauAvg = TRUE, groupNum, UandV){
  
  n <- dim(A)[1]
  DR <- rowSums(A)
  DC <- colSums(A)
  DR2 <- diag(DR)
  DC2 <- diag(DC)
  tau2R <- ifelse(tauAvg == TRUE, sum(DR2)/n, 0)
  tau2C <- ifelse(tauAvg == TRUE, sum(DC2)/n, 0)
  
  
  newDR <- DR+ tau2R
  newDC <- DC + tau2C
  
  DR_sqrt <- diag(1/sqrt(newDR))
  DC_sqrt <- diag(1/sqrt(newDC))
  
  L <- DR_sqrt %*% A %*% DC_sqrt
  
  ## Find eigen values and vectors of L
  eigenDecomp        <- irlba(L, nu = groupNum, nv = groupNum)
  
  
  
  communityClustU <- kmeans(eigenDecomp$u, centers = groupNum)
  
  
  
  
  memberShipU <- as.factor(communityClustU$cluster)
  memberShipU <- as.data.frame(memberShipU)
  
  
  
  memberModelU <- model.matrix(~memberShipU-1, memberShipU)
  U_hat <- data.matrix(memberModelU)
  colSums(U_hat)
  if(UandV == TRUE){
    communityClustV <- kmeans(eigenDecomp$v, centers = groupNum)
    memberShipV <- as.factor(communityClustV$cluster)
    memberShipV <- as.data.frame(memberShipV)
    
    memberModelV<- model.matrix(~memberShipV-1, memberShipV)
    V_hat <- data.matrix(memberModelV)
    results <- list(U_hat = U_hat, V_hat = V_hat) 
  }
  else{
    results <- list(U_hat = U_hat) 
  }
  
  return(results)
}


# function to estimate U using probit regression residuals
# this function was very specific to this example
# vectorize Y and all covariates, give it a k x k lambda
#' @param y_vect vectorized Y matrix
#' @param vectX1 vectorized  covariate X1 matrix
#' @param vectX1_2 vectorized  covariate X1 matrix
#' @param vectX2 vectorized  covariate X2 matrix
#' @param vectX2_2 vectorized  covariate X2 matrix
#' @param lambda kxk matrix to specify k and to later make ULUT
probitUEstimation <- function(y_vect, vectX1, vectX1_2,
                              vectX2,
                              vectX2_2, lambda){
  
  n <- sqrt(length(y_vect))
  k <- dim(lambda)[1]
  df <- cbind(y_vect, vectX1, vectX1_2, vectX2, vectX2_2)
  
  df <- as.data.frame(df)
  #### 
  proModel <- glm(y_vect ~., data = df, family = binomial(link = "probit"))
  
  
  
  residMatrix <- matrix(proModel$residuals, nrow = n, ncol = n)
  eigenDecomp <- eigen(residMatrix)
  eigenValues <- eigenDecomp$values
  eigenL_valuesOrder <- order(abs(eigenValues), decreasing = TRUE)
  
  
  # order the eigen vectors based on ordered eigen values (abs value)
  eigen_vectorsOrdered <- eigenDecomp$vectors[,eigenL_valuesOrder]
  
  # get the k important vectors
  kEigenVectors <- eigen_vectorsOrdered[,1:k]
  
  # let X be the only important ones
  X <- kEigenVectors
  
  
  X <- X[,1:k]
  
  communityClust <- kmeans(X, centers = k)
  
  
  
  memberShip <- as.factor(communityClust$cluster)
  memberShip <- as.data.frame(memberShip)
  
  
  
  memberModel <- model.matrix(~memberShip-1, memberShip)
  U_hatR5 <- data.matrix(memberModel)
  
  ULU <- U_hatR5 %*% lambda %*% t(U_hatR5)
  return(ULU) 
}


# not used, additional code for note
fullSimulation_amen_master <-function(data, anyIndep, family,odmax=0,burnin,iter, clusterNum, dataName){
  set.seed(1312)
  
  numGroup <- 3
  
  
  data_prior     <- data$priorSet
  data_initial   <- data$initialSet
  tic()
  set.seed(14)
  
  
  tempTimeo <- toc()
  
  ## dep amen
  if(clusterNum == 1){
    if(anyIndep == "all"){
      Model <- "depALL"
      amen_dep <- amen_master(data$X_r, data$X_c, data$Y,iter = iter,numGroup=3,
                              prior_beta_mu= data_prior$prior_beta_mu,
                              prior_beta_var = data_prior$prior_beta_var,
                              prior_alpha_mu= data_prior$prior_alpha_mu,
                              prior_alpha_var= data_prior$prior_alpha_var,
                              start_beta_c= matrix(rep(1,6), nrow = 2, ncol = 3),
                              start_beta_r = matrix(rep(1,6), nrow = 2, ncol = 3),
                              start_alpha= data_initial$start_alpha, keep.U = TRUE, dcor =TRUE,
                              symmetric = FALSE,
                              model = family, odmax = odmax,indepBeta=1, "all", burnin)
      int <-list(betaC1 = amen_dep$results$beta_cMat1, betaR1 = amen_dep$results$beta_rMat1,betaR2 = amen_dep$results$beta_rMat2,
                 betaC2 = amen_dep$results$beta_cMat2)
      
    }else{
      Model <- "depSome"
      amen_dep <-amen_master(data$X_r, data$X_c, data$Y,iter = iter,numGroup=3,
                             prior_beta_mu= data_prior$prior_beta_mu,
                             prior_beta_var = data_prior$prior_beta_var,
                             prior_alpha_mu= data_prior$prior_alpha_mu,
                             prior_alpha_var= data_prior$prior_alpha_var,
                             start_beta_c= data_initial$start_beta_c,
                             start_beta_r = data_initial$start_beta_r,
                             start_alpha= data_initial$start_alpha, keep.U = TRUE, dcor =TRUE,
                             symmetric = FALSE,
                             model = family, odmax = odmax,indepBeta=1,"some", burnin)
      int <-list(betaC1 = amen_dep$results$beta_cMat1, betaR1 = amen_dep$results$beta_rMat1,betaR2 = amen_dep$results$beta_rMatI1,
                 betaC2 = amen_dep$results$beta_cMatI1)
    }
    
  }
  
  
  # R = 3
  if(clusterNum == 2){
    Model <- "R3"
    amen_2_R_3_a <-  ame(data$Y, Xdyad = NULL, Xrow = data$df ,Xcol = data$df, dcor = TRUE,
                         R=3, cvar = TRUE,rvar = TRUE, nscan = iter, symmetric = FALSE, 
                         model = family,burn = burnin,odens = 10,plot = FALSE)
    colnames(amen_2_R_3_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")
    
    int <- mcmc_intervals_data(as.mcmc(amen_2_R_3_a$BETA), prob_outer = .95, point_est = "mean")[,-1]
    
    
  }
  
  
  
  # R = 0
  if(clusterNum==3){
    Model <- "R0"
    amen_2_R_0_a <-  ame(data$Y, Xdyad = NULL, Xrow =data$df , Xcol = data$df, dcor = TRUE,
                         R=0, rvar = TRUE, cvar = TRUE, nscan = iter, symmetric = FALSE, 
                         model = family,burn = burnin, odens = 10, plot = FALSE) 
    
    colnames(amen_2_R_0_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")
    int <- mcmc_intervals_data(as.mcmc(amen_2_R_0_a$BETA), prob_outer = .95, point_est = "mean")[,-1]
    
    
  }
  
  # run amen with indep beta and discrete U update step
  
  if(clusterNum == 4){
    Model <- "MCMCLAMBDAU"
    
    cluster_symm     <- matrix(specClust_out(data$Y, data$df, tauAvg = TRUE, FALSE, 3), ncol = numGroup)
    
    
    amen_2_updateU_L_a <- myAME_estimateLambda(data$Y, Xdyad = NULL, Xrow = data$df, Xcol = data$df,
                                               dcor = TRUE,
                                               R=3, rvar = TRUE, cvar = TRUE, nscan =iter, symmetric = FALSE,
                                               model = family,burn = burnin, initialU = cluster_symm,
                                               initialV = cluster_symm,odens = 10, plot = FALSE)
    
    colnames(amen_2_updateU_L_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")
    int <- mcmc_intervals_data(as.mcmc(amen_2_updateU_L_a$BETA), prob_outer = .95, point_est = "mean")[,-1]
    
    
    
  }
  
  
  if(clusterNum == 5){
    Model <- "SepAmen"
    X_r_1 <- data$X_r[1:50, 1:50,]
    X_r_2 <- data$X_r[51:100, 51:100,]
    X_r_3 <- data$X_r[101:150, 101:150,]
    
    X_c_1 <- data$X_c[1:50, 1:50,]
    X_c_2 <- data$X_c[51:100, 51:100,]
    X_c_3 <- data$X_c[101:150, 101:150,]
    
    Y_1 <- data$Y[1:50, 1:50]
    Y_2 <- data$Y[51:100, 51:100]
    Y_3 <- data$Y[101:150, 101:150]
    
    amen_2_R_0_1_a <-  ame(Y_1, Xdyad = NULL, Xrow =data$df[1:50,] , Xcol = data$df[1:50,], dcor = TRUE,
                           R=0, rvar = TRUE, cvar = TRUE, nscan = iter, symmetric = FALSE, model = "bin",burn = burnin, odens = 10,
                           plot = FALSE)
    
    amen_2_R_0_2_a <-  ame(Y_2, Xdyad = NULL, Xrow =data$df[51:100,] , Xcol = data$df[51:100,], dcor = TRUE,
                           R=0, rvar = TRUE, cvar = TRUE, nscan = iter, symmetric = FALSE, model = "bin",burn = burnin, odens = 10,
                           plot = FALSE)
    
    amen_2_R_0_3_a <-  ame(Y_3, Xdyad = NULL, Xrow =data$df[101:150,] , Xcol = data$df[101:150,], dcor = TRUE,
                           R=0, rvar = TRUE, cvar = TRUE, nscan = iter, symmetric = FALSE, model = "bin",burn = burnin, odens = 10,
                           plot = FALSE)
    
    
    # comparing amen models
    colnames(amen_2_R_0_1_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")
    colnames(amen_2_R_0_2_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")
    colnames(amen_2_R_0_3_a$BETA) <- c("intercept", "X1.r", "X2.r", "X1.c", "X2.c")
    
    
    
    
    int1 <- mcmc_intervals_data(as.mcmc(amen_2_R_0_1_a$BETA), prob_outer = .95, point_est = "mean")[,-1]
    int2 <- mcmc_intervals_data(as.mcmc(amen_2_R_0_2_a$BETA), prob_outer = .95, point_est = "mean")[,-1]
    int3 <- mcmc_intervals_data(as.mcmc(amen_2_R_0_3_a$BETA), prob_outer = .95, point_est = "mean")[,-1]
    
    
    int <- as.data.frame(rbind(int1, int2, int3))
    
  }
  
  
  return(list(int = int, model = Model, data = dataName))
}


# Function to calc. probit likelihood for the basic adjusted amen models
probitLikelihood_myame <- function(U, V, lambda, XB, Z, Y,numU,rho,td,to){
  ULV <- U%*%lambda%*%t(V)
  n <- dim(Z)[1]
  meanMat <- XB + ULV
  
  lh <-0
  # what we care about...
  
  for(i in numU){
    for(j in 1:n){
      if(i !=j){
        if(Y[i,j]==1 && Y[j,i] == 1){
          
          
          tmp3 <-  pbivnorm(meanMat[i,j], y = meanMat[j,i], rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }   
          lh <- lh + log(tmp3)
          #lh<-lh+log(tmp1) +log(tmp2)
        }else if(Y[i,j]==1 && Y[j,i] == 0){
          
          
          tmp3 <-  pbivnorm(meanMat[i,j], y = -meanMat[j,i], -rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }          
          lh <- lh + log(tmp3)
          
          # lh<-lh+log(tmp1) +log(1-tmp2)
          
        }else if(Y[i,j]==0 && Y[j,i] == 1){
          
          tmp3 <-  pbivnorm(-meanMat[i,j], y = meanMat[j,i], -rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }           
          
          lh <- lh + log(tmp3)
          
        }else{
          tmp3 <-  pbivnorm(-meanMat[i,j], y = -meanMat[j,i], rho) 
          if(tmp3 < .000000001){tmp3 <- .0000000001}
          if((1-tmp3) < abs(.00001)){
            tmp3 <- .9999999999
          }   
          lh <- lh + log(tmp3)
          
          # lh <- lh+log(1-tmp1) +log(1-tmp2)
          
        }
        # print(tmp3);print(i);print(j); print(Y[i,j])
      }
    }
    
  }
  
  tmp1 <- pnorm(diag(meanMat)[numU], sd =sqrt(1 + rho))
  
  tmp1 <- ifelse(tmp1 == 0, .000000000001, tmp1)
  tmp1 <- ifelse(tmp1 == 1, .999999999999, tmp1)
  lh <-lh+ sum(log(1-tmp1))
  
  return(lh)
}


# Function to draw U, U/V
# Function to calculate dep X beta that is NOT all added together (returns array of XB)
#' @param lambda: current lambda (k x k matrix)
#' @param U: current U (n x k matrix)
#' @param V: current V (n x k matrix)
#' @param XB: current XB value
#' @param a,b: n vectors with current row/col random effects
#' @param Y: Observed network
#' @param UandV: TRUE/FALSE (are both estimated or just U)
#' @param counter: counts total iterations with burnin
#' @param MH: counts acceptance of MH step
#' @param X_r: n x n x p_dep array of row covar
#' @param X_c: n x n x p_dep array of col covar
#' @param X_d: n x n x p_dep_dyad array of dyadic covar
#' @param beta_r: current beta_r, p_dep x k matrix
#' @param beta_c: current beta_c, p_dep x k matrix
#' @param alpha: current alpha value
#' @param beta_dr: current beta_dr, p_dep_dyad x k matrix
#' @param beta_dc: current beta_dc, p_dep_dyad x k matrix
#' @param h_r_past: Previous calculation of h_r
#' @param h_c_past: Previous calculation of h_c
#' @param Z: current Z value (n x n matrix)
#' @param td, to: current de-correlation values
#' @param Sab: covar from a,b
#' @param rho: current rho value

drawUV_probitlLikelihood_amen_lambda <- function (lambda, U, V,  XB, a,b,Y, UandV,counter, MH,
                                                  h_r_past, h_c_past,Z,td,to,Sab,rho) {
  z <- dim(U)[2]
  
  #y.post    <- Y[lower.tri(Y, diag = FALSE)]
  
  y.post <-Y
  Z2 <- Z
  u.post <- U
  v.post <- V
  
  
  ####################################################################### Metropolis
  
  
  
  
  if(UandV == TRUE){
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    numV  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    oldV <-  apply(V[numV,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- prop_groupV <- c()
    currU       <- u.post[numU,]
    currV       <- v.post[numV,]
    curr_groupU <- curr_groupV <- c()
    currRowU <- currRowV    <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      currRowV   <- currV[k,]
      curr_groupV[k] <- which(currRowV == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1)
      temp_propV     <- sample(1:z,1)
      
      currRowU       <- currU[j,]
      currRowV       <- currV[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      curr_groupV[j] <- which(currRowV ==1)
      prop_groupV[j] <- temp_propV
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    temp_v.post                <- v.post
    
    temp_v.post[numV,]          <- rep(0,z)
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      temp_v.post[numV[h],prop_groupV[h]] <- 1
      
      
    }
    
    
    accept_prob      <-  probitLikelihood_myame(temp_u.post,temp_v.post, lambda, 
                                                XB+outer(a,b, "+") ,Z2,Y,c(numU,numV),rho,td,to)
    
    accept_prob_past <- probitLikelihood_myame(u.post,v.post, lambda,
                                               XB+ outer(a,b,"+"), Z2,Y,c(numU,numV),rho,td,to)
    
    
    
    
    runi <- runif(1)
    ratio_pro <- accept_prob -accept_prob_past
    
    print(ratio_pro)
    ###################################################3 accept reject
    if( ratio_pro > log(runi)) {
      
      ###################################################3 accept reject
      u.post <- temp_u.post
      v.post <- temp_v.post
      updateU <- TRUE
      swapPeep <- c(numU, numV)
      MH <- MH + 1
    }else{
      u.post <- u.post
      v.post <- v.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    counter <- counter + 1
  }else{
    numU  <- sample(1:dim(U)[1],2, replace = FALSE)
    
    oldU <-  apply(U[numU,], 1, function(x) which(x ==1))
    
    
    prop_groupU <- c()
    currU       <- u.post[numU,]
    curr_groupU <-  c()
    currRowU   <- c()
    
    for(k in 1:length(numU)){
      currRowU       <- currU[k,]
      curr_groupU[k] <- which(currRowU == 1)
      
      
    }
    
    
    ##### proposing membership
    for(j in 1:length(numU)){
      
      temp_propU     <- sample(1:z,1)
      
      currRowU       <- currU[j,]
      
      curr_groupU[j] <- which(currRowU ==1)
      prop_groupU[j] <- temp_propU
      
      
      
    }
    
    
    
    
    #################################### Proposal
    
    temp_u.post                <- u.post
    
    temp_u.post[numU,]          <- rep(0,z)
    
    
    
    for(h in 1:length(numU)){
      
      temp_u.post[numU[h],prop_groupU[h]] <- 1
      
      
    }
    
    
    
    ########################################### Compare
    
    
    accept_prob      <-  probitLikelihood_myame(temp_u.post,temp_u.post, lambda, 
                                                XB+outer(a,b, "+") ,Z2,Y,numU,rho,td,to)
    
    accept_prob_past <- probitLikelihood_myame(u.post,u.post, lambda,
                                               XB+ outer(a,b,"+"), Z2,Y,numU,rho,td,to)
    
    
    
    
    runi <- runif(1)
    ratio_pro <- accept_prob -accept_prob_past
    
    print(ratio_pro)
    ###################################################3 accept reject
    if( ratio_pro > log(runi)) {
      u.post <- temp_u.post
      v.post <- temp_u.post
      updateU <- TRUE
      swapPeep <- c(numU)
      MH <- MH + 1
      print("RATIOS")
      
    }else{
      u.post <- u.post
      v.post <- u.post
      updateU <- FALSE
      swapPeep <- NULL
    }
    
    
    U <- u.post
    V <- v.post
    
    
    counter <- counter + 1
  }
  ## pick random people to switch
  
  
  return(list(U = U, V = V, updateU = updateU,
              swapPeep = swapPeep, 
              MH =MH, counter = counter))
}


# Plotting Functions
# plot a matrix with given colors
matrixPlotter <- function(matrix, col2="darkblue", col1= "lightcyan2",title = "", alpha_num = .8){
  
  
  n <- dim(matrix)[1]
  p <-  ggplot(melt(matrix), aes(x = X1, y = X2, fill = as.factor(value))) +
    geom_raster() +ggtitle(title)+ xlab("") + ylab("")+
    theme_bw() +
    labs(fill = "")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
      aspect.ratio = 1,
      panel.grid = element_blank(), 
      plot.title = element_text(hjust = 0.5),
      # axis.line.x = element_line(),
      panel.border = element_rect(size=2,linetype="solid",color="gray"),
      plot.margin = rep(unit(0,"null"),4),
      panel.spacing = unit(0,"null")
      #plot.margin=unit(c(1,1,1,1), "cm")
      
      #panel.border=element_rect(fill=FALSE),axis.line.y = element_line()
    )+  scale_fill_manual(values = c(col1, col2))+
    scale_x_continuous(limits = c(0,n), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,n), expand = c(0, 0)) +guides(fill = FALSE)
  return(p) 
}


# additional clustering functions that take in mcmc estimates of beta and the true beta values.
# cluster based on the posterior mmeans of individual $U\tilde{beta}$ and then separate and 
# create intervals based on those estimated memberships as outlined in paper

ubeta_posterCI_cluster_some <- function(exp1,
                                        beta_c1, beta_c2, beta_r1, beta_r2,n,burnin = 10, indepDep = FALSE){
  
  
  if(burnin !=0){
    exp1$results$U<-exp1$results$U[-(1:burnin),]
    exp1$results$beta_rMat1 <-  exp1$results$beta_rMat1[-(1:burnin),]
    exp1$results$beta_cMat1 <-  exp1$results$beta_cMat1[-(1:burnin),]
    
  }
  
  
  
  currUBR1 <-  currUBR2 <-  currUBC1 <- currUBC2 <-membership <- matrix(0, nrow = dim(exp1$results$U)[1], ncol = n)
  for(l in 1:dim(exp1$results$U)[1]){
    currUBR2[l,] <-  c(matrix(exp1$results$U[l, ],nrow = n, ncol = 3) %*% exp1$results$beta_rMat1[l,])
    
    currUBC2[l,] <-  c(matrix(exp1$results$U[l, ],nrow = n, ncol = 3) %*% exp1$results$beta_cMat1[l,])
    membership[l, ] <- apply(matrix(exp1$results$U[l, ],nrow = n, ncol = 3), 1, function(x) which(x==1))
  }
  
  
  currUBC2_df <- as.data.frame( currUBC2)
  colnames(currUBC2_df) <- paste0("V", seq(1, n,1))
  intC2 <- as.data.frame(mcmc_intervals_data( currUBC2_df, prob_outer = .95, point_est = "mean"))
  classC2<- kmeans(intC2$m, centers=3)$cluster
  classC2 <- make_labels(classC2)
  
  
  
  
  currUBR2_df <- as.data.frame( currUBR2)
  colnames(currUBR2_df) <- paste0("V", seq(1, n,1))
  intR2 <- as.data.frame(mcmc_intervals_data( currUBR2_df, prob_outer = .95, point_est = "mean"))
  classR2<- kmeans(intR2$m, centers=3)$cluster
  classR2 <- make_labels(classR2)
  
  
  ## getoverall  means
  intR2_mean <- as.data.frame(rbind(colMeans(intR2[which(classR2 == "A"), c(5,7,9)]),
                                    colMeans(intR2[which(classR2 == "B"), c(5,7,9)]),
                                    colMeans(intR2[which(classR2 == "C"), c(5,7,9)])))
  intR2_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  
  
  intC2_mean <- as.data.frame(rbind(colMeans(intC2[which(classC2 == "A"), c(5,7,9)]),
                                    colMeans(intC2[which(classC2 == "B"), c(5,7,9)]),
                                    colMeans(intC2[which(classC2 == "C"), c(5,7,9)])))
  intC2_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  
  allInt <- rbind(intC2_mean, intR2_mean)
  
  allInt$covar <-c(rep("Column",3), rep("Row", 3))
  allInt$which <- c(rep(1,3), rep(2,3))
  trueLabels <- c(rep(1, n/3), rep(2, n/3), rep(3,n/3))
  #allInt$variable <- as.factor(c(rep(trueLabels, 4)))
  
  
  
  colorOrder <- function(m, trueBeta){
    
    tmp <- c()
    ordering <- c()
    for(j in 1:length(m)){
      for(i in 1:length(trueBeta)){
        tmp[i] <- abs(m[j]-trueBeta[i])
      }
      ordering[j] <- trueBeta[which.min(tmp)]
      
    }
    missing <- which(!(trueBeta %in% ordering)==TRUE)
    print(missing)
    if(length(missing) ==1){
      dups <- which(duplicated(ordering)==TRUE)
      ordering[dups[1]] <- trueBeta[missing[1]]
    }
    
    return(ordering)
  }
  
  beta_r2 <-colorOrder(intR2_mean$m, beta_r2)
  beta_c2 <-colorOrder(intC2_mean$m, beta_c2)
  
  
  # print(beta_c1)
  # print(beta_r1)
  # print(beta_c2)
  # print(beta_r2)
  
  
  dummy22 <- data.frame(X = c("Column", "Row"), one = c(beta_c2[1],beta_r2[1]), 
                        two = c(beta_c2[2],beta_r2[2]), three = c(beta_c2[3],beta_r2[3]) )
  dummy22$which <- c(2,2)
  
  dummyAll <-rbind(dummy22)
  theme_set(theme_grey())
  maxB <- max(c(allInt$hh)) + .05
  minB <- min(c(allInt$ll)) - .05
  
  #
  
  dummyMelt <- melt(dummyAll, id = c("X", "which"))
  
  if(indepDep == FALSE){
    levels(dummyMelt$which) <- c("Covariate 1", "Covariate 2")
    levels(allInt$which) <- c("Covariate 1", "Covariate 2")
    
  }else{
    levels(dummyMelt$which) <- c("Indep Covar", "Dep Covar")
    levels(allInt$which) <- c("Indep Covar", "Dep Covar")
    
    
  }
  
  levels(dummyMelt$variable) <- c("Community 1", "Community 2", "Community 3")
  
  gg1 <- ggplot(subset(allInt, covar == "Column"), aes(x= community, y= m))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh,col = community), 
                  position = position_dodge(width = .5))+ 
    #facet_grid(.~which) +
    geom_point(aes(col = community),size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta} U^T "))+
    ggtitle("95% CI for Column Covariates")+geom_hline(data = subset(dummyMelt, X == "Column"),
                                                       size = 1, 
                                                       aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    #scale_color_manual(values = c("maroon", "orange3", "cyan3"))+
    scale_color_manual(values = c("black", "dimgray", "gray"))+    theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                                                                                      axis.text.x=element_blank(),
                                                                                      axis.title.y = element_text(size = 16),
                                                                                      axis.text.y = element_text(size = 16),
                                                                                      title = element_text(size = 14),
                                                                                      strip.text = element_text(size = 16))
  
  
  
  gg2 <- ggplot(subset(allInt, covar == "Row"), aes(x= community, y= m))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh,col = community), 
                  position = position_dodge(width = .5))+ 
    
    #facet_grid(.~which) +
    geom_point(aes(col = community),size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta} U^T "))+
    ggtitle("95% CI for Row Covariates")+geom_hline(data = subset(dummyMelt, X == "Row"),
                                                    size = 1, 
                                                    aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    scale_color_manual(values = c("black", "dimgray", "gray"))+
    theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                       axis.text.x=element_blank(),
                       axis.title.y = element_text(size = 16),
                       axis.text.y = element_text(size = 16),
                       title = element_text(size = 14),
                       strip.text = element_text(size = 16))
  
  
  
  grid_arrange_shared_legend(gg1,gg2)
  return(list(allInt= allInt, dummyMelt = dummyMelt))
  
}


# make labeels for classification
make_labels<-function(class1){
  first_lab <- class1[1]
  second_lab <- class1[-which(class1==first_lab)][1]
  third_lab <- class1[-which(class1%in% c(first_lab,second_lab))][1]
  
  class_new <- ifelse(class1 == first_lab, "A",
                      ifelse(class1 == second_lab, "B", "C"))
  return(class_new)
}


# ubeta (cluster our Ubetas) to get posterior credible interval
# mainn plotting function for simulations
# @param: exp1 is model input coming from our mmodel after running amen master
# @param: beta_c1,..., beta_r2 are all true values (since this is for plotting simulation)
# @paramm: n is number of indiviiduals
# @param: burnin: specify number of iterations for burnin
# @param indepDep set to false, specifies if there were comm. indep variables (for future adding
# into this funnction, sep. plotting function has been created for indep variables)
# @param prob_choice: posterior intervals level
ubeta_posterCI_cluster_may2020<- function(exp1,
                                          beta_c1, beta_c2, beta_r1, beta_r2,n,burnin = 0, 
                                          indepDep = FALSE,
                                          prob_choice = .95){
  
  
  if(burnin !=0){
    exp1$U<-exp1$U[-c(1:burnin),]
    exp1$beta_rMat1 <- exp1$beta_rMat1[-c(1:burnin),]
    exp1$beta_cMat1 <- exp1$beta_cMat1[-c(1:burnin),]
    exp1$beta_rMat2 <- exp1$beta_rMat2[-c(1:burnin),]
    exp1$beta_cMat2 <- exp1$beta_cMat2[-c(1:burnin),]
    
  }
  
  
  currUBR1 <-  currUBR2 <-  currUBC1 <- currUBC2 <-membership <- matrix(0, nrow = dim(exp1$U)[1], ncol = n)
  for(l in 1:(dim(exp1$U)[1])){
    currUBR1[l,] <- c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_rMat1[l,])
    currUBR2[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_rMat2[l,])
    
    currUBC1[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_cMat1[l,])
    currUBC2[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_cMat2[l,])
    membership[l, ] <- apply(matrix(exp1$U[l, ],nrow = n, ncol = 3), 1, function(x) which(x==1))
  }
  
  currUBC1_df <- as.data.frame( currUBC1)
  colnames(currUBC1_df) <- paste0("V", seq(1, n,1))
  intC1 <- as.data.frame(mcmc_intervals_data( currUBC1_df, prob_outer = prob_choice, point_est = "mean"))
  classC1<- kmeans(intC1$m, centers=3)$cluster
  classC1 <- make_labels(classC1)
  
  currUBC2_df <- as.data.frame( currUBC2)
  colnames(currUBC2_df) <- paste0("V", seq(1, n,1))
  intC2 <- as.data.frame(mcmc_intervals_data( currUBC2_df, prob_outer = prob_choice, point_est = "mean"))
  classC2<- kmeans(intC2$m, centers=3)$cluster
  classC2 <- make_labels(classC2)
  
  
  currUBR1_df <- as.data.frame( currUBR1)
  colnames(currUBR1_df) <- paste0("V", seq(1, n,1))
  intR1 <- as.data.frame(mcmc_intervals_data( currUBR1_df, prob_outer = prob_choice, point_est = "mean"))
  
  classR1<- kmeans(intR1$m, centers=3)$cluster
  classR1 <- make_labels(classR1)
  
  currUBR2_df <- as.data.frame( currUBR2)
  colnames(currUBR2_df) <- paste0("V", seq(1, n,1))
  intR2 <- as.data.frame(mcmc_intervals_data( currUBR2_df, prob_outer = prob_choice, point_est = "mean"))
  classR2<- kmeans(intR2$m, centers=3)$cluster
  classR2 <- make_labels(classR2)
  
  
  ## getoverall  means
  intR2_mean <- as.data.frame(rbind(colMeans(intR2[which(classR2 == "A"), c(5,7,9)]),
                                    colMeans(intR2[which(classR2 == "B"), c(5,7,9)]),
                                    colMeans(intR2[which(classR2 == "C"), c(5,7,9)])))
  intR2_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  
  
  intR1_mean <- as.data.frame(rbind(colMeans(intR1[which(classR1 == "A"), c(5,7,9)]),
                                    colMeans(intR1[which(classR1 == "B"), c(5,7,9)]),
                                    colMeans(intR1[which(classR1 == "C"), c(5,7,9)])))
  intR1_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  
  intC1_mean <- as.data.frame(rbind(colMeans(intC1[which(classC1 == "A"), c(5,7,9)]),
                                    colMeans(intC1[which(classC1 == "B"), c(5,7,9)]),
                                    colMeans(intC1[which(classC1 == "C"), c(5,7,9)])))
  intC1_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  intC2_mean <- as.data.frame(rbind(colMeans(intC2[which(classC2 == "A"), c(5,7,9)]),
                                    colMeans(intC2[which(classC2 == "B"), c(5,7,9)]),
                                    colMeans(intC2[which(classC2 == "C"), c(5,7,9)])))
  intC2_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  
  allInt <- rbind(intC1_mean, intR1_mean, intC2_mean, intR2_mean)
  
  allInt$covar <-c(rep("Column",3), rep("Row", 3), rep("Column",3), rep("Row", 3))
  allInt$which <- c(rep(1,3*2), rep(2,3*2))
  trueLabels <- c(rep(1, n/3), rep(2, n/3), rep(3,n/3))
  #allInt$variable <- as.factor(c(rep(trueLabels, 4)))
  
  
  
  colorOrder <- function(m, trueBeta){
    
    tmp <- c()
    ordering <- c()
    for(j in 1:length(m)){
      for(i in 1:length(trueBeta)){
        tmp[i] <- abs(m[j]-trueBeta[i])
      }
      ordering[j] <- trueBeta[which.min(tmp)]
      
    }
    missing <- which(!(trueBeta %in% ordering)==TRUE)
    print(missing)
    if(length(missing) ==1){
      dups <- which(duplicated(ordering)==TRUE)
      ordering[dups[1]] <- trueBeta[missing[1]]
    }
    
    return(ordering)
  }
  
  beta_r1 <-colorOrder(intR1_mean$m, beta_r1)
  beta_c1 <-colorOrder(intC1_mean$m, beta_c1)
  beta_r2 <-colorOrder(intR2_mean$m, beta_r2)
  beta_c2 <-colorOrder(intC2_mean$m, beta_c2)
  
  
  # print(beta_c1)
  # print(beta_r1)
  # print(beta_c2)
  # print(beta_r2)
  
  
  dummy2 <- data.frame(X = c("Column", "Row"), one = c(beta_c1[1],beta_r1[1]), 
                       two = c(beta_c1[2],beta_r1[2]), three = c(beta_c1[3],beta_r1[3]) )
  dummy22 <- data.frame(X = c("Column", "Row"), one = c(beta_c2[1],beta_r2[1]), 
                        two = c(beta_c2[2],beta_r2[2]), three = c(beta_c2[3],beta_r2[3]) )
  dummy2$which  <- c(1,1)
  dummy22$which <- c(2,2)
  
  dummyAll <-rbind(dummy2, dummy22)
  theme_set(theme_grey())
  maxB <- max(c(allInt$hh)) + .05
  minB <- min(c(allInt$ll)) - .05
  
  #
  
  dummyMelt <- melt(dummyAll, id = c("X", "which"))
  
  if(indepDep == FALSE){
    levels(dummyMelt$which) <- c("Covariate 1", "Covariate 2")
    levels(allInt$which) <- c("Covariate 1", "Covariate 2")
    
  }else{
    levels(dummyMelt$which) <- c("Indep Covar", "Dep Covar")
    levels(allInt$which) <- c("Indep Covar", "Dep Covar")
    
    
  }
  
  levels(dummyMelt$variable) <- c("Community 1", "Community 2", "Community 3")
  
  gg1 <- ggplot(subset(allInt, covar == "Column"), aes(x= community, y= m))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh,col = community), 
                  position = position_dodge(width = .5))+  facet_grid(.~which) +
    geom_point(aes(col = community),size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta} U^T "))+
    # ggtitle("95% CI for Column Covariates")+
    geom_hline(data = subset(dummyMelt, X == "Column"),
               size = 1, 
               aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    scale_color_manual(values = c("black", "dimgray", "gray"))+     theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                                                                                       axis.text.x=element_blank(),
                                                                                       axis.title.y = element_text(size = 16),
                                                                                       axis.text.y = element_text(size = 16),
                                                                                       title = element_text(size = 14),
                                                                                       strip.text = element_text(size = 16))
  
  
  
  gg2 <- ggplot(subset(allInt, covar == "Row"), aes(x= community, y= m))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh,col = community), 
                  position = position_dodge(width = .5))+  facet_grid(.~which) +
    geom_point(aes(col = community),size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta} U^T "))+
    #ggtitle("95% CI for Row Covariates")+
    geom_hline(data = subset(dummyMelt, X == "Row"),
               size = 1, 
               aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    scale_color_manual(values = c("black", "dimgray", "gray"))+theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                                                                                  axis.text.x=element_blank(),
                                                                                  axis.title.y = element_text(size = 16),
                                                                                  axis.text.y = element_text(size = 16),
                                                                                  title = element_text(size = 14),
                                                                                  strip.text = element_text(size = 16))
  
  
  
  grid_arrange_shared_legend(gg1,gg2)
  return(allInt)
}


# more plotting functions
# for model with indep and dep covariates.

# @param depResults: model output from amen_master, this function will plot for a specified
# covariate that was specific for our plots, not generalized for other use
# @param trueBeta_1R,trueBeta_1C: true value of indep. covariate beta
# @paramm burnin: amount of burnin to allow
plotBeta_giveAll_uni_INDEP <- function(depResults, trueBeta_1R,trueBeta_1C, burnin = 40){
  
  
  
  beta_cMat1 <- depResults$beta_cMatI1
  beta_rMat1 <- depResults$beta_rMatI1
  
  
  estCol1 <- as.mcmc(beta_cMat1[-c(1:burnin),])
  
  estRow1 <- as.mcmc(beta_rMat1[-c(1:burnin), ])
  
  intC1 <- as.data.frame(mcmc_intervals_data(as.matrix(estCol1), prob_outer = .95, point_est = "mean"))
  intR1 <- as.data.frame(mcmc_intervals_data(as.matrix(estRow1), prob_outer = .95, point_est = "mean"))
  
  intC1$covar <- rep("Column", 1)
  intR1$covar <- rep("Row", 1)
  
  intC1$which <-rep(1,1)
  intR1$which <-rep(1,1)
  
  
  allInt <- rbind(intC1, intR1)
  
  
  beta_r1 <- trueBeta_1R
  beta_c1 <- trueBeta_1C
  
  
  ##################### get indep
  
  theme_set(theme_grey())
  maxB <- max(c(allInt$hh)) + .05
  minB <- min(c(allInt$ll)) - .05
  
  
  allInt$which <- as.factor(allInt$which)
  levels(allInt$which) <- c("Covariate 1")
  
  
  truth <- as.data.frame(rbind(beta_r1, beta_c1))
  truth$covar <- c("Row", "Column")
  dummyMelt   <- melt(truth)
  
  
  
  ggplot(allInt, aes(x= covar, y= m, col =covar))+ theme_bw()+
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh), 
                  position = position_dodge(width = .5))+ 
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\beta"))+
    ggtitle("95% CI for Independent Covariate")+geom_hline(data = dummyMelt,size = 2, 
                                                           aes(yintercept = value,
                                                               col = covar)) + ylim(c(minB,maxB))+
    guides(col = FALSE)+ scale_color_manual(values = c("turquoise2", "black"))+
    theme( axis.title.x =  element_text(size = 16),axis.text.x =  element_text(size = 16),
           axis.title.y = element_text(size = 16),
           axis.text.y = element_text(size = 16),
           title = element_text(size = 16),
           strip.text = element_text(size = 16)) 
  
  
  
  return(allInt)
  
}


# specific for cluster
plotBeta_giveAll_uniCLUSTER <- function( beta_cMat1,beta_rMat1,
                                         trueBeta_1R,trueBeta_1C){
  
  
  
  estCol1 <- as.mcmc(beta_cMat1)
  
  colnames(estCol1) <- c("1", "2", "3")
  
  
  estRow1 <- as.mcmc(beta_rMat1)
  colnames(estRow1) <- c("1", "2", "3")

  intC1 <- as.data.frame(mcmc_intervals_data(estCol1, prob_outer = .95, point_est = "mean"))
  intR1 <- as.data.frame(mcmc_intervals_data(estRow1, prob_outer = .95, point_est = "mean"))
  
  intC1$covar <- rep("Column", 3)
  intR1$covar <- rep("Row", 3)
  
  intC1$which <-rep(1,3)
  intR1$which <-rep(1,3)
  
  
  allInt <- rbind(intC1, intR1)
  
  
  beta_r1 <- trueBeta_1R
  beta_c1 <- trueBeta_1C
  
  colorOrder <- function(m, trueBeta){
    
    tmp <- c()
    ordering <- c()
    for(j in 1:length(m)){
      for(i in 1:length(trueBeta)){
        tmp[i] <- abs(m[j]-trueBeta[i])
      }
      ordering[j] <- trueBeta[which.min(tmp)]
    }
    
    return(ordering)
  }
  
  beta_r1 <-colorOrder(intR1$m, beta_r1)
  beta_c1 <-colorOrder(intC1$m, beta_c1)
  
  dummy2 <- data.frame(X = c("Column", "Row"), one = c(beta_c1[1],beta_r1[1]), 
                       two = c(beta_c1[2],beta_r1[2]), three = c(beta_c1[3],beta_r1[3]) )
  dummy2$which  <- c(1,1)
  
  dummyAll <-rbind(dummy2)
  
  ##################### get indep
  
  theme_set(theme_grey())
  maxB <- max(c(allInt$hh)) + .05
  minB <- min(c(allInt$ll)) - .05
  
  
  allInt$parameter <- as.factor(allInt$parameter)
  levels(allInt$parameter) <- c("Community 1", "Community 2", "Community 3")
  allInt$which <- as.factor(allInt$which)
  # levels(allInt$which) <- c("Covariate 1")
  levels(allInt$which) <- c("Dep Covar")
  
  dummyMelt <- melt(dummyAll, id = c("X", "which"))
  levels(dummyMelt$variable) <- levels(allInt$parameter)
  dummyMelt$which <- as.factor(dummyMelt$which)
  # levels(dummyMelt$which) <- c("Covariate 1")
  levels(dummyMelt$which) <- c("Dep Covar")
  
  
  
  
  colnames(allInt)[1] <- "Community"
  gg1 <- ggplot(subset(allInt, covar == "Column"), aes(x= covar, y= m, col = Community, group = Community))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh, group = Community), 
                  position = position_dodge(width = .5))+ theme_bw()+
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    ggtitle("95% CI for Dep. Column Covariate")+geom_hline(data = subset(dummyMelt, X == "Column"),size = 1, 
                                                           aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    facet_grid(.~which) + 
    #scale_color_manual(values = c("maroon", "orange3", "cyan3"))+ theme( axis.title.x =  element_text(size = 14),
    scale_color_manual(values = c("orange3", "cyan3", "maroon"))+ theme( axis.title.x =  element_text(size = 14),
                                                                         axis.text.x=element_blank(),
                                                                         axis.title.y = element_text(size = 14),
                                                                         axis.text.y = element_text(size = 14),
                                                                         title = element_text(size = 14),
                                                                         # strip.text = element_text(size = 14))
                                                                         strip.text = element_blank())
  
  
  
  
  gg2 <- ggplot(subset(allInt, covar == "Row"), aes(x= covar, y= m, col = Community, group = Community))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5,
                  aes(ymin= ll, ymax = hh, group = Community), 
                  position = position_dodge(width = .5))+ guides(col=guide_legend(title="Estimate"))+theme_bw()+
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    ggtitle("95% CI for Dep. Row Covariate")+geom_hline(data = subset(dummyMelt, X == "Row"),size = 1, 
                                                        aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    facet_grid(.~which)+ 
    scale_color_manual(values = c("orange3", "cyan3", "maroon"))+ theme( axis.title.x =  element_text(size = 14),
                                                                         axis.text.x=element_blank(),
                                                                         axis.title.y = element_text(size = 14),
                                                                         axis.text.y = element_text(size = 14),
                                                                         title = element_text(size = 14),
                                                                         #  strip.text = element_text(size = 14))
                                                                         strip.text = element_blank())
  
  
  
  grid_arrange_shared_legend(gg1,gg2)
  
  
  
  
}

# plot ci for 2 dep covariates when u is not being estimated (do not have to worry
# about label switching in this case)
# @param beta_cMat1... beta_rMat2: MCMC output from model (either the initialize at U and stay, 
# or if you specified true community labels)
# @param   trueBeta_1C,trueBeta_2C,trueBeta_1R,trueBeta_2R: true values
# @param seet indepDep to FALSE for this functio
plot_2Dep <- function(beta_cMat1,beta_cMat2,beta_rMat1,beta_rMat2,
                      trueBeta_1C,trueBeta_2C,trueBeta_1R,trueBeta_2R,indepDep = FALSE){
  
  
  estCol1 <- as.mcmc(beta_cMat1)
  estCol2 <- as.mcmc(beta_cMat2)
  
  colnames(estCol1) <- c("1", "2", "3")
  
  colnames(estCol2) <- c("1", "2", "3")
  
  estRow1 <- as.mcmc(beta_rMat1)
  colnames(estRow1) <- c("1", "2", "3")
  
  
  estRow2 <- as.mcmc(beta_rMat2)
  colnames(estRow2) <- c("1", "2", "3")
  
  
  
  ## here we get label switching
  
  intC1 <- as.data.frame(mcmc_intervals_data(estCol1, prob_outer = .95, point_est = "mean"))
  intR1 <- as.data.frame(mcmc_intervals_data(estRow1, prob_outer = .95, point_est = "mean"))
  intC2 <- as.data.frame(mcmc_intervals_data(estCol2, prob_outer = .95, point_est = "mean"))
  intR2 <- as.data.frame(mcmc_intervals_data(estRow2, prob_outer = .95, point_est = "mean"))
  
  
  intC1$covar <- rep("Column", 3)
  intR1$covar <- rep("Row", 3)
  intC2$covar <- rep("Column", 3)
  intR2$covar <- rep("Row", 3)
  
  intC1$which <-rep(1,3)
  intC2$which <-rep(2,3)
  intR1$which <-rep(1,3)
  intR2$which <-rep(2,3)
  
  
  allInt <- rbind(intC1, intR1, intC2, intR2)
  
  
  beta_r1 <- trueBeta_1R
  beta_c1 <- trueBeta_1C
  beta_r2 <- trueBeta_2R
  beta_c2 <- trueBeta_2C
  
  
  print(beta_c1)
  print(beta_r1)
  print(beta_c2)
  print(beta_r2)
  
  
  colorOrder <- function(m, trueBeta){
    
    tmp <- c()
    ordering <- c()
    for(j in 1:length(m)){
      for(i in 1:length(trueBeta)){
        tmp[i] <- abs(m[j]-trueBeta[i])
      }
      ordering[j] <- trueBeta[which.min(tmp)]
      
    }
    missing <- which(!(trueBeta %in% ordering)==TRUE)
    print(missing)
    if(length(missing) ==1){
      dups <- which(duplicated(ordering)==TRUE)
      ordering[dups[1]] <- trueBeta[missing[1]]
    }
    
    return(ordering)
  }
  
  beta_r1 <-colorOrder(intR1$m, beta_r1)
  beta_c1 <-colorOrder(intC1$m, beta_c1)
  beta_r2 <-colorOrder(intR2$m, beta_r2)
  beta_c2 <-colorOrder(intC2$m, beta_c2)
  
  
  print(beta_c1)
  print(beta_r1)
  print(beta_c2)
  print(beta_r2)
  
  
  dummy2 <- data.frame(X = c("Column", "Row"), one = c(beta_c1[1],beta_r1[1]), 
                       two = c(beta_c1[2],beta_r1[2]), three = c(beta_c1[3],beta_r1[3]) )
  dummy22 <- data.frame(X = c("Column", "Row"), one = c(beta_c2[1],beta_r2[1]), 
                        two = c(beta_c2[2],beta_r2[2]), three = c(beta_c2[3],beta_r2[3]) )
  dummy2$which  <- c(1,1)
  dummy22$which <- c(2,2)
  
  dummyAll <-rbind(dummy2, dummy22)
  theme_set(theme_grey())
  maxB <- max(c(allInt$hh)) + .05
  minB <- min(c(allInt$ll)) - .05
  
  
  allInt$parameter <- as.factor(allInt$parameter)
  levels(allInt$parameter) <- c("Community 1", "Community 2", "Community 3")
  allInt$which <- as.factor(allInt$which)
  #
  
  dummyMelt <- melt(dummyAll, id = c("X", "which"))
  levels(dummyMelt$variable) <- levels(allInt$parameter)
  dummyMelt$which <- as.factor(dummyMelt$which)
  
  if(indepDep == FALSE){
    levels(dummyMelt$which) <- c("Covariate 1", "Covariate 2")
    levels(allInt$which) <- c("Covariate 1", "Covariate 2")
    
  }else{
    levels(dummyMelt$which) <- c("Indep Covar", "Dep Covar")
    levels(allInt$which) <- c("Indep Covar", "Dep Covar")
    
    
  }
  
  colnames(allInt)[1] <- "Community"
  gg1 <- ggplot(subset(allInt, covar == "Column"), aes(x= covar, y= m, col = Community, group = Community))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh, group = Community), 
                  position = position_dodge(width = .5))+ 
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    ggtitle("95% CI for Column Covariates")+geom_hline(data = subset(dummyMelt, X == "Column"),size = 1, 
                                                       aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    facet_grid(.~which) + scale_color_manual(values = c("maroon", "orange3", "cyan3"))+theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                                                                                                          axis.text.x=element_blank(),
                                                                                                          axis.title.y = element_text(size = 16),
                                                                                                          axis.text.y = element_text(size = 16),
                                                                                                          title = element_text(size = 14),
                                                                                                          strip.text = element_text(size = 16))
  
  
  
  
  gg2 <- ggplot(subset(allInt, covar == "Row"), aes(x= covar, y= m, col = Community, group = Community))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5,
                  aes(ymin= ll, ymax = hh, group = Community), 
                  position = position_dodge(width = .5))+ guides(col=guide_legend(title="Estimate"))+
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    ggtitle("95% CI for Row Covariates")+geom_hline(data = subset(dummyMelt, X == "Row"),size = 1, 
                                                    aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    facet_grid(.~which)+ scale_color_manual(values = c("maroon", "orange3", "cyan3"))+theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                                                                                                         axis.text.x=element_blank(),
                                                                                                         axis.title.y = element_text(size = 16),
                                                                                                         axis.text.y = element_text(size = 16),
                                                                                                         title = element_text(size = 14),
                                                                                                         strip.text = element_text(size = 16))
  
  grid_arrange_shared_legend(gg1,gg2)
  return(allInt)
  
}



# plot ci for standard amen results usinig ggplot, this generates plots for the standard model
# as well as the results from the partitioned amen models (we ran AMEN models on subnetworks
# partitioned by community)
# Reference the figure generating code for this to be a bit more intuitive
# @param allAMEN: output from 3 standard amen models that looked at R = 0, R = 3 , and our version
# of multiplicative effects 
# @param sepAMEN: output from the partitioned models
# @param beta_c1, beta_r1, beta_c2, beta_r2: true beta values
plotAMEN <- function(allAMEN,sepAMEN, beta_c1, beta_r1, beta_c2, beta_r2, allIndep){
  # beta_c1S <- beta_c1; beta_c2S <- beta_c2; beta_r1S <- beta_r1; beta_r2S <- beta_r2
  # 
  # 
  # colorOrder <- function(m, trueBeta){
  #   
  #   tmp <- c()
  #   ordering <- c()
  #   for(j in 1:length(m)){
  #     for(i in 1:length(trueBeta)){
  #       tmp[i] <- abs(m[j]-trueBeta[i])
  #     }
  #     ordering[j] <- trueBeta[which.min(tmp)]
  #   }
  #   
  #   return(ordering)
  # }
  # 
  # 
  # beta_c1 <- colorOrder(sepAMEN$m[c(1,5,9)], beta_c1)
  # beta_r1 <- colorOrder(sepAMEN$m[c(2,6,10)], beta_r1)
  # beta_c2 <- colorOrder(sepAMEN$m[c(3,7,11)], beta_c2)
  # beta_r2 <- colorOrder(sepAMEN$m[c(4,8,12)], beta_r2)
  # 
  # 
  # 
  
  allAMEN$coef <- as.factor(allAMEN$coef)
  ## get avg beta value
  dummy2 <- data.frame(X = c("Column", "Row"), one = c(beta_c1[1],beta_r1[1]), 
                       two = c(beta_c1[2],beta_r1[2]), three = c(beta_c1[3],beta_r1[3]) )
  dummy22 <- data.frame(X = c("Column", "Row"), one = c(beta_c2[1],beta_r2[1]), 
                        two = c(beta_c2[2],beta_r2[2]), three = c(beta_c2[3],beta_r2[3]) )
  dummy2$which  <- c(1,1)
  dummy22$which <- c(2,2)
  
  dummyAll <-rbind(dummy2, dummy22)
  dummyAll$mean <- rowMeans(dummyAll[,2:4])
  dummyAll_amen <- dummyAll[,c(1,5,6)]
  dummyAll_indiv <- dummyAll[,c(1,2,3,4,5)]
  
  dummyMelt <- melt(dummyAll_amen, id = c("X", "which"))
  dummyMelt$which <- as.factor(dummyMelt$which)
  if(allIndep==TRUE){
    # levels(dummyMelt$which) <- c("Covariate 1", "Covariate 2")
    dummyMelt$coef <- c("Covariate 1(Col)","Covariate 1(Row)", "Covariate 2(Col)","Covariate 2(Row)")
    
  }else{
    dummyMelt$coef <- c("Indep Covar (Col)","Indep Covar (Row)", "Dep Covar (Col)","Dep Covar (Row)")
    
  }
  
  
  dummyMelt_ind <- melt(dummyAll_indiv, id = c("X", "which"))
  dummyMelt_ind$which <- as.factor(dummyMelt_ind$which)
  
  if(allIndep==TRUE){
    # levels(dummyMelt$which) <- c("Covariate 1", "Covariate 2")
    dummyMelt_ind$coef <- rep(c("Covariate 1(Col)","Covariate 1(Row)", "Covariate 2(Col)","Covariate 2(Row)"),3)
    
    levels(allAMEN$coef) <- c("Covariate 1(Col)","Covariate 1(Row)", "Covariate 2(Col)","Covariate 2(Row)")
    
  }else{
    dummyMelt_ind$coef <- rep(c("Indep Covar (Col)","Indep Covar (Row)", "Dep Covar (Col)","Dep Covar (Row)"),3)
    
  }
  
  
  
  
  
  g1 <- ggplot(allAMEN, aes(x = model, y = m, col = model)) +facet_grid(.~coef)+
    geom_point(aes(), size = 3.7)+ 
    geom_errorbar(aes(ymin = ll, ymax =hh), width = .5,size = 2)+
    ggtitle(TeX("Posterior $95\\%$ CI for AMEN, $\\beta$")) + 
    ylab(TeX("$\\beta ")) + xlab("Model")+ labs(col = "Model")+scale_color_manual("",  values =  c("cyan3", "violetred4", "dimgray"), 
                                                                                  labels = unname(c( TeX("  MCMC: Update $\\Lambda$"),"AMEN: R = 0", "AMEN: R = 3"
                                                                                  )))
  
  
  amenPlot <- g1+geom_hline(data = dummyMelt,size = 1.5, 
                            aes(yintercept = value, group = coef))+theme_bw()+ theme(legend.position = "bottom",
                                                                                     legend.title = element_blank(),
                                                                                     legend.text = element_text(size=16),
                                                                                     legend.box = "horizontal", 
                                                                                     axis.title.x =  element_text(size = 16),
                                                                                     axis.text.x=element_blank(),
                                                                                     axis.title.y = element_text(size = 20),
                                                                                     axis.text.y = element_text(size = 16),
                                                                                     title = element_text(size = 16),
                                                                                     strip.text = element_text(size = 16))
  model2 <- c(
    rep("Community 1", 4), 
    rep("Community 2", 4),
    rep("Community 3", 4))
  
  
  
  colnames(dummyMelt_ind)[3] <- "model"
  
  dummyMelt_ind$model <-model2
  
  sepAMEN$model <- c(
    rep("Community 1", 4), 
    rep("Community 2", 4),
    rep("Community 3", 4))
  
  if(allIndep==TRUE){
    sepAMEN$coef <- rep(c("Covariate 1(Col)","Covariate 1(Row)", "Covariate 2(Col)","Covariate 2(Row)"),3)
    
  }else{
    sepAMEN$coef <- rep(c("Indep Covar (Col)","Indep Covar (Row)", "Dep Covar (Col)","Dep Covar (Row)"),3)
    
  }
  
  
  sepAMEN <- as.data.frame(sepAMEN)
  
  
  g1 <- ggplot(sepAMEN, aes(x = model, y = m, col = model)) +facet_grid(.~coef)+
    geom_point(aes(), size = 3.7)+ 
    geom_errorbar(aes(ymin = ll, ymax =hh), width = .5,size = 2)+
    ggtitle(TeX("Posterior $95\\%$ CI for Separate Community AMEN, $\\beta$")) + 
    ylab(TeX("$\\beta$")) + xlab("Community Model")+
    scale_color_manual(values = c("cyan3", "magenta3", "hotpink2"))+ labs(col = "Community")
  
  sepPlot <- g1+geom_hline(data = dummyMelt_ind,size = 1.5, 
                           aes(yintercept = value, group = coef))+theme_bw()+ theme(legend.position = "bottom",
                                                                                    legend.title = element_blank(),
                                                                                    legend.text = element_text(size=16),
                                                                                    legend.box = "horizontal", 
                                                                                    axis.title.x =  element_text(size = 16),
                                                                                    axis.text.x=element_blank(),
                                                                                    axis.title.y = element_text(size = 20),
                                                                                    axis.text.y = element_text(size = 16),
                                                                                    title = element_text(size = 16),
                                                                                    strip.text = element_text(size = 16))
  return(list(amenPlot = amenPlot, sepPlot = sepPlot, dummyMelt_sep= dummyMelt_ind,
              dummyMelt_ind = dummyMelt, sepAMEN = sepAMEN,
              allAMEN = allAMEN))
}






# creates a data frame for U \tilde{beta}$ for up to 2 covariates

# @param: exp1 is model input coming from our mmodel after running amen master
# @param: beta_c1,..., beta_r2 are all true values (since this is for plotting simulation)
# @paramm: n is number of indiviiduals
# @param: burnin: specify number of iterations for burnin
# @param indepDep: Set to FALSE for now
# @param prob_choice: desired credible itnerval percent
getUB_df <- function(exp1,
                     beta_c1, beta_c2, beta_r1, beta_r2,n,burnin = 0, indepDep = FALSE,
                     prob_choice = .95){
  
  
  
  if(burnin !=0){
    exp1$U<-exp1$U[-c(1:burnin),]
    exp1$beta_rMat1 <- exp1$beta_rMat1[-c(1:burnin),]
    exp1$beta_cMat1 <- exp1$beta_cMat1[-c(1:burnin),]
    exp1$beta_rMat2 <- exp1$beta_rMat2[-c(1:burnin),]
    exp1$beta_cMat2 <- exp1$beta_cMat2[-c(1:burnin),]
    
  }
  
  currUBR1 <-  currUBR2 <-  currUBC1 <- currUBC2 <-membership <- matrix(0, nrow = dim(exp1$U)[1], ncol = n)
  for(l in 1:(dim(exp1$U)[1])){
    #currUBR1[l,] <- c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_rMat1[l,])
    #currUBR2[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_rMat2[l,])
    
    currUBC1[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_cMat1[l,])
    #currUBC2[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_cMat2[l,])
    #membership[l, ] <- apply(matrix(exp1$U[l, ],nrow = n, ncol = 3), 1, function(x) which(x==1))
  }
  
  currUBC1_df <- as.data.frame( currUBC1)
  colnames(currUBC1_df) <- paste0("V", seq(1, n,1))
  intC1 <- as.data.frame(mcmc_intervals_data( currUBC1_df, prob_outer = prob_choice, point_est = "mean"))
  #classC1<- kmeans(intC1$m, centers=3)$cluster
  #classC1 <- make_labels(classC1)
  
  return(currUBC1_df)
}




# function for intermediate saves output from the cluster, helper function if ever using intermediate
# saves
ubeta_posterCI_cluster_intermediate_saves <- function(exp1,
                                                      beta_c1, beta_c2, beta_r1, beta_r2,n,burnin = 0, indepDep = FALSE,
                                                      prob_choice = .95){
  
  
  narows <- which(is.na(exp1$beta_rMat1[,1]))
  if(length(narows)!=0){
    exp1$U<-exp1$U[-narows,]
    
  }
  exp1$U<-exp1$U[-c(1:burnin),]
  exp1$beta_rMat1 <- exp1$beta_rMat1[-c(1:burnin),]
  exp1$beta_cMat1 <- exp1$beta_cMat1[-c(1:burnin),]
  exp1$beta_rMat2 <- exp1$beta_rMat2[-c(1:burnin),]
  exp1$beta_cMat2 <- exp1$beta_cMat2[-c(1:burnin),]
  
  
  currUBR1 <-  currUBR2 <-  currUBC1 <- currUBC2 <-membership <- matrix(0, nrow = dim(exp1$U)[1], ncol = n)
  for(l in 1:(dim(exp1$U)[1])){
    currUBR1[l,] <- c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_rMat1[l,])
    currUBR2[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_rMat2[l,])
    
    currUBC1[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_cMat1[l,])
    currUBC2[l,] <-  c(matrix(exp1$U[l, ],nrow = n, ncol = 3) %*% exp1$beta_cMat2[l,])
    membership[l, ] <- apply(matrix(exp1$U[l, ],nrow = n, ncol = 3), 1, function(x) which(x==1))
  }
  
  currUBC1_df <- as.data.frame( currUBC1)
  colnames(currUBC1_df) <- paste0("V", seq(1, n,1))
  intC1 <- as.data.frame(mcmc_intervals_data( currUBC1_df, prob_outer = prob_choice, point_est = "mean"))
  classC1<- kmeans(intC1$m, centers=3)$cluster
  classC1 <- make_labels(classC1)
  
  currUBC2_df <- as.data.frame( currUBC2)
  colnames(currUBC2_df) <- paste0("V", seq(1, n,1))
  intC2 <- as.data.frame(mcmc_intervals_data( currUBC2_df, prob_outer = prob_choice, point_est = "mean"))
  classC2<- kmeans(intC2$m, centers=3)$cluster
  classC2 <- make_labels(classC2)
  
  
  currUBR1_df <- as.data.frame( currUBR1)
  colnames(currUBR1_df) <- paste0("V", seq(1, n,1))
  intR1 <- as.data.frame(mcmc_intervals_data( currUBR1_df, prob_outer = prob_choice, point_est = "mean"))
  
  classR1<- kmeans(intR1$m, centers=3)$cluster
  classR1 <- make_labels(classR1)
  
  currUBR2_df <- as.data.frame( currUBR2)
  colnames(currUBR2_df) <- paste0("V", seq(1, n,1))
  intR2 <- as.data.frame(mcmc_intervals_data( currUBR2_df, prob_outer = prob_choice, point_est = "mean"))
  classR2<- kmeans(intR2$m, centers=3)$cluster
  classR2 <- make_labels(classR2)
  
  
  ## getoverall  means
  intR2_mean <- as.data.frame(rbind(colMeans(intR2[which(classR2 == "A"), c(5,7,9)]),
                                    colMeans(intR2[which(classR2 == "B"), c(5,7,9)]),
                                    colMeans(intR2[which(classR2 == "C"), c(5,7,9)])))
  intR2_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  
  
  intR1_mean <- as.data.frame(rbind(colMeans(intR1[which(classR1 == "A"), c(5,7,9)]),
                                    colMeans(intR1[which(classR1 == "B"), c(5,7,9)]),
                                    colMeans(intR1[which(classR1 == "C"), c(5,7,9)])))
  intR1_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  
  intC1_mean <- as.data.frame(rbind(colMeans(intC1[which(classC1 == "A"), c(5,7,9)]),
                                    colMeans(intC1[which(classC1 == "B"), c(5,7,9)]),
                                    colMeans(intC1[which(classC1 == "C"), c(5,7,9)])))
  intC1_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  intC2_mean <- as.data.frame(rbind(colMeans(intC2[which(classC2 == "A"), c(5,7,9)]),
                                    colMeans(intC2[which(classC2 == "B"), c(5,7,9)]),
                                    colMeans(intC2[which(classC2 == "C"), c(5,7,9)])))
  intC2_mean$community <- c("Community 1", "Community 2", "Community 3")
  
  
  allInt <- rbind(intC1_mean, intR1_mean, intC2_mean, intR2_mean)
  
  allInt$covar <-c(rep("Column",3), rep("Row", 3), rep("Column",3), rep("Row", 3))
  allInt$which <- c(rep(1,3*2), rep(2,3*2))
  trueLabels <- c(rep(1, n/3), rep(2, n/3), rep(3,n/3))
  #allInt$variable <- as.factor(c(rep(trueLabels, 4)))
  
  
  
  colorOrder <- function(m, trueBeta){
    
    tmp <- c()
    ordering <- c()
    for(j in 1:length(m)){
      for(i in 1:length(trueBeta)){
        tmp[i] <- abs(m[j]-trueBeta[i])
      }
      ordering[j] <- trueBeta[which.min(tmp)]
      
    }
    missing <- which(!(trueBeta %in% ordering)==TRUE)
    print(missing)
    if(length(missing) ==1){
      dups <- which(duplicated(ordering)==TRUE)
      ordering[dups[1]] <- trueBeta[missing[1]]
    }
    
    return(ordering)
  }
  
  beta_r1 <-colorOrder(intR1_mean$m, beta_r1)
  beta_c1 <-colorOrder(intC1_mean$m, beta_c1)
  beta_r2 <-colorOrder(intR2_mean$m, beta_r2)
  beta_c2 <-colorOrder(intC2_mean$m, beta_c2)
  
  
  # print(beta_c1)
  # print(beta_r1)
  # print(beta_c2)
  # print(beta_r2)
  
  
  dummy2 <- data.frame(X = c("Column", "Row"), one = c(beta_c1[1],beta_r1[1]), 
                       two = c(beta_c1[2],beta_r1[2]), three = c(beta_c1[3],beta_r1[3]) )
  dummy22 <- data.frame(X = c("Column", "Row"), one = c(beta_c2[1],beta_r2[1]), 
                        two = c(beta_c2[2],beta_r2[2]), three = c(beta_c2[3],beta_r2[3]) )
  dummy2$which  <- c(1,1)
  dummy22$which <- c(2,2)
  
  dummyAll <-rbind(dummy2, dummy22)
  theme_set(theme_grey())
  maxB <- max(c(allInt$hh)) + .05
  minB <- min(c(allInt$ll)) - .05
  
  #
  
  dummyMelt <- melt(dummyAll, id = c("X", "which"))
  
  if(indepDep == FALSE){
    levels(dummyMelt$which) <- c("Covariate 1", "Covariate 2")
    levels(allInt$which) <- c("Covariate 1", "Covariate 2")
    
  }else{
    levels(dummyMelt$which) <- c("Indep Covar", "Dep Covar")
    levels(allInt$which) <- c("Indep Covar", "Dep Covar")
    
    
  }
  
  levels(dummyMelt$variable) <- c("Community 1", "Community 2", "Community 3")
  
  gg1 <- ggplot(subset(allInt, covar == "Column"), aes(x= community, y= m))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh,col = community), 
                  position = position_dodge(width = .5))+  facet_grid(.~which) +
    geom_point(aes(col = community),size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta} U^T "))+
    # ggtitle("95% CI for Column Covariates")+
    geom_hline(data = subset(dummyMelt, X == "Column"),
               size = 1, 
               aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    scale_color_manual(values = c("black", "dimgray", "gray"))+     theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                                                                                       axis.text.x=element_blank(),
                                                                                       axis.title.y = element_text(size = 16),
                                                                                       axis.text.y = element_text(size = 16),
                                                                                       title = element_text(size = 14),
                                                                                       strip.text = element_text(size = 16))
  
  
  
  gg2 <- ggplot(subset(allInt, covar == "Row"), aes(x= community, y= m))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh,col = community), 
                  position = position_dodge(width = .5))+  facet_grid(.~which) +
    geom_point(aes(col = community),size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta} U^T "))+
    #ggtitle("95% CI for Row Covariates")+
    geom_hline(data = subset(dummyMelt, X == "Row"),
               size = 1, 
               aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    scale_color_manual(values = c("black", "dimgray", "gray"))+theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                                                                                  axis.text.x=element_blank(),
                                                                                  axis.title.y = element_text(size = 16),
                                                                                  axis.text.y = element_text(size = 16),
                                                                                  title = element_text(size = 14),
                                                                                  strip.text = element_text(size = 16))
  
  
  
  grid_arrange_shared_legend(gg1,gg2)
  return(allInt)
}




##### Additional functionns
# plot 2 dep covariates that have more than 2 communities
plot_2DepManyGROUPS <- function(beta_cMat1,beta_cMat2,beta_rMat1,beta_rMat2,
                                trueBeta_1C,trueBeta_2C,trueBeta_1R,trueBeta_2R,indepDep = FALSE,numGroup){
  
  seqToTry <- c(1,4,5,3,6,2)
  beta_cMat1 <- beta_cMat1[,seqToTry]
  
  beta_rMat1 <- beta_rMat1[,seqToTry]
  beta_rMat2 <- beta_rMat2[,seqToTry]
  beta_cMat2 <- beta_cMat2[,seqToTry]
  
  
  estCol1 <- as.mcmc(beta_cMat1)
  estCol2 <- as.mcmc(beta_cMat2)
  
  estRow1 <- as.mcmc(beta_rMat1)
  
  
  estRow2 <- as.mcmc(beta_rMat2)
  
  
  colnames(estCol1) <- colnames(estRow1) <- colnames(estCol2) <- colnames(estRow2) <- as.character(seq(1, numGroup))
  
  
  
  
  
  ## here we get label switching
  
  intC1 <- as.data.frame(mcmc_intervals_data(estCol1, prob_outer = .95, point_est = "mean"))
  intR1 <- as.data.frame(mcmc_intervals_data(estRow1, prob_outer = .95, point_est = "mean"))
  intC2 <- as.data.frame(mcmc_intervals_data(estCol2, prob_outer = .95, point_est = "mean"))
  intR2 <- as.data.frame(mcmc_intervals_data(estRow2, prob_outer = .95, point_est = "mean"))
  
  
  intC1$covar <- rep("Column", numGroup)
  intR1$covar <- rep("Row", numGroup)
  intC2$covar <- rep("Column", numGroup)
  intR2$covar <- rep("Row", numGroup)
  
  intC1$which <-rep(1,numGroup)
  intC2$which <-rep(2,numGroup)
  intR1$which <-rep(1,numGroup)
  intR2$which <-rep(2,numGroup)
  
  
  allInt <- rbind(intC1, intR1, intC2, intR2)
  
  
  beta_r1 <- trueBeta_1R
  beta_c1 <- trueBeta_1C
  beta_r2 <- trueBeta_2R
  beta_c2 <- trueBeta_2C
  
  
  
  
  
  colorOrder <- function(m, trueBeta){
    
    tmp <- c()
    ordering <- c()
    for(j in 1:length(m)){
      for(i in 1:length(trueBeta)){
        tmp[i] <- abs(m[j]-trueBeta[i])
      }
      ordering[j] <- trueBeta[which.min(tmp)]
      
    }
    missing <- which(!(trueBeta %in% ordering)==TRUE)
    print(missing)
    if(length(missing) ==1){
      dups <- which(duplicated(ordering)==TRUE)
      ordering[dups[1]] <- trueBeta[missing[1]]
    }
    
    return(ordering)
  }
  trueGroups <- length(beta_r1)
  minSeq <- c(1,2,3)
  beta_r1 <-colorOrder(intR1$m[minSeq], beta_r1)
  beta_c1 <-colorOrder(intC1$m[minSeq], beta_c1)
  beta_r2 <-colorOrder(intR2$m[minSeq], beta_r2)
  beta_c2 <-colorOrder(intC2$m[minSeq], beta_c2)
  
  
  print(beta_c1)
  print(beta_r1)
  print(beta_c2)
  print(beta_r2)
  
  
  dummy2 <- data.frame(X = c("Column", "Row"), 
                       one = c(beta_c1[1],beta_r1[1]), 
                       two = c(beta_c1[2],beta_r1[2]), 
                       three = c(beta_c1[3],beta_r1[3]) )
  
  dummy22 <- data.frame(X = c("Column", "Row"), 
                        one = c(beta_c2[1],beta_r2[1]), 
                        two = c(beta_c2[2],beta_r2[2]), 
                        three = c(beta_c2[3],beta_r2[3]) )
  dummy2$which  <- c(1,1)
  dummy22$which <- c(2,2)
  
  dummyAll <-rbind(dummy2, dummy22)
  theme_set(theme_grey())
  maxB <- max(c(allInt$hh)) + .05
  minB <- min(c(allInt$ll)) - .05
  
  
  allInt$parameter <- as.factor(allInt$parameter)
  #Slevels(allInt$parameter) <- c("Community 1", "Community 2", "Community 3")
  allInt$which <- as.factor(allInt$which)
  #
  
  dummyMelt <- melt(dummyAll, id = c("X", "which"))
  levels(dummyMelt$variable) <- levels(allInt$parameter)
  dummyMelt$which <- as.factor(dummyMelt$which)
  
  if(indepDep == FALSE){
    levels(dummyMelt$which) <- c("Covariate 1", "Covariate 2")
    levels(allInt$which) <- c("Covariate 1", "Covariate 2")
    
  }else{
    levels(dummyMelt$which) <- c("Indep Covar", "Dep Covar")
    levels(allInt$which) <- c("Indep Covar", "Dep Covar")
    
    
  }
  
  colnames(allInt)[1] <- "Community"
  gg1 <- ggplot(subset(allInt, covar == "Column"), aes(x= covar, y= m, col = Community, group = Community))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh, group = Community), 
                  position = position_dodge(width = .5))+ 
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    ggtitle("95% CI for Column Covariates")+geom_hline(data = subset(dummyMelt, X == "Column"),size = 1, 
                                                       aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    facet_grid(.~which) +
    # + scale_color_manual(values = c( "orange3","maroon","blue",
    #                                                       "cyan3",
    #                                                       "forestgreen", "gray2"))+
    theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                       axis.text.x=element_blank(),
                       axis.title.y = element_text(size = 16),
                       axis.text.y = element_text(size = 16),
                       title = element_text(size = 14),
                       strip.text = element_text(size = 16))
  
  
  
  
  gg2 <- ggplot(subset(allInt, covar == "Row"), aes(x= covar, y= m, col = Community, group = Community))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5,
                  aes(ymin= ll, ymax = hh, group = Community), 
                  position = position_dodge(width = .5))+ guides(col=guide_legend(title="Estimate"))+
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    ggtitle("95% CI for Row Covariates")+geom_hline(data = subset(dummyMelt, X == "Row"),size = 1, 
                                                    aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    facet_grid(.~which)+ theme_bw()+ theme( axis.title.x =  element_text(size = 16),
                                            axis.text.x=element_blank(),
                                            axis.title.y = element_text(size = 16),
                                            axis.text.y = element_text(size = 16),
                                            title = element_text(size = 14),
                                            strip.text = element_text(size = 16))
  
  
  
  
  grid_arrange_shared_legend(gg1,gg2)
  
}




plotAllCovar_dep <- function(results, numGroup,burnin=0){
  
  for (i in seq_along(results)){
    colnames(results[[i]]) <- as.character(seq(1, numGroup))
    
  } 
  
  allInt <- c()
  for(i in 1:length(results)){
    coeffname <- names(results)[i]
    
    nam <- paste("est", coeffname, sep = "")
    assign(nam, as.mcmc(results[[i]][-c(1:burnin),]))
    
    
    
    
    
    
    tmpDat <- as.data.frame(mcmc_intervals_data(eval(parse(text=nam)), prob_outer = .95, point_est = "mean"))
    nam2 <- paste("intMCMC", coeffname, sep = "")  
    assign(nam2, tmpDat)
    
    
    allInt <- rbind(allInt, get(nam2))
    
    
  }
  
  covar <- c()
  for(i in 1:length(results)){
    covarT <- rep(names(results)[i], numGroup)
    covar <- c(covar, covarT)
  }
  
  
  
  allInt$covar <- covar
  
  allInt$parameter <- as.factor(allInt$parameter)
  levels(allInt$parameter) <- c("Community 1", "Community 2", "Community 3")
  
  ggplot(allInt, aes(x= covar, y= m, col =parameter, group = parameter))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh, group = parameter), 
                  position = position_dodge(width = .5))+ 
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    geom_hline(aes(yintercept = 0)) 
  #guides(col = FALSE)+
  
  
}

# additioal cluster function that take in mcmc estimates of beta and the true beta values
plotBeta_giveAll_uniCLUSTER <- function( beta_cMat1,beta_rMat1,
                                         trueBeta_1R,trueBeta_1C, burnin = 1000){
  
  
  if(burnin !=0){
    estCol1 <- as.mcmc(beta_cMat1[-c(1:burnin),])
    
    colnames(estCol1) <- c("1", "2", "3")
    
    
    estRow1 <- as.mcmc(beta_rMat1[-c(1:burnin), ])
    colnames(estRow1) <- c("1", "2", "3")
  }else{
    estCol1 <- as.mcmc(beta_cMat1)
    
    colnames(estCol1) <- c("1", "2", "3")
    
    
    estRow1 <- as.mcmc(beta_rMat1)
    colnames(estRow1) <- c("1", "2", "3")
  }
  
  
  
  
  ## here we get label switching
  
  intC1 <- as.data.frame(mcmc_intervals_data(estCol1, prob_outer = .95, point_est = "mean"))
  intR1 <- as.data.frame(mcmc_intervals_data(estRow1, prob_outer = .95, point_est = "mean"))
  
  intC1$covar <- rep("Column", 3)
  intR1$covar <- rep("Row", 3)
  
  intC1$which <-rep(1,3)
  intR1$which <-rep(1,3)
  
  
  allInt <- rbind(intC1, intR1)
  
  
  beta_r1 <- trueBeta_1R
  beta_c1 <- trueBeta_1C
  
  colorOrder <- function(m, trueBeta){
    
    tmp <- c()
    ordering <- c()
    for(j in 1:length(m)){
      for(i in 1:length(trueBeta)){
        tmp[i] <- abs(m[j]-trueBeta[i])
      }
      ordering[j] <- trueBeta[which.min(tmp)]
    }
    
    return(ordering)
  }
  
  beta_r1 <-colorOrder(intR1$m, beta_r1)
  beta_c1 <-colorOrder(intC1$m, beta_c1)
  
  dummy2 <- data.frame(X = c("Column", "Row"), one = c(beta_c1[1],beta_r1[1]), 
                       two = c(beta_c1[2],beta_r1[2]), three = c(beta_c1[3],beta_r1[3]) )
  dummy2$which  <- c(1,1)
  
  dummyAll <-rbind(dummy2)
  
  ##################### get indep
  
  theme_set(theme_grey())
  maxB <- max(c(allInt$hh)) + .05
  minB <- min(c(allInt$ll)) - .05
  
  
  allInt$parameter <- as.factor(allInt$parameter)
  levels(allInt$parameter) <- c("Community 1", "Community 2", "Community 3")
  allInt$which <- as.factor(allInt$which)
  # levels(allInt$which) <- c("Covariate 1")
  levels(allInt$which) <- c("Dep Covar")
  
  dummyMelt <- melt(dummyAll, id = c("X", "which"))
  levels(dummyMelt$variable) <- levels(allInt$parameter)
  dummyMelt$which <- as.factor(dummyMelt$which)
  # levels(dummyMelt$which) <- c("Covariate 1")
  levels(dummyMelt$which) <- c("Dep Covar")
  
  
  
  
  colnames(allInt)[1] <- "Community"
  gg1 <- ggplot(subset(allInt, covar == "Column"), aes(x= covar, y= m, col = Community, group = Community))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5, 
                  aes(ymin= ll, ymax = hh, group = Community), 
                  position = position_dodge(width = .5))+ theme_bw()+
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    ggtitle("Mixed Simulation: 95% CI for Dep. Column Covariate")+geom_hline(data = subset(dummyMelt, X == "Column"),size = 1, 
                                                                             aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    #guides(col = FALSE)+
    facet_grid(.~which) + scale_color_manual(values = c("maroon", "orange3", "cyan3"))+ theme( axis.title.x =  element_text(size = 18),
                                                                                               axis.text.x=element_blank(),
                                                                                               axis.title.y = element_text(size = 18),
                                                                                               axis.text.y = element_text(size = 18),
                                                                                               title = element_text(size = 16),
                                                                                               # strip.text = element_text(size = 18))
                                                                                               strip.text = element_blank())
  
  
  
  
  gg2 <- ggplot(subset(allInt, covar == "Row"), aes(x= covar, y= m, col = Community, group = Community))+ 
    #+ facet_grid(Community~covar,scales = "free")+
    geom_errorbar(width = .8, size = 2.5,
                  aes(ymin= ll, ymax = hh, group = Community), 
                  position = position_dodge(width = .5))+ guides(col=guide_legend(title="Estimate"))+theme_bw()+
    geom_point(size = 6, position=position_dodge(width=.5))+xlab("")+ylab(TeX("$\\tilde{\\beta}_k "))+
    ggtitle("95% CI for Dep. Row Covariate")+geom_hline(data = subset(dummyMelt, X == "Row"),size = 1, 
                                                        aes(yintercept = value, col = variable)) + ylim(c(minB,maxB))+
    facet_grid(.~which)+ scale_color_manual(values = c("maroon", "orange3", "cyan3"))+ theme( axis.title.x =  element_text(size = 18),
                                                                                              axis.text.x=element_blank(),
                                                                                              axis.title.y = element_text(size = 18),
                                                                                              axis.text.y = element_text(size = 18),
                                                                                              title = element_text(size = 16),
                                                                                              #  strip.text = element_text(size = 18))
                                                                                              strip.text = element_blank())
  
  
  
  grid_arrange_shared_legend(gg1,gg2)
  
  return(allInt)
  
  
}





############### grid_arrange_shared_legend is from Shaun Jackmman
# public github
# Create a grid arrangement of ggplot2 plots with a shared legend
# Original code by Shaun Jackman (http://rpubs.com/sjackman/grid_arrange_shared_legend)
# Source: https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs


library(gridExtra)
library(grid)
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position, legend.title =element_blank(), legend.text = element_text(size=18)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}


