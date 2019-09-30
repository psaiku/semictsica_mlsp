rm(list=ls())
library(ica)
library(data.table)
library(truncnorm)
library(truncdist)

rsemicts <- function(n, pzero = 0.5, r.func = NA, cts.density = "truncnorm",
                     cts.param=list(a=0, b=Inf, mean = 0, sd=1)) {
  if(pzero==0 & is.na(r.func)) {
    stop("Pr(Y=0) cannot be 0!")
  }
  if(!is.na(r.func)) {
    r.func <- match.fun(r.func)
    pzero = r.func(cts.param)
  }
  u <- runif(n)
  x <- u
  x <- (u>pzero)*u

  idx <- which(x>0)

  fapp <- function(t, f, argss) {
    return(do.call(f, c(list(p=(u[t]-pzero)/(1-pzero)), argss)))
  }

  if(length(idx)>0) {
    f <- match.fun(paste0("q", cts.density))
    x[idx] <- apply(as.matrix(idx), 1, fapp, f, cts.param)
  }
  return(x)
}

semictsica <- function(X, maxit = 500, step_size = 0.01, eps = 0.0001, eps_x = 0.0001) {
  tmp <- dim(X)
  N <- tmp[1]
  T <- tmp[2]

  X_orig <- copy(X)
  tmp <- pre_processing(X)
  X <- tmp$X
  P <- tmp$P

  W = matrix(runif(N*N), nrow=N, ncol=N)
  W = W/P
  W_new = matrix(NA, nrow=N, ncol=N)

  Cost = matrix(0, maxit, 1)
  Grad= matrix(0, maxit, 1)

  for(iter in  1:maxit) {
    #print(sprintf("Iter: %s", iter))
    entropy <- matrix(0, N, 1)
    grad_mon <- matrix(0, N, 1)
    Cost[iter] <- -log(abs(det(W)))
    Wten = array(NA, c(N-1, N, N))
    for(kk in 1:N) {
      nout <- matrix(seq(1, N)[-c(kk)])
      Wten[,,kk] <- W[nout,]
    }
    for(n in 1:N) {
      #print(n)
      temp1 <- matrix(rnorm(N), nrow=N, ncol=1)
      temp2 <- Wten[,,n]
      h <- temp1 - t(temp2) %*% (solve((temp2 %*% t(temp2)), temp2)) %*% temp1
      w <- matrix(W[n,])
      y <- t(w) %*% X
      y <- (y>eps_x)*y
      delJ <- get_semicts_entropy_der(y)
      entropy[n] <- delJ$entropy
      phi <- delJ$phi
      grad = (X %*% phi)/T + h/as.numeric(t(h) %*% w)
      grad = grad - as.numeric(t(w) %*%grad) * w
      grad_mon[n] = norm_vec(grad)
      grad = grad / norm_vec(grad)
      w1 = w + step_size * grad
      w1 = w1 / norm_vec(w1)
      W_new[n,] <- t(w1)
    }

    Grad[iter] = sum(grad_mon)/N^2
    Cost[iter] = Cost[iter] + sum(entropy)
    #print(Grad[iter])
    #print(Cost[iter])

    if(iter > 1) {
      if(abs(Cost[iter]-Cost[iter-1]) < eps) {
        break
      } else {
        W <- W_new
      }
    } else {
      W <- W_new
    }
  }
  W <- W %*% P
  S_est <- W %*% X_orig
  return(list("W" = W, "Cost" = Cost, "Grad" = Grad, "S" = t(S_est)))
}

icaexp <- function(X, maxit = 500, step_size = 0.01, eps = 0.0001, eps_x = 0.0001) {
  tmp <- dim(X)
  N <- tmp[1]
  T <- tmp[2]
  
  X_orig <- copy(X)
  tmp <- pre_processing(X)
  X <- tmp$X
  P <- tmp$P
  
  W = matrix(runif(N*N), nrow=N, ncol=N)
  W = W/P
  W_new = matrix(NA, nrow=N, ncol=N)
  
  Cost = matrix(0, maxit, 1)
  Grad= matrix(0, maxit, 1)
  
  for(iter in  1:maxit) {
    #print(sprintf("Iter: %s", iter))
    entropy <- matrix(0, N, 1)
    grad_mon <- matrix(0, N, 1)
    Cost[iter] <- -log(abs(det(W)))
    Wten = array(NA, c(N-1, N, N))
    for(kk in 1:N) {
      nout <- matrix(seq(1, N)[-c(kk)])
      Wten[,,kk] <- W[nout,]
    }
    for(n in 1:N) {
      temp1 <- matrix(rnorm(N), nrow=N, ncol=1)
      temp2 <- Wten[,,n]
      h <- temp1 - t(temp2) %*% (solve((temp2 %*% t(temp2)), temp2)) %*% temp1
      w <- matrix(W[n,])
      y <- t(w) %*% X
      y <- (y>eps_x)*y
      delJ <- get_exp_der(y)
      entropy[n] <- delJ$entropy
      phi <- delJ$phi
      grad = (X %*% phi)/T + h/as.numeric(t(h) %*% w)
      grad = grad - as.numeric(t(w) %*%grad) * w
      grad_mon[n] = norm_vec(grad)
      grad = grad / norm_vec(grad)
      w1 = w + step_size * grad
      w1 = w1 / norm_vec(w1)
      W_new[n,] <- t(w1)
    }
    
    Grad[iter] = sum(grad_mon)/N^2
    Cost[iter] = Cost[iter] + sum(entropy)
    
    if(iter > 1) {
      if(abs(Cost[iter]-Cost[iter-1]) < eps) {
        break
      } else {
        W <- W_new
      }
    } else {
      W <- W_new
    }
  }
  W <- W %*% P
  S_est <- W %*% X_orig
  return(list("W" = W, "Cost" = Cost, "Grad" = Grad, "S" = t(S_est)))
}

ccorr1 <- function(i, T, S, S_hat, N) {
  crr1 <- c()
  x1 <- S[i,]
  for(j in 1:N) {
    x2 <- S_hat[j,]
    crr1[j] <- cor(x1, x2)
  }
  return(crr1)
}

inv_sqrtmH <- function(B) {
  eig <- eigen(B)
  d <- matrix(eig$values, ncol=1)
  d <- 1/sqrt(d)
  V <- eig$vectors
  A <- V %*% diag(c(d)) %*% t(V);
}

pre_processing <- function(X) {
  tmp <- dim(X)
  T <- tmp[2]
  Xmean <- matrix(rowMeans(X))
  X <- X - Xmean%*%matrix(1, nrow=1, ncol=T)
  R <- (X %*% t(X))/T
  P <- inv_sqrtmH(R)
  X <- P %*% X
  return(list("X" = X, "P" = P))
}

norm_vec <- function(x) sqrt(sum(x^2))

### We choose two measuring functions for the positive data: x and log(x). This choice, along
### with the proportion of zeroes in the data, will give us the two-part gamma distribution
### as the MaxEnt distribution.
### Given the vector y of size T, we first estimate the MaxEnt semi-cts gamma distr.
### We then calculate the entropy of this estimated distribution.
### Finally, we calculate the derivative vector.
get_semicts_entropy_der <- function(y) {
  T <- dim(t(y))[1]
  # Estimate the pdf using entropy maximization. Fix the constraints for semi-cts gamma
  entr <- 0
  der <- matrix(0, T, 1)
  idx <- which(y>0)
  gmma <- length(which(y==0))/length(y)
  mean_x <- mean(y[idx])
  mean_logx <- mean(log(y[idx]))

  f1 <- function(x) {
    mean_logx - log(x) - digamma(mean_x / x)
  }

  thet <- uniroot(f1, c(0.00001, 10000))$root
  k <- mean_x / thet
  entr_g <- k + log(thet) +log(gamma(k)) + (1-k)*digamma(k)
  if(gmma > 0) {
    entr <- -gmma * log(gmma) - (1-gmma)*log(1-gmma) + (1-gmma)*entr_g
  } else {
    entr <- entr_g
  }
  # Calculate the derivative
  f2 <- function(x) {
    if(x == 0) {
      return(0)
    } else {
      return(dgamma(x, shape = k, scale = thet) * (thet*k - thet - x)/(x*thet))
    }
  }
  der <- matrix(apply(y, 2, f2), ncol=1)
  return(list("entropy" = entr, "phi" = der))
}

get_exp_der <- function(y) {
  T <- dim(t(y))[1]
  # Approximate the pdf with an exponential distribution
  entr <- 0
  der <- matrix(0, T, 1)
  lambda_hat <- T/mean(y)
  entr <- 1-log(lambda_hat)
  der <- matrix(-lambda_hat, T, 1)
  return(list("entropy" = entr, "phi" = der))
}

simICA_scem <- function(X, S, N, T) {
  epsilon <- 0.0001
  corMat <- matrix(NA, N, N)
  ic <- semictsica(X)
  S_est <- t(ic$S)
  S_est <- (S_est > epsilon)*S_est
  for(i in 1:N) {
    corMat[i,] <- ccorr1(i, T, S, S_est, N)
  }
  return("dcor" = matrix(apply(corMat, 1, max), nrow=1, ncol=N))
}

simICA_exp <- function(X, S, N, T) {
  epsilon <- 0.0001
  corMat <- matrix(NA, N, N)
  ic <- icaexp(X)
  S_est <- t(ic$S)
  S_est <- (S_est > epsilon)*S_est
  for(i in 1:N) {
    corMat[i,] <- ccorr1(i, T, S, S_est, N)
  }
  return("dcor" = matrix(apply(corMat, 1, max), nrow=1, ncol=N))
}

simICA_fastica <- function(X, S, N, T) {
  epsilon <- 0.0001
  corMat <- matrix(NA, N, N)
  ic <- icafast(t(X), nc = N, maxit = 500, tol = 0.0001)
  S_est <- t(ic$S)
  S_est <- (S_est > epsilon)*S_est
  for(i in 1:N) {
    corMat[i,] <- ccorr1(i, T, S, S_est, N)
  }
  return("dcor" = matrix(apply(corMat, 1, max), nrow=1, ncol=N))
}

simICA_infomax <- function(X, S, N, T) {
  epsilon <- 0.0001
  corMat <- matrix(NA, N, N)
  ic <- icaimax(t(X), nc = N, maxit = 500, tol = 0.0001)
  S_est <- t(ic$S)
  S_est <- (S_est > epsilon)*S_est
  for(i in 1:N) {
    corMat[i,] <- ccorr1(i, T, S, S_est, N)
  }
  return("dcor" = matrix(apply(corMat, 1, max), nrow=1, ncol=N))
}

simICA_jade <- function(X, S, N, T) {
  epsilon <- 0.0001
  corMat <- matrix(NA, N, N)
  ic <- icajade(t(X), nc = N, maxit = 500, tol = 0.0001)
  S_est <- t(ic$S)
  S_est <- (S_est > epsilon)*S_est
  for(i in 1:N) {
    corMat[i,] <- ccorr1(i, T, S, S_est, N)
  }
  return("dcor" = matrix(apply(corMat, 1, max), nrow=1, ncol=N))
}

#### Arguments:
####  sim_scenarios - accepts "tpgamma", "tplognormal", and "mixed"
####  M - number of replications
####  N - number of latent sources
####  T - number of samples or realizations per source
####  barplot_file - full path to the png file to save the barplot
runSim <- function(sim_scenario = "tpgamma", M = 100, N = 5, T = 1000, barplot_file) {
  k <- 1
  dcorrMat_scem <- matrix(NA, M, N)
  dcorrMat_fastica <- matrix(NA, M, N)
  dcorrMat_infomax <- matrix(NA, M, N)
  dcorrMat_jade <- matrix(NA, M, N)
  print(sprintf("Running simulations for scenario %s with %d replications, %d sources, and %d samples per source", sim_scenario, M, N, T))

  while(k <= M) {
    print(sprintf("Iter: %d", k))
    set.seed(k)
    S <- matrix(NA, nrow=N, ncol=T)
    A <- matrix(runif(N*N), nrow=N, ncol=N)
    switch (sim_scenario,
            "tpgamma" = S <- matrix(rsemicts(N*T, pzero=0.6, cts.density="gamma", cts.param = list(shape = 1.5, rate = 1)), nrow=N, ncol=T),
            "tplognormal" = S <- matrix(rsemicts(N*T, pzero=0.6, cts.density="lnorm", cts.param = list(meanlog = 0, sdlog = 1)), nrow=N, ncol=T),
            "mixed" = {S[1,] <- matrix(rsemicts(T, pzero=0.6, cts.density="gamma", cts.param = list(shape = 1.5, rate = 2)), nrow=1, ncol=T)
            #S[2,] <- matrix(rsemicts(T, pzero=0.5, cts.density="gamma", cts.param = list(shape = 1, rate = 1)), nrow=1, ncol=T)
            S[2,] <- matrix(rgamma(T, shape = 2, rate = 1.5), nrow=1, ncol=T)
            S[3,] <- matrix(rsemicts(T, pzero=0.6, cts.density="lnorm", cts.param = list(meanlog = 0, sdlog = 1)), nrow=1, ncol=T)
            S[4,] <- matrix(rsemicts(T, pzero=0.5, cts.density="lnorm", cts.param = list(meanlog = 0.5, sdlog = 1.5)), nrow=1, ncol=T)
            #S[5,] <- matrix(rsemicts(T, pzero=0.5, cts.density="lnorm", cts.param = list(meanlog = 1, sdlog = 2)), nrow=1, ncol=T)}
            S[5,] <- matrix(rlnorm(T, meanlog=1, sdlog=2), nrow=1, ncol=T)}
    )
    X <- A%*%S
    dcorrMat_scem[k,] <- matrix(simICA_scem(X, S, N, T), nrow=1, ncol=N)
    dcorrMat_fastica[k,] <- matrix(simICA_fastica(X, S, N, T), nrow=1, ncol=N)
    dcorrMat_infomax[k,] <- matrix(simICA_infomax(X, S, N, T), nrow=1, ncol=N)
    dcorrMat_jade[k,] <- matrix(simICA_jade(X, S, N, T), nrow=1, ncol=N)
    #dcorrMat_jade[k,] <- matrix(simICA_exp(X, N, T), nrow=1, ncol=N)
    k <- k+1
  }

  ans <- matrix(NA, nrow=4, ncol=N)
  ans[1,] <- matrix(apply(dcorrMat_scem, 2, mean, na.rm=T), ncol=N)
  ans[2,] <- matrix(apply(dcorrMat_fastica, 2, mean, na.rm=T), ncol=N)
  ans[3,] <- matrix(apply(dcorrMat_infomax, 2, mean, na.rm=T), ncol=N)
  ans[4,] <- matrix(apply(dcorrMat_jade, 2, mean, na.rm=T), ncol=N)
  colnames(ans) <- seq(1:N)
  ylim <- c(0, 1)
  angle1 <- rep(c(45,45,135), length.out=4)
  angle2 <- rep(c(45,135,135), length.out=4)
  density1 <- seq(N,35,length.out=4)
  density2 <- seq(N,35,length.out=4)
  col <- 1 # rainbow(7)

  png(barplot_file)
  op <- par(mar=c(3,3,1,1))
  barplot(ans, ylim = ylim, beside=TRUE, col=col, angle=angle1, density=density1, ylab="Correlation", xlab="Sources")
  barplot(ans, add=TRUE, beside=TRUE, col=col, angle=angle2, density=density2, ylab="Correlation", xlab="Sources")
  legend("top", legend=c("ICA-SCEM", "FastICA", "ICA Info-Max", "ICA-JADE"), ncol=4, fill=TRUE, col=col, angle=angle1, density=density1)
  par(bg="transparent")
  legend("top", legend=c("ICA-SCEM", "FastICA", "ICA Info-Max", "ICA-JADE"), ncol=4, fill=TRUE, col=col, angle=angle2, density=density2)
  par(op)
  dev.off()

  return(ans)
}
