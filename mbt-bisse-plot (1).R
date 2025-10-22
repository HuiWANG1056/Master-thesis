############# lambda##############
## 参数
pars_sym <- c(lambda0=0.10, lambda1=0.10, mu0=0.03, mu1=0.03, q01=0.01, q10=0.01)
pars_asym_lambda <- c(lambda0=0.10, lambda1=0.20, mu0=0.03, mu1=0.03, q01=0.01, q10=0.01)

set.seed(1)
n <- 100
tips <- 300

sim_sym <- vector("list", n)
i <- 1L; attempts <- 0L
while (i <= n) {
  attempts <- attempts + 1L
  phy <- tree.bisse(pars_sym, max.taxa=tips, max.t=Inf, include.extinct=FALSE, x0=NA)
  if (!is.null(phy) && length(phy$tip.label)==tips) { sim_sym[[i]] <- phy; i <- i+1L }
  if (attempts %% 50L == 0L) message(sprintf("[sym] Kept %d/%d (attempts=%d)", i-1L, n, attempts))
}

sim_asym <- vector("list", n)
i <- 1L; attempts <- 0L
while (i <= n) {
  attempts <- attempts + 1L
  phy <- tree.bisse(pars_asym_lambda, max.taxa=tips, max.t=Inf, include.extinct=FALSE, x0=NA)
  if (!is.null(phy) && length(phy$tip.label)==tips) { sim_asym[[i]] <- phy; i <- i+1L }
  if (attempts %% 50L == 0L) message(sprintf("[asym] Kept %d/%d (attempts=%d)", i-1L, n, attempts))
}


par_names <- c("lambda0","lambda1","mu0","mu1","q01","q10")
est_sym  <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))
est_asym <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))

message("Fitting symmetric set...")
for (k in seq_len(n)) {
  phy <- sim_sym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  
  lik <- make.bisse.mbt(phy, states)              
  start <- starting.point.bisse(phy)          

  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_sym[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}
save.image("/Volumes/LaCie/hui wang/result/mbt/speciation/res_sym_notraitconst_nosurv")
message("Fitting asymmetric-λ set...")
for (k in seq_len(n)) {
  phy <- sim_asym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  
  lik <- make.bisse.mbt(phy, states)             
  start <- starting.point.bisse(phy)

  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_asym[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}
save.image("/Volumes/LaCie/hui wang/result/mbt/speciation/res_asym_notraitconst_nosurv")
save(est_asym,file="/Volumes/LaCie/hui wang/result/mbt/speciation/asym_notraitconst_nosurv")

est_sym  <- est_sym[complete.cases(est_sym), , drop=FALSE]
est_asym <- est_asym[complete.cases(est_asym), , drop=FALSE]


idx_sim <- which((sapply(1:100,function (i) sum(sim_sym[[i]]$tip.state))/300>=0.2)*(sapply(1:100,function (i) sum(sim_sym[[i]]$tip.state))/300<=0.8)==0)
color_sim <- rep("black",100)
color_sim[idx_sim] <- "red"

idx_asym <- which((sapply(1:100,function (i) sum(sim_asym[[i]]$tip.state))/300>=0.2)*(sapply(1:100,function (i) sum(sim_asym[[i]]$tip.state))/300<=0.8)==0)
color_asym <- rep("black",100)
color[idx_asym] <- "red"

par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[,1], est_sym[,2],
     xlim=c(0,0.4), ylim=c(0,0.4), asp=1,
     xlab=expression(estimated~lambda[0]),
     ylab=expression(estimated~lambda[1]),
     pch=1, cex=0.7, col=color)
points(est_asym[,1], est_asym[,2], pch=17, cex=0.7,col=color)

abline(v=0.10, lty=1, lwd=1)
abline(h=0.10, lty=1, lwd=1)
abline(h=0.20, lty=2, lwd=1)

m_sym  <- colMeans(est_sym[-idx,]);  m_asym <- colMeans(est_asym[-idx,])
points(m_sym[1],  m_sym[2],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym[1], m_asym[2], pch=24, bg="white", cex=1.6, lwd=1.2)

m_sym  <- colMeans(est_sym[idx,]);  m_asym <- colMeans(est_asym[idx,])
points(m_sym[1],  m_sym[2],  pch=21, bg="white", cex=1.6, lwd=1.2, col="red")
points(m_asym[1], m_asym[2], pch=24, bg="white", cex=1.6, lwd=1.2, col="red")


mtext("Speciation rates", side=3, adj=1, line=-1, cex=0.95)

## -------- Extinction rates--------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[,3], est_sym[,4],
     xlim=c(0,0.15), ylim=c(0,0.15), asp=1,
     xlab=expression(estimated~mu[0]),
     ylab=expression(estimated~mu[1]),
     pch=1, cex=0.7, col=color)
points(est_asym[,3], est_asym[,4], pch=17, cex=0.7,col=color)

abline(v=0.03, lty=1, lwd=1)  # mu0 (symmetric)
abline(h=0.03, lty=1, lwd=1)  # mu1 (symmetric)

m_sym  <- colMeans(est_sym[-idx,]);  m_asym <- colMeans(est_asym[-idx,])
points(m_sym["mu0"],  m_sym["mu1"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["mu0"], m_asym["mu1"], pch=24, bg="white", cex=1.6, lwd=1.2)

m_sym  <- colMeans(est_sym[idx,]);  m_asym <- colMeans(est_asym[idx,])
points(m_sym["mu0"],  m_sym["mu1"],  pch=21, bg="white", cex=1.6, lwd=1.2, col="red")
points(m_asym["mu0"], m_asym["mu1"], pch=24, bg="white", cex=1.6, lwd=1.2, col="red")


mtext("Extinction rates", side=3, adj=1, line=-1, cex=0.95)

## -------- State change rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[, "q01"], est_sym[, "q10"],
     xlim=c(0, 0.04), ylim=c(0, 0.04), asp=1,
     xlab=expression(estimated~q[01]),
     ylab=expression(estimated~q[10]),
     pch=1, cex=0.7, col="black")
points(est_asym[, "q01"], est_asym[, "q10"], pch=17, cex=0.7,col=color)

abline(v=0.01, lty=1, lwd=1)
abline(h=0.01, lty=1, lwd=1)

m_sym  <- colMeans(est_sym[-idx,]);  m_asym <- colMeans(est_asym[-idx,])
points(m_sym["q01"],  m_sym["q10"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["q01"], m_asym["q10"], pch=24, bg="white", cex=1.6, lwd=1.2)

m_sym  <- colMeans(est_sym[idx,]);  m_asym <- colMeans(est_asym[idx,])
points(m_sym["q01"],  m_sym["q10"],  pch=21, bg="white", cex=1.6, lwd=1.2, col="red")
points(m_asym["q01"], m_asym["q10"], pch=24, bg="white", cex=1.6, lwd=1.2, col="red")


mtext("State change rates", side=3, adj=1, line=-1, cex=0.95)


##BISSE
## -------- fit --------
par_names <- c("lambda0","lambda1","mu0","mu1","q01","q10")
est_sym_b  <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))
est_asym_b <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))

message("Fitting symmetric set...")
for (k in seq_len(n)) {
  phy <- sim_sym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  
  lik <- make.bisse(phy, states)              
  start <- starting.point.bisse(phy)          

  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_sym_b[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}

message("Fitting asymmetric-λ set...")
for (k in seq_len(n)) {
  phy <- sim_asym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  
  lik <- make.bisse(phy, states)             
  start <- starting.point.bisse(phy)

  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_asym_b[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}

est_sym_b  <- est_sym_b[complete.cases(est_sym_b), , drop=FALSE]
est_asym_b <- est_asym_b[complete.cases(est_asym_b), , drop=FALSE]

## -------- Speciation rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[,1], est_sym_b[,2],
     xlim=c(0,0.4), ylim=c(0,0.4), asp=1,
     xlab=expression(estimated~lambda[0]),
     ylab=expression(estimated~lambda[1]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[,1], est_asym_b[,2], pch=17, cex=0.7)

abline(v=0.10, lty=1, lwd=1)
abline(h=0.10, lty=1, lwd=1)
abline(h=0.20, lty=2, lwd=1)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym[1],  m_sym[2],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym[1], m_asym[2], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("Speciation rates", side=3, adj=1, line=-1, cex=0.95)

## -------- Extinction rates--------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[,3], est_sym_b[,4],
     xlim=c(0,0.15), ylim=c(0,0.15), asp=1,
     xlab=expression(estimated~mu[0]),
     ylab=expression(estimated~mu[1]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[,3], est_asym_b[,4], pch=17, cex=0.7)


abline(v=0.03, lty=1, lwd=1)  # mu0 (symmetric)
abline(h=0.03, lty=1, lwd=1)  # mu1 (symmetric)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym["mu0"],  m_sym["mu1"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["mu0"], m_asym["mu1"], pch=24, bg="white", cex=1.6, lwd=1.2)

mtext("Extinction rates", side=3, adj=1, line=-1, cex=0.95)

## -------- State change rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[, "q01"], est_sym_b[, "q10"],
     xlim=c(0, 0.04), ylim=c(0, 0.04), asp=1,
     xlab=expression(estimated~q[01]),
     ylab=expression(estimated~q[10]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[, "q01"], est_asym_b[, "q10"], pch=17, cex=0.7)


abline(v=0.01, lty=1, lwd=1)
abline(h=0.01, lty=1, lwd=1)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym["q01"],  m_sym["q10"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["q01"], m_asym["q10"], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("State change rates", side=3, adj=1, line=-1, cex=0.95)


#################mu#############
pars_sym <- c(lambda0=0.10, lambda1=0.10, mu0=0.03, mu1=0.03, q01=0.01, q10=0.01)
pars_asym_mu <- c(lambda0=0.10, lambda1=0.10, mu0=0.06, mu1=0.03, q01=0.01, q10=0.01)

pars_sym <- c(lambda0=0.10, lambda1=0.10, mu0=0.08, mu1=0.08, q01=0.01, q10=0.01)
pars_asym_mu <- c(lambda0=0.10, lambda1=0.10, mu0=0.10, mu1=0.08, q01=0.01, q10=0.01)

set.seed(1)
n <- 100
tips <- 300

sim_sym <- vector("list", n)
i <- 1L; attempts <- 0L
while (i <= n) {
  attempts <- attempts + 1L
  phy <- tree.bisse(pars_sym, max.taxa=tips, max.t=Inf, include.extinct=FALSE, x0=NA)
  if (!is.null(phy) && length(phy$tip.label)==tips) { sim_sym[[i]] <- phy; i <- i+1L }
  if (attempts %% 50L == 0L) message(sprintf("[sym] Kept %d/%d (attempts=%d)", i-1L, n, attempts))
}

sim_asym <- vector("list", n)
i <- 1L; attempts <- 0L
while (i <= n) {
  attempts <- attempts + 1L
  phy <- tree.bisse(pars_asym_mu, max.taxa=tips, max.t=Inf, include.extinct=FALSE, x0=NA)
  if (!is.null(phy) && length(phy$tip.label)==tips) { sim_asym[[i]] <- phy; i <- i+1L }
  if (attempts %% 50L == 0L) message(sprintf("[asym-μ] Kept %d/%d (attempts=%d)", i-1L, n, attempts))
}

par_names <- c("lambda0","lambda1","mu0","mu1","q01","q10")
est_sym  <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))
est_asym <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))

message("Fitting symmetric set...")
for (k in seq_len(n)) {
  phy <- sim_sym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  lik <- make.bisse.mbt(phy, states)
  start <- starting.point.bisse(phy)
  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_sym[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}
save.image("/Volumes/LaCie/hui wang/result/mbt/highextinction/res_sym_notraitconst_nosurv")
message("Fitting asymmetric-μ set...")
for (k in seq_len(n)) {
  phy <- sim_asym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  lik <- make.bisse.mbt(phy, states)
  start <- starting.point.bisse(phy)
  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_asym[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}
save.image("/Volumes/LaCie/hui wang/result/mbt/highextinction/res_asym_notraitconst_nosurv")
save(est_asym,file="/Volumes/LaCie/hui wang/result/mbt/highextinction/asym_notraitconst_nosurv")

est_sym  <- est_sym[complete.cases(est_sym), , drop=FALSE]
est_asym <- est_asym[complete.cases(est_asym), , drop=FALSE]

par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[,1]-est_sym[,3], est_sym[,2]-est_sym[,4],
     xlim=c(0,0.4), ylim=c(0,0.4), asp=1,
     xlab=expression(estimated~lambda[0]-mu[0]),
     ylab=expression(estimated~lambda[1]-mu[1]),
     pch=1, cex=0.7, col="black")
points(est_asym[,1]-est_asym[,3], est_asym[,2]-est_asym[,4], pch=17, cex=0.7)

abline(v=0.07, lty=1, lwd=1)  # lambda0-mu0 (symmetric)
abline(h=0.07, lty=1, lwd=1)  # lambda1-mu1 (symmetric)
abline(v=0.04, lty=2, lwd=1)  # lambda0-mu0 (asymmetric)

points(mean(est_sym[,1]-est_sym[,3]),  mean(est_sym[,2]-est_sym[,4]),  pch=21, bg="white", cex=1.6, lwd=1.2)
points(mean(est_asym[,1]-est_asym[,3]),  mean(est_sym[,2]-est_sym[,4]), pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("Diversification rates", side=3, adj=1, line=-1, cex=0.95)


## -------- Speciation rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[,1], est_sym[,2],
     xlim=c(0,0.4), ylim=c(0,0.4), asp=1,
     xlab=expression(estimated~lambda[0]),
     ylab=expression(estimated~lambda[1]),
     pch=1, cex=0.7, col="black")
points(est_asym[,1], est_asym[,2], pch=17, cex=0.7)

abline(v=0.10, lty=1, lwd=1)
abline(h=0.10, lty=1, lwd=1)

m_sym  <- colMeans(est_sym);  m_asym <- colMeans(est_asym)
points(m_sym[1],  m_sym[2],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym[1], m_asym[2], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("Speciation rates", side=3, adj=1, line=-1, cex=0.95)

## -------- Extinction rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[,3], est_sym[,4],
     xlim=c(0,0.15), ylim=c(0,0.15), asp=1,
     xlab=expression(estimated~mu[0]),
     ylab=expression(estimated~mu[1]),
     pch=1, cex=0.7, col="black")
points(est_asym[,3], est_asym[,4], pch=17, cex=0.7)


abline(v=0.03, lty=1, lwd=1)  # mu0 (symmetric)
abline(h=0.03, lty=1, lwd=1)  # mu1 (symmetric)
abline(v=0.06, lty=2, lwd=1)  # mu0 (asymmetric)

m_sym  <- colMeans(est_sym);  m_asym <- colMeans(est_asym)
points(m_sym["mu0"],  m_sym["mu1"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["mu0"], m_asym["mu1"], pch=24, bg="white", cex=1.6, lwd=1.2)

mtext("Extinction rates", side=3, adj=1, line=-1, cex=0.95)

## -------- State change rates  --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[, "q01"], est_sym[, "q10"],
     xlim=c(0, 0.04), ylim=c(0, 0.04), asp=1,
     xlab=expression(estimated~q[01]),
     ylab=expression(estimated~q[10]),
     pch=1, cex=0.7, col="black")
points(est_asym[, "q01"], est_asym[, "q10"], pch=17, cex=0.7)

abline(v=0.01, lty=1, lwd=1)
abline(h=0.01, lty=1, lwd=1)

m_sym  <- colMeans(est_sym);  m_asym <- colMeans(est_asym)
points(m_sym["q01"],  m_sym["q10"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["q01"], m_asym["q10"], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("State change rates", side=3, adj=1, line=-1, cex=0.95)

##BISSE
par_names <- c("lambda0","lambda1","mu0","mu1","q01","q10")
est_sym_b  <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))
est_asym_b <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))

message("Fitting symmetric set...")
for (k in seq_len(n)) {
  phy <- sim_sym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  lik <- make.bisse(phy, states)
  start <- starting.point.bisse(phy)
  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_sym_b[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}

message("Fitting asymmetric-μ set...")
for (k in seq_len(n)) {
  phy <- sim_asym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  lik <- make.bisse(phy, states)
  start <- starting.point.bisse(phy)
  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_asym_b[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}

est_sym_b  <- est_sym_b[complete.cases(est_sym_b), , drop=FALSE]
est_asym_b <- est_asym_b[complete.cases(est_asym_b), , drop=FALSE]

par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[,1]-est_sym_b[,3], est_sym_b[,2]-est_sym_b[,4],
     xlim=c(0,0.4), ylim=c(0,0.4), asp=1,
     xlab=expression(estimated~lambda[0]-mu[0]),
     ylab=expression(estimated~lambda[1]-mu[1]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[,1]-est_asym_b[,3], est_asym_b[,2]-est_asym_b[,4], pch=17, cex=0.7)

abline(v=0.07, lty=1, lwd=1)  # lambda0-mu0 (symmetric)
abline(h=0.07, lty=1, lwd=1)  # lambda1-mu1 (symmetric)
abline(v=0.04, lty=2, lwd=1)  # lambda0-mu0 (asymmetric)

points(mean(est_sym_b[,1]-est_sym_b[,3]),  mean(est_sym_b[,2]-est_sym_b[,4]),  pch=21, bg="white", cex=1.6, lwd=1.2)
points(mean(est_asym_b[,1]-est_asym_b[,3]),  mean(est_sym_b[,2]-est_sym_b[,4]), pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("Diversification rates", side=3, adj=1, line=-1, cex=0.95)

## -------- Speciation rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[,1], est_sym_b[,2],
     xlim=c(0,0.4), ylim=c(0,0.4), asp=1,
     xlab=expression(estimated~lambda[0]),
     ylab=expression(estimated~lambda[1]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[,1], est_asym_b[,2], pch=17, cex=0.7)

abline(v=0.10, lty=1, lwd=1)
abline(h=0.10, lty=1, lwd=1)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym[1],  m_sym[2],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym[1], m_asym[2], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("Speciation rates", side=3, adj=1, line=-1, cex=0.95)

## -------- Extinction rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[,3], est_sym_b[,4],
     xlim=c(0,0.3), ylim=c(0,0.3), asp=1,
     xlab=expression(estimated~mu[0]),
     ylab=expression(estimated~mu[1]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[,3], est_asym_b[,4], pch=17, cex=0.7)

abline(v=0.08, lty=1, lwd=1)  # mu0 (symmetric)
abline(h=0.08, lty=1, lwd=1)  # mu1 (symmetric)
abline(v=0.10, lty=2, lwd=1)  # mu0 (asymmetric)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym["mu0"],  m_sym["mu1"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["mu0"], m_asym["mu1"], pch=24, bg="white", cex=1.6, lwd=1.2)

mtext("Extinction rates", side=3, adj=1, line=-1, cex=0.95)

## -------- State change rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[, "q01"], est_sym_b[, "q10"],
     xlim=c(0, 0.04), ylim=c(0, 0.04), asp=1,
     xlab=expression(estimated~q[01]),
     ylab=expression(estimated~q[10]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[, "q01"], est_asym_b[, "q10"], pch=17, cex=0.7)

abline(v=0.01, lty=1, lwd=1)
abline(h=0.01, lty=1, lwd=1)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym["q01"],  m_sym["q10"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["q01"], m_asym["q10"], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("State change rates", side=3, adj=1, line=-1, cex=0.95)


##############q##############
pars_sym <- c(lambda0=0.10, lambda1=0.10, mu0=0.03, mu1=0.03, q01=0.01, q10=0.01)
pars_asym_q <- c(lambda0=0.10, lambda1=0.10, mu0=0.03, mu1=0.03, q01=0.01, q10=0.005)

set.seed(1)
n <- 100
tips <- 300

sim_sym <- vector("list", n)
i <- 1L; attempts <- 0L
while (i <= n) {
  attempts <- attempts + 1L
  phy <- tree.bisse(pars_sym, max.taxa=tips, max.t=Inf, include.extinct=FALSE, x0=NA)
  if (!is.null(phy) && length(phy$tip.label)==tips && sum(phy$tip.state)/tips <= 0.8 && sum(phy$tip.state)/tips >= 0.2) { sim_sym[[i]] <- phy; i <- i+1L }
  if (attempts %% 50L == 0L) message(sprintf("[sym] Kept %d/%d (attempts=%d)", i-1L, n, attempts))
}

sim_asym <- vector("list", n)
i <- 1L; attempts <- 0L
while (i <= n) {
  attempts <- attempts + 1L
  phy <- tree.bisse(pars_asym_q, max.taxa=tips, max.t=Inf, include.extinct=FALSE, x0=NA)
  if (!is.null(phy) && length(phy$tip.label)==tips && sum(phy$tip.state)/tips <= 0.8 && sum(phy$tip.state)/tips >= 0.2) { sim_asym[[i]] <- phy; i <- i+1L }
  if (attempts %% 50L == 0L) message(sprintf("[asym-q] Kept %d/%d (attempts=%d)", i-1L, n, attempts))
}

par_names <- c("lambda0","lambda1","mu0","mu1","q01","q10")
est_sym  <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))
est_asym <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))

message("Fitting symmetric set...")
for (k in seq_len(n)) {
  phy <- sim_sym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  
  lik <- make.bisse.mbt(phy, states)
  start <- starting.point.bisse(phy)
  
  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_sym[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}
save.image("/Volumes/LaCie/hui wang/result/mbt/transition/res_sym_notraitconst_nosurv")
message("Fitting asymmetric-q set...")
for (k in seq_len(n)) {
  phy <- sim_asym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  
  lik <- make.bisse.mbt(phy, states)
  start <- starting.point.bisse(phy)
  
  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_asym[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}
save.image("/Volumes/LaCie/hui wang/result/mbt/transition/res_asym_notraitconst_nosurv")
save(est_asym,file="/Volumes/LaCie/hui wang/result/mbt/transition/asym_notraitconst_nosurv")
#
est_sym  <- est_sym[complete.cases(est_sym), , drop=FALSE]
est_asym <- est_asym[complete.cases(est_asym), , drop=FALSE]

## -------- Speciation rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[,1], est_sym[,2],
     xlim=c(0,0.4), ylim=c(0,0.4), asp=1,
     xlab=expression(estimated~lambda[0]),
     ylab=expression(estimated~lambda[1]),
     pch=1, cex=0.7, col="black")
points(est_asym[,1], est_asym[,2], pch=17, cex=0.7)

abline(v=0.10, lty=1, lwd=1)
abline(h=0.10, lty=1, lwd=1)

m_sym  <- colMeans(est_sym);  m_asym <- colMeans(est_asym)
points(m_sym[1],  m_sym[2],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym[1], m_asym[2], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("Speciation rates", side=3, adj=1, line=-1, cex=0.95)

## -------- Extinction rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[,3], est_sym[,4],
     xlim=c(0,0.15), ylim=c(0,0.15), asp=1,
     xlab=expression(estimated~mu[0]),
     ylab=expression(estimated~mu[1]),
     pch=1, cex=0.7, col="black")
points(est_asym[,3], est_asym[,4], pch=17, cex=0.7)

abline(v=0.03, lty=1, lwd=1)  # mu0 (symmetric)
abline(h=0.03, lty=1, lwd=1)  # mu1 (symmetric)

m_sym  <- colMeans(est_sym);  m_asym <- colMeans(est_asym)
points(m_sym["mu0"],  m_sym["mu1"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["mu0"], m_asym["mu1"], pch=24, bg="white", cex=1.6, lwd=1.2)

mtext("Extinction rates", side=3, adj=1, line=-1, cex=0.95)

## -------- State change rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym[, "q01"], est_sym[, "q10"],
     xlim=c(0, 0.06), ylim=c(0, 0.06), asp=1,
     xlab=expression(estimated~q[01]),
     ylab=expression(estimated~q[10]),
     pch=1, cex=0.7, col="black")
points(est_asym[, "q01"], est_asym[, "q10"], pch=17, cex=0.7)

abline(v=0.01, lty=1, lwd=1)
abline(h=0.01, lty=1, lwd=1)
abline(h=0.005, lty=2, lwd=1)

m_sym  <- colMeans(est_sym);  m_asym <- colMeans(est_asym)
points(m_sym["q01"],  m_sym["q10"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["q01"], m_asym["q10"], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("State change rates", side=3, adj=1, line=-1, cex=0.95)

##BISSE
par_names <- c("lambda0","lambda1","mu0","mu1","q01","q10")
est_sym_b  <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))
est_asym_b <- matrix(NA_real_, n, 6, dimnames=list(NULL, par_names))

message("Fitting symmetric set...")
for (k in seq_len(n)) {
  phy <- sim_sym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  
  lik <- make.bisse(phy, states)
  start <- starting.point.bisse(phy)
  
  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_sym_b[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}

message("Fitting asymmetric-q set...")
for (k in seq_len(n)) {
  phy <- sim_asym[[k]]
  states <- phy$tip.state; if (is.factor(states)) states <- as.numeric(as.character(states))
  
  lik <- make.bisse(phy, states)
  start <- starting.point.bisse(phy)
  
  if (requireNamespace("subplex", quietly=TRUE)) {
    fit <- try(find.mle(lik, start, method="subplex",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  } else {
    fit <- try(find.mle(lik, start, method="Nelder-Mead",
                        control=list(maxit=1000, reltol=1e-6)), TRUE)
  }
  est_asym_b[k,] <- if (inherits(fit,"try-error")) rep(NA_real_, 6) else coef(fit)[par_names]
  if (k %% 10 == 0) message(sprintf("  %d/%d", k, n))
}

est_sym_b  <- est_sym_b[complete.cases(est_sym_b), , drop=FALSE]
est_asym_b <- est_asym_b[complete.cases(est_asym_b), , drop=FALSE]

## -------- Speciation rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[,1], est_sym_b[,2],
     xlim=c(0,0.4), ylim=c(0,0.4), asp=1,
     xlab=expression(estimated~lambda[0]),
     ylab=expression(estimated~lambda[1]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[,1], est_asym_b[,2], pch=17, cex=0.7)

abline(v=0.10, lty=1, lwd=1)
abline(h=0.10, lty=1, lwd=1)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym[1],  m_sym[2],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym[1], m_asym[2], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("Speciation rates", side=3, adj=1, line=-1, cex=0.95)

## --------Extinction rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[,3], est_sym_b[,4],
     xlim=c(0,0.15), ylim=c(0,0.15), asp=1,
     xlab=expression(estimated~mu[0]),
     ylab=expression(estimated~mu[1]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[,3], est_asym_b[,4], pch=17, cex=0.7)

abline(v=0.03, lty=1, lwd=1)  # mu0 (symmetric)
abline(h=0.03, lty=1, lwd=1)  # mu1 (symmetric)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym["mu0"],  m_sym["mu1"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["mu0"], m_asym["mu1"], pch=24, bg="white", cex=1.6, lwd=1.2)

mtext("Extinction rates", side=3, adj=1, line=-1, cex=0.95)

## -------- State change rates --------
par(mar=c(4.2,4.2,1.2,1))
plot(est_sym_b[, "q01"], est_sym_b[, "q10"],
     xlim=c(0, 0.06), ylim=c(0, 0.06), asp=1,
     xlab=expression(estimated~q[01]),
     ylab=expression(estimated~q[10]),
     pch=1, cex=0.7, col="black")
points(est_asym_b[, "q01"], est_asym_b[, "q10"], pch=17, cex=0.7)

abline(v=0.01, lty=1, lwd=1)
abline(h=0.01, lty=1, lwd=1)
abline(h=0.005, lty=2, lwd=1)

m_sym  <- colMeans(est_sym_b);  m_asym <- colMeans(est_asym_b)
points(m_sym["q01"],  m_sym["q10"],  pch=21, bg="white", cex=1.6, lwd=1.2)
points(m_asym["q01"], m_asym["q10"], pch=24, bg="white", cex=1.6, lwd=1.2)
mtext("State change rates", side=3, adj=1, line=-1, cex=0.95)
