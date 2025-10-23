suppressPackageStartupMessages(library(deSolve))

if (!exists("tree") || !exists("states")) {
  set.seed(1)
  pars0  <- c(lambda0=0.10, lambda1=0.10, mu0=0.03, mu1=0.03, q01=0.01, q10=0.01)
  tree   <- tree.bisse(pars0, max.taxa=100, max.t=Inf, include.extinct=FALSE, x0=NA)
  states <- setNames(as.integer(tree$tip.state), tree$tip.label)
}

lik <- make.bisse.mbt(tree, states)
p0     <- starting.point.bisse(tree)
fit_ml <- find.mle(lik, p0)                        
hat    <- coef(fit_ml)                             
ll_hat <- as.numeric(logLik(fit_ml))     


mu0_max <- max(0.30, 2*hat["mu0"])
mu1_max <- max(0.30, 2*hat["mu1"])
mu0_grid <- seq(0, mu0_max, length.out = 51)
mu1_grid <- seq(0, mu1_max, length.out = 51)




free_names <- c("lambda0","lambda1","q01","q10")
Lprof <- matrix(NA_real_, length(mu0_grid), length(mu1_grid),
                dimnames=list(mu0=round(mu0_grid,4), mu1=round(mu1_grid,4)))

start_free <- hat[free_names]   
pb <- txtProgressBar(min=0, max=length(mu0_grid)*length(mu1_grid), style=3)
k  <- 0L
for (i in seq_along(mu0_grid)) {
  for (j in seq_along(mu1_grid)) {
    mu0_val <- mu0_grid[i]; mu1_val <- mu1_grid[j]
    lik_fix <- constrain(lik, mu0 ~ mu0_val, mu1 ~ mu1_val)  
    fit_ij  <- try(find.mle(lik_fix, start_free), silent=TRUE) 
    if (inherits(fit_ij, "try-error")) {                        
      fit_ij <- try(find.mle(lik_fix, hat[free_names]), silent=TRUE)
    }
    if (!inherits(fit_ij, "try-error")) {
      Lprof[i, j] <- as.numeric(logLik(fit_ij))                 
      start_free  <- coef(fit_ij)[free_names]                   
    }
    k <- k + 1L; setTxtProgressBar(pb, k)
  }
}
close(pb)


op <- par(no.readonly=TRUE)

persp(x=mu0_grid, y=mu1_grid, z=Lprof, theta=40, phi=25,
      xlab=expression(mu[0]), ylab=expression(mu[1]),
      zlab=expression(ell[prof](mu[0],mu[1])),
      ticktype="detailed", expand=0.75, shade=0.2,
      col="lightgrey", border="grey40")

par(op)

Delta <- 2*(ll_hat - Lprof)            
lev95 <- 5.991                      
image(mu0_grid, mu1_grid, Delta,
      xlab=expression(mu[0]), ylab=expression(mu[1]),
      col = hcl.colors(60, "YlOrRd", rev=TRUE))


contour(mu0_grid, mu1_grid, Delta, levels=lev95, add=TRUE, drawlabels=TRUE, lwd=2)
points(hat["mu0"], hat["mu1"], pch=8, cex=1.1); text(hat["mu0"], hat["mu1"], " MLE", pos=4, cex=0.9)

prof_out <- list(mu0_grid=mu0_grid, mu1_grid=mu1_grid,
                 Lprof=Lprof, Delta=Delta,
                 par_global=hat, lnLik_global=ll_hat)