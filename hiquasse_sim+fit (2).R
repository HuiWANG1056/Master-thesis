# sim_hiquasse.R =================================
library(diversitree)
set.seed(123)
Ntips_set <- c(50, 100, 200)
R <- 10
sig2 <- 0.1
lambda_logistic <- function(x) 0.1 + 0.2 / (1 + exp(-5 * (x - 0)))  # xmid=0
mu_const        <- function(x) rep(0.05, length(x))

Tcut      <- 1        
sigma_obs <- 0.25     
link      <- "probit" 

outdir <- "out/sim/hiquasse"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (N in Ntips_set) {
  accepted <- 0                                
  while (accepted < R) {                       
    cat("Simulating HiQuaSSE: Ntips =", N, "target rep =", accepted + 1, "\n")  

    tree_ok <- FALSE
    for (attempt in 1:10) {
      seed <- sample.int(.Machine$integer.max, 1)
      set.seed(seed)
      
      pars <- list(
        lambda_logistic,
        mu_const,
        make.brownian.with.drift(drift = 0, diffusion = sig2)
      )
      
      sim <- try(tree.quasse(pars = pars, max.taxa = N,
                             x0 = 0, include.extinct = FALSE,
                             single.lineage = FALSE, verbose = FALSE),
                 silent = TRUE)
      if (inherits(sim, "phylo")) { tree_ok <- TRUE; break }
    }
    if (!tree_ok) { 
      warning("Simulation failed for N=", N, " (no valid phylo after 10 attempts)"); 
      next                                
    }
    
    x_vals <- setNames(as.numeric(sim$tip.state), sim$tip.label)
    p1 <- if (link == "probit") pnorm((x_vals - Tcut)/sigma_obs) else plogis((x_vals - Tcut)/sigma_obs)
    y_vals <- setNames(rbinom(length(p1), 1, p1), names(p1))
    
    prop0 <- mean(y_vals == 0)                 
    if (prop0 < 0.20 || prop0 > 0.80) {         
      cat("  Discarded (prop0 =", round(prop0, 3), "outside [0.20, 0.80])\n")  
      next                                      
    }
    

    accepted <- accepted + 1                    
    rep <- accepted                              
    
    states_sd <- setNames(rep(sigma_obs, length(p1)), names(p1))  
    out <- list(
      tree = sim,
      tip_trait_x = x_vals,
      tip_states_01 = as.integer(y_vals),
      states_sd = states_sd,
      meta = list(
        true_model="hiquasse", Ntips=N, rep=rep, sig2=sig2,
        lambda_type="logistic(xmid=0)", mu=0.05, Tcut=Tcut,
        sigma_obs=sigma_obs, link=link, seed=seed,
        prop0=prop0,                              
        diversitree_version=as.character(packageVersion("diversitree"))
      )
    )
    
    saveRDS(out, sprintf("%s/hiquasse_N%d_rep%03d.rds", outdir, N, rep))
    cat("  Saved: rep =", rep, "prop0 =", round(prop0, 3), "\n")  
  }
}



# fit_hiquasse.R =================================
files <- list.files("out/sim/hiquasse", pattern="\\.rds$", full.names=TRUE)
stopifnot(length(files) > 0)

for (f in files) {
  x <- readRDS(f)
  phy <- x$tree
  states <- setNames(as.integer(x$tip_states_01), phy$tip.label)
  n <- length(states)
  
  ## ------- BiSSE -------
  lik.b <- make.bisse(phy, states)
  p0.b  <- starting.point.bisse(phy)
  fit.b <- try(find.mle(lik.b, p0.b), silent = TRUE)
  ll.b  <- if (inherits(fit.b, "try-error")) NA_real_ else as.numeric(logLik(fit.b))
  k.b   <- if (inherits(fit.b, "try-error")) NA_integer_ else length(coef(fit.b))
  AIC.b <- if (is.na(ll.b)) NA_real_ else 2*k.b - 2*ll.b
  AICc.b<- if (is.na(AIC.b) || n <= k.b + 1) NA_real_ else AIC.b + (2*k.b*(k.b+1))/(n - k.b - 1)
  
  ## ------- CID-like -------
  lik.c <- constrain(lik.b, lambda1 ~ lambda0, mu1 ~ mu0)
  p0.c  <- if (inherits(fit.b, "try-error")) p0.b else coef(fit.b)
  p0.c  <- p0.c[c("lambda0","mu0","q01","q10")]
  fit.c <- try(find.mle(lik.c, p0.c), silent = TRUE)
  ll.c  <- if (inherits(fit.c, "try-error")) NA_real_ else as.numeric(logLik(fit.c))
  k.c   <- if (inherits(fit.c, "try-error")) NA_integer_ else length(coef(fit.c))
  AIC.c <- if (is.na(ll.c)) NA_real_ else 2*k.c - 2*ll.c
  AICc.c<- if (is.na(AIC.c) || n <= k.c + 1) NA_real_ else AIC.c + (2*k.c*(k.c+1))/(n - k.c - 1)
  
  ## ------- hiquasse-------
  lambda <- sigmoid.x
  mu     <- constant.x
  
  states.num <- as.numeric(states)
  names(states.num) <- names(states)
  states.sd  <- if (!is.null(x$states_sd)) x$states_sd else setNames(rep(0.25, n), names(states))
  
  ctrl.h <- list(
    obs   = list(threshold = x$meta$Tcut, link = "probit", sd = states.sd),
    prior = list(sd = 1),
    nx = 1024, r = 4, method = "fftC",
    xmid = 0, dx=0.01
  )
  lik.h0 <- make.hiquasse(phy, states.num, states.sd, lambda, mu, control = ctrl.h)
  lik.h  <- constrain(lik.h0, drift ~ 0)
  
  p.const <- starting.point.quasse(phy, states.num)
  p.start <- c(p.const["lambda"], p.const["lambda"], mean(states.num), 1, p.const["mu"], p.const["diffusion"])
  names(p.start) <- argnames(lik.h)
  
  lower <- c(0, 0, min(states.num), -Inf, 0, 0)
  names(lower) <- argnames(lik.h)
  
  fit.h <- try(find.mle(lik.h, p.start, lower = lower, control = list(parscale = .1)), silent = TRUE)
  ll.h  <- if (inherits(fit.h, "try-error")) NA_real_ else as.numeric(logLik(fit.h))
  k.h   <- if (inherits(fit.h, "try-error")) NA_integer_ else length(coef(fit.h))
  AIC.h <- if (is.na(ll.h)) NA_real_ else 2*k.h - 2*ll.h
  AICc.h<- if (is.na(AIC.h) || n <= k.h + 1) NA_real_ else AIC.h + (2*k.h*(k.h+1))/(n - k.h - 1)
  
  res <- data.frame(model = c("BiSSE","CID-like","hiquasse"),
                    k = c(k.b, k.c, k.h),
                    logLik = c(ll.b, ll.c, ll.h),
                    AIC = c(AIC.b, AIC.c, AIC.h),
                    AICc = c(AICc.b, AICc.c, AICc.h),
                    check.names = FALSE)
  
  d <- res$AICc - min(res$AICc, na.rm = TRUE)
  w <- exp(-0.5 * d)
  w <- w / sum(w[is.finite(w)])
  res$dAICc <- round(d, 3)
  res$weight <- round(w, 3)
  
  ord <- order(res$AICc, na.last = TRUE)
  best <- ord[1]
  second <- if (length(ord) >= 2) ord[2] else NA_integer_
  strong <- if (!is.na(second) && is.finite(res$AICc[best]) && is.finite(res$AICc[second]))
    (res$AICc[second] - res$AICc[best] > 2) else NA
  
  cat("\n=== File:", basename(f), "| true:", x$meta$true_model, "| Ntips:", n, "===\n")
  print(res[order(res$AICc), ], row.names = FALSE)
  cat("Best by AICc:", res$model[best],
      "; ΔAICc(second-best) =", if (is.na(second)) NA_real_ else round(res$AICc[second]-res$AICc[best], 3),
      "; Strong support (ΔAICc>2)?", strong, "\n")
}
