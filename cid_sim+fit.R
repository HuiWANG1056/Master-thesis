# sim_cidlike.R ==========================================
library(diversitree)

set.seed(123)
Ntips_set <- c(50, 100, 200)
R <- 10


lambda_A <- 0.18; lambda_B <- 0.28
mu_A     <- 0.04; mu_B     <- 0.04


q01 <- 0.02  # observed 0<->1 within same hidden class
qAB <- 0.01  # hidden A<->B within same observed state

# States: 1=0A, 2=1A, 3=0B, 4=1B
lambdas <- c(lambda_A, lambda_A, lambda_B, lambda_B)
mus     <- c(mu_A,     mu_A,     mu_B,     mu_B)
q <- c(
  q12 = q01, q13 = qAB, q14 = 0,
  q21 = q01, q23 = 0,   q24 = qAB,
  q31 = qAB, q32 = 0,   q34 = q01,
  q41 = 0,   q42 = qAB, q43 = q01
)
pars <- c(lambdas, mus, unname(q))  

outdir <- "out/sim/cidlike"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (N in Ntips_set) {
  for (rep in 1:R) {
    cat("Simulating CID-like (MuSSE 4-state): Ntips =", N, "rep =", rep, "\n")
    
    ok <- FALSE
    for (attempt in 1:10) {
      sim <- try(
        tree.musse(pars = pars, max.taxa = N, x0 = 1, include.extinct = FALSE),
        silent = TRUE
      )
      if (inherits(sim, "phylo")) { ok <- TRUE; break }
    }
    if (!ok) { warning("Simulation failed for N=", N, " rep=", rep); next }
    
    # Map 4 hidden states -> observed binary: {0A,0B}->0 ; {1A,1B}->1
    raw4  <- setNames(as.integer(sim$tip.state), sim$tip.label)
    obs01 <- setNames(ifelse(raw4 %in% c(1, 3), 0L, 1L), names(raw4))
    
    if (length(sim$tip.label) != N)
      warning("Got ", length(sim$tip.label), " tips; expected ", N)
    
    out <- list(
      tree           = sim,
      tip_states_01  = obs01,   
      tip_states_4   = raw4,   
      meta = list(
        true_model = "CIDlike",
        Ntips = N, rep = rep,
        lambda_A = lambda_A, lambda_B = lambda_B,
        mu_A = mu_A, mu_B = mu_B,
        q01 = q01, qAB = qAB,
        transitions = "no dual (1<->4, 2<->3 set 0)",
        diversitree_version = as.character(packageVersion("diversitree"))
      )
    )
    saveRDS(out, sprintf("%s/cidlike_N%d_rep%03d.rds", outdir, N, rep))
  }
}

# fit_cidlike.R =================================
files <- list.files("out/sim/cidlike", pattern="\\.rds$", full.names=TRUE)
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
  states.sd  <- if (!is.null(x$states_sd)) x$states_sd else setNames(rep(0.3, n), names(states))
  
  ctrl.h <- list(
    obs   = list(threshold = x$meta$Tcut, link = "probit", sd = if (!is.null(x$states_sd)) NA_real_ else 0.3),
    prior = list(sd = 1),
    nx = 1024, r = 4, method = "fftC"
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
