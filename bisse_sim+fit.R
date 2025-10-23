# sim_bisse.R =====================================
library(diversitree)
set.seed(123)

Ntips_set <- c(50, 100, 200)
R         <- 10
pars      <- c(lambda0=0.20, lambda1=0.30, mu0=0.05, mu1=0.05, q01=0.02, q10=0.02)

# ---- latent grid ----
Tcut      <- 0.0     
sigma_obs <- 0.25   
Nx   <- 1024          
rpad <- 4             
L    <- 8             
xmin <- Tcut - L * sigma_obs
xmax <- Tcut + L * sigma_obs
dx   <- (xmax - xmin) / (Nx - 1)
xmid <- 0.5 * (xmin + xmax)

# ---- output ----
outdir <- "out/sim/bisse"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (N in Ntips_set) {
  for (rep in 1:R) {
    cat("Simulating BiSSE: Ntips =", N, "rep =", rep, "\n")
    
    tree_ok <- FALSE
    for (attempt in 1:20) {
      seed <- sample.int(.Machine$integer.max, 1)
      set.seed(seed)
      sim <- try(tree.bisse(pars, max.taxa = N, x0 = NA, include.extinct = FALSE),
                 silent = TRUE)
      if (inherits(sim, "phylo")) { tree_ok <- TRUE; break }
    }
    
    if (!tree_ok) {
      warning("Simulation failed for N = ", N, " rep = ", rep)
      next
    }
    
    ntips     <- length(sim$tip.label)
    states01  <- as.integer(sim$tip.state)             
    states_sd <- setNames(rep(sigma_obs, ntips), sim$tip.label)
    
    out <- list(
      tree = sim,
      tip_states_01 = states01,
      states_sd = states_sd,
      meta = list(
        true_model = "BiSSE",
        Ntips = N, rep = rep, pars = pars, seed = seed,
        Tcut = Tcut, sigma_obs = sigma_obs,
        grid = list(Nx = Nx, r = rpad, xmin = xmin, xmax = xmax,
                    dx = dx, xmid = xmid)
      )
    )
    
    outfile <- sprintf("%s/bisse_N%03d_rep%03d.rds", outdir, N, rep)
    saveRDS(out, file = outfile)
  }
}


# fit_bisse.R =================================
files <- list.files("out/sim/bisse", pattern="\\.rds$", full.names=TRUE)
stopifnot(length(files) > 0)

for (f in files) {
  x   <- readRDS(f)
  phy <- x$tree
  
  # ---------- discrete states for BiSSE / CID ----------
  states01 <- setNames(as.integer(x$tip_states_01), phy$tip.label)
  n        <- length(states01)
  
  # ---------- BiSSE ----------
  lik.b <- make.bisse(phy, states01)
  p0.b  <- starting.point.bisse(phy)
  fit.b <- try(find.mle(lik.b, p0.b), silent = TRUE)
  
  ll.b  <- if (inherits(fit.b, "try-error")) NA_real_ else as.numeric(logLik(fit.b))
  k.b   <- if (inherits(fit.b, "try-error")) NA_integer_ else length(coef(fit.b))
  AIC.b <- if (is.na(ll.b)) NA_real_ else 2*k.b - 2*ll.b
  AICc.b<- if (is.na(AIC.b) || n <= k.b + 1) NA_real_ else AIC.b + (2*k.b*(k.b+1))/(n - k.b - 1)
  
  # ---------- CID-like (lambda0=lambda1, mu0=mu1) ----------
  lik.c <- constrain(lik.b, lambda1 ~ lambda0, mu1 ~ mu0)
  p0.c  <- if (inherits(fit.b, "try-error")) p0.b else coef(fit.b)
  p0.c  <- p0.c[c("lambda0","mu0","q01","q10")]
  fit.c <- try(find.mle(lik.c, p0.c), silent = TRUE)
  
  ll.c  <- if (inherits(fit.c, "try-error")) NA_real_ else as.numeric(logLik(fit.c))
  k.c   <- if (inherits(fit.c, "try-error")) NA_integer_ else length(coef(fit.c))
  AIC.c <- if (is.na(ll.c)) NA_real_ else 2*k.c - 2*ll.c
  AICc.c<- if (is.na(AIC.c) || n <= k.c + 1) NA_real_ else AIC.c + (2*k.c*(k.c+1))/(n - k.c - 1)
  
  # ---------- HiQuaSSE on the same dataset ----------
  lambda <- sigmoid.x
  mu     <- constant.x
  
  states01   <- setNames(as.integer(x$tip_states_01), phy$tip.label)
  states.num <- setNames(as.numeric(states01),          phy$tip.label)
  
  states.sd  <- setNames(as.numeric(x$states_sd[phy$tip.label]), phy$tip.label)
  stopifnot(identical(names(states.num), phy$tip.label),
            identical(names(states.sd),  phy$tip.label),
            all(is.finite(states.sd)), all(states.sd > 0))
  
  
  Tcut      <- if (!is.null(x$meta$Tcut)) x$meta$Tcut else 0.0
  sigma_obs <- if (!is.null(x$meta$sigma_obs)) x$meta$sigma_obs else 0.25
  
  if (!is.null(x$meta$grid)) {
    Nx   <- x$meta$grid$Nx
    rpad <- x$meta$grid$r
    dx   <- x$meta$grid$dx
    xmid <- x$meta$grid$xmid
  } else {
    Nx <- 1024; rpad <- 4
    L  <- 12
    xmin <- Tcut - L*sigma_obs
    xmax <- Tcut + L*sigma_obs
    dx   <- (xmax - xmin)/(Nx - 1)
    xmid <- 0.5*(xmin + xmax)
  }

  ctrl.h <- list(
    obs   = list(threshold = Tcut, link = "probit", sd = NA_real_),
    prior = list(sd = 1),
    nx = Nx, r = rpad, dx = dx, xmid = xmid,
    tips.combined = FALSE,
    method = "fftC"
  )
  
  lik.h0 <- make.hiquasse(phy, states.num, states.sd, lambda, mu, control = ctrl.h)
  lik.h  <- if ("drift" %in% argnames(lik.h0)) constrain(lik.h0, drift ~ 0) else lik.h0
  
  nm0   <- argnames(lik.h0)
  nm    <- argnames(lik.h)
  
  p.const <- suppressWarnings(starting.point.quasse(phy, states.num))
  safe    <- function(z, alt) if (is.finite(z)) z else alt
  mu_bar  <- if (!inherits(fit.b, "try-error")) {
    safe(mean(coef(fit.b)[c("mu0","mu1")]), p.const["mu"])
  } else p.const["mu"]
  
  seed_known <- c(
    "l.y0"      = if (!inherits(fit.b, "try-error"))
      safe(coef(fit.b)["lambda0"], p.const["lambda"]) else p.const["lambda"],
    "l.y1"      = if (!inherits(fit.b, "try-error"))
      safe(coef(fit.b)["lambda1"], p.const["lambda"]) else p.const["lambda"],
    "l.xmid"    = xmid,
    "l.r"       = 1,
    "m.c"       = mu_bar,
    "diffusion" = max(safe(p.const["diffusion"], 0.05), 1e-6),
    "drift"     = 0
  )
  seed_fallback <- function(name) {
    if (grepl("^l\\.", name))        return(0.1)
    if (grepl("^m\\.", name))        return(mu_bar)
    if (name == "diffusion")         return(0.05)
    if (name == "drift")             return(0)
    return(0.1)
  }
  p.start <- setNames(numeric(length(nm)), nm)
  for (k in nm) p.start[k] <- if (k %in% names(seed_known)) seed_known[k] else seed_fallback(k)
  p.start[p.start <= 0 & names(p.start) %in% c("l.y0","l.y1","m.c","diffusion")] <- 1e-8
  lower <- rep(-Inf, length(nm))
  names(lower) <- nm
  upper <- rep( Inf, length(nm))
  names(upper) <- nm
  for (nm_pos in intersect(c("l.y0","l.y1","m.c","diffusion"), nm)) lower[nm_pos] <- 1e-8
  if ("l.r"    %in% nm) { lower["l.r"] <- 1e-6; upper["l.r"] <- 50 }
  if ("l.xmid" %in% nm) { lower["l.xmid"] <- xmid - 5*sigma_obs; upper["l.xmid"] <- xmid + 5*sigma_obs }
  
  ll0 <- tryCatch(lik.h(p.start),
                  error = function(e) {cat("[HiQuaSSE] start error:", e$message, "\n"); NA_real_})
  if (!is.finite(ll0)) {
    cat("[HiQuaSSE] Non-finite logLik at start for", basename(f), "\n")
    next
  }
  
  lik.h.safe <- protect(lik.h, fail.value.default = -1e100)
  fit.h <- try(find.mle(lik.h.safe, p.start,
                        lower = lower, upper = upper,
                        method = "subplex",
                        control = list(parscale = pmax(abs(p.start), 1)),
                        finite = TRUE),
               silent = TRUE)
  
  ll.h  <- if (inherits(fit.h, "try-error")) NA_real_ else as.numeric(logLik(fit.h))
  k.h   <- if (inherits(fit.h, "try-error")) NA_integer_ else length(coef(fit.h))
  AIC.h <- if (is.na(ll.h)) NA_real_ else 2*k.h - 2*ll.h
  AICc.h<- if (is.na(AIC.h) || n <= k.h + 1) NA_real_ else AIC.h + (2*k.h*(k.h+1))/(n - k.h - 1)
  
  # ---------- summary ----------
  res <- data.frame(model = c("BiSSE","CID-like","HiQuaSSE"),
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
      "; \u0394AICc(second-best) =", if (is.na(second)) NA_real_ else round(res$AICc[second]-res$AICc[best], 3),
      "; Strong support (ΔAICc>2)?", strong, "\n")
}
