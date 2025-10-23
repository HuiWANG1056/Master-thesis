## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## Common other functions include:
##   stationary.freq
##   starting.point
##   branches
library(diversitree)
## 1: make
make.hiquasse <- function(tree, states, states.sd, lambda, mu,
                        control=NULL, sampling.f=NULL) {
  cache <- diversitree:::make.cache.quasse(tree, states, states.sd, lambda, mu,
                             control, sampling.f)
  all_branches <- make.all_branches.hiquasse(cache, cache$control)
  rootfunc <- diversitree:::make.rootfunc.quasse(cache)
  f.pars <- diversitree:::make.pars.quasse(cache)
  
  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.f=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all_branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.f, intermediates)
  }
  
  class(ll) <- c("quasse", "dtlik", "function")
  ll
}

initial.tip.hiquasse <- function(cache, control, x) {
  nx <- control$nx * control$r
  npad <- nx - length(x)
  e0 <- 1 - cache$sampling.f
  dx <- control$dx
  
  if (!exists("%||%", mode = "function"))
    `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
  
  obs   <- control$obs %||% list()
  Tcut  <- obs$threshold %||% 0
  link  <- obs$link %||%"probit"
  k     <- obs$k %||% 1
  sd_default <- obs$sd %||% NA_real_
  prior_ctrl <- control$prior %||% list()
  prior_fun  <- prior_ctrl$fn   #p(x)
  if (!is.function(prior_fun)) {
    s0 <- prior_ctrl$sd   
    prior_fun <- function(xx, cache_) stats::dnorm(xx, mean = control$xmid, sd = s0)
  }
  prior_vec <- prior_fun(x, cache)
  
  p_y1_given_x <- function(xx, sigma_i) {    #observed: p(Y=1∣x)
    if (link == "probit") {
      stats::pnorm((xx - Tcut)/sigma_i)    #probit
    } else {
      stats::plogis(k * (xx - Tcut) / sigma_i)  #logit
    }
  }
  
  tip_from_Y <- function(Yi, sigma_i) {     #p(x|Y)
    sig <- if (is.na(sigma_i)) sd_default else sigma_i
    if (is.na(sig) || sig <= 0) stop("Observation-layer sigma must be > 0 ...")
    py1 <- p_y1_given_x(x, sig)
    w   <- if (Yi == 1) prior_vec * py1 else prior_vec * (1 - py1)
	w / (1-stats::pnorm(Tcut/sqrt(sigma_i^2+1),mean=0,sd=1)) 
  }
  
  if (control$tips.combined) {
    tips <- cache$tips
    t    <- cache$len[tips]
    i    <- order(t)
    target    <- tips[i]
    states    <- cache$states[i]        # now Y (0/1) at tips
    states.sd <- cache$states.sd[i]    
    y <- mapply(function(Yi, sdi)
      c(tip_from_Y(Yi, sdi), rep(0, npad)),
      states, states.sd, SIMPLIFY = FALSE)
    y <- matrix(c(rep(e0, nx), unlist(y)), nx, length(target) + 1)
    list(target = target, y = y, t = t[i])
  } else {
    y <- mapply(function(Yi, sdi)
      c(rep(e0, nx), tip_from_Y(Yi, sdi), rep(0, npad)),
      cache$states, cache$states.sd, SIMPLIFY = FALSE)
    diversitree:::dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
  }
}


## Extra core stuff:
make.all_branches.hiquasse <- function(cache, control) {
  branches <- diversitree:::make.branches.quasse(cache, control)
  initial.conditions <- diversitree:::make.initial.conditions.quasse(control)
  ## TODO: This is where tips.combined goes, *not* in the likelihood
  ## function...
  function(pars, intermediates, preset=NULL) {
    cache$y <- initial.tip.hiquasse(cache, cache$control, pars[[1]]$x)
    diversitree:::all_branches.list(pars, cache, initial.conditions,
                      branches, preset)
  }
}
