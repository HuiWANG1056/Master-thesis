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

## 1: make
make.bisse.mbt <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10, strict=TRUE, control=list()) {
  cache <- make.cache.bisse.mbt(tree, states, unresolved, sampling.f,
                            nt.extra, strict)
  unresolved <- cache$unresolved

  rootfunc <- diversitree:::rootfunc.musse
  
  ll <- function(pars, condition.surv=FALSE, root=ROOT.FLAT,
                 root.p=NULL, intermediates=FALSE) {
    diversitree:::check.pars.bisse(pars)
    ans <- all.branches.matrix.mbt(pars, cache, control,
                                          initial.conditions.bisse.mbt)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse", "dtlik", "function")
  ll
}

## 3: make.cache (& initial.tip)
make.cache.bisse.mbt <- function(tree, states, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10,
                             strict=TRUE) {
  tree <- diversitree:::check.tree(tree)
  states <- diversitree:::check.states(tree, states,
                         strict=strict && is.null(unresolved),
                         strict.vals=0:1)
  
  if ( inherits(tree, "clade.tree") ) {
    if ( !is.null(unresolved) )
      stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
    unresolved <- diversitree:::make.unresolved.bisse(tree$clades, states)
    states <- states[tree$tip.label]
  }
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  
  cache <- diversitree:::make.cache(tree)
  cache$info <- diversitree:::make.info.bisse(tree)
  cache$states     <- states
  cache$sampling.f <- diversitree:::check.sampling.f(sampling.f, 2)
  
  cache$unresolved <- diversitree:::check.unresolved(cache, unresolved, nt.extra)
  if ( !is.null(cache$unresolved) ) {
    cache$tips   <- cache$tips[-cache$unresolved$i]
    cache$states <- cache$states[-cache$unresolved$i]
  }
  
  if ( strict && !is.null(unresolved) ) {
    seen <- ((0:1) %in% unique(na.omit(cache$states)) |
               sapply(cache$unresolved[c("n0", "n1")], sum) > 0)
    if ( !all(seen) )
      stop("Because strict state checking requested, all (and only) ",
           "states in 0 1, are allowed")
  }
  
  if ( length(cache$tips) > 0 ) # in case all tips are unresolved.
    cache$y <- initial.tip.mbt(cache)
  
  cache
}

## 4: initial.conditions:
## Note that we ignore both 't' and 'idx'.
initial.tip.mbt <- function(cache) {
  tips <- cache$tips
  states <- cache$states
  t <- rep(0, length(tips))
  y_list <- vector("list", length(tips))
  for (i in seq_along(tips)) {
    idx <- tips[i]
    st <- states[i]  # 0 or 1
    if (st == 0) {
      v1 <- c(0,0,1,0)
      v2 <- c(0,0,0,0)
    } else {
      v1 <- c(0,0,0,0)
      v2 <- c(0,0,0,1)
    }
    y_mat <- cbind(v1, v2)
    y_list[[i]] <- list(y=y_mat, t=0, target=idx)
  }
  y_list
}

initial.conditions.bisse.mbt <- function(init, pars, t, idx){
  # init: 4×2 matrix，rows=E0,E1,D0,D1; cols=2 sets of ODE 
  E0 <- init[1,1]
  E1 <- init[2,1]
  D0 <- 2 * pars[1] * init[3,1] * init[3,2]
  D1 <- 2 * pars[2] * init[4,1] * init[4,2]
  q <- D0+D1
  D0 <- D0/q
  D1 <- D1/q
  c( E0, E1, D0, D1, q)
}

                                          
all.branches.matrix.mbt <- function(pars, cache, control, initial.conditions,
                                preset=NULL) {
  control <- diversitree:::check.control.ode(control)                     	
  info <- cache$info
  
  n.var <- info$ny
  n.par <- info$np 
  rtol <- atol <- control$tol
  stepper <- control$gsl.stepper
  
  derivs <- info$derivs
  ode   <- new(diversitree:::GslOdeR,        derivs, environment(derivs), n.var)
  
  ode$set_control(list(atol=atol, rtol=rtol, algorithm=stepper,
                       hini=1e-4))
                          	
  len <- cache$len
  depth <- cache$depth
  children <- cache$children
  order <- cache$order[-length(cache$order)]
  root <- cache$root
  
  n <- length(len)
  lq <- rep(0,n)
  n.tip <- cache$n.tip
  
  y <- cache$y
  branch.init <- branch.base <- matrix(NA, cache$info$ny, n)
  G_list <- vector("list", n)
  
  ## TODO: It might be an idea here to check preset is OK:
  ## Must have names target, lq, base
  ## must be of same length.
  if ( !is.null(preset) ) {
    lq[preset$target] <- preset$lq
    branch.base[,preset$target] <- preset$base
  }
  
  if (is.null(names(y))) {
    for (x in y) {
      idx <- x$target
      y_mat <- x$y          

      out_mat <- matrix(NA, nrow=4, ncol=2)
      for (j in 1:2) {
        init_vec <- y_mat[,j]
        if (sum(init_vec[3:4]) > 0) {
          out_mat[,j] <- ode$run(c(0,len[idx]), init_vec, pars)
        } else {
          out_mat[,j] <- init_vec
        }
      }
      G_list[[idx]] <- out_mat
      branch.base[,idx] <- out_mat[,which.max(colSums(out_mat[3:4,]))]

    }
  } else {
    tip.t <- y$t
    tip.target <- y$target
    y_list <- y$y
    for (k in seq_along(tip.t)) {
      idx <- tip.target[k]
      y_mat <- y_list[[k]]          
      out_mat <- matrix(NA, nrow=4, ncol=2)
      for (j in 1:2) {
        init_vec <- y_mat[,j]
        if (sum(init_vec[3:4]) > 0) {
          out_mat[,j] <- ode$run(c(0,len[idx]), init_vec, pars)
        } else {
          out_mat[,j] <- init_vec
        }
      }
      G_list[[idx]] <- out_mat
	  branch.base[,idx] <- out_mat[,which.max(colSums(out_mat[3:4,]))]
    }
  }
  
  # Internal nodes
  for (i in order) {
    y.in <- initial.conditions(branch.base[,children[i,]], pars,
                               depth[i], i)
    branch.base[,i] <- y.in[1:4]
    lq[i] <- log(y.in[5])
    out1 <- ode$run(c(0,len[i]), c(y.in[1:2],1,0), pars)
    out2 <- ode$run(c(0,len[i]), c(y.in[1:2],0,1), pars)
    G_list[[i]] <- cbind(out1, out2)  
    branch.base[1:2,i] <- out1[1:2]
    branch.base[3:4,i] <- G_list[[i]][3:4,] %*% y.in[c(3,4)]
  }
  
  y.in <- initial.conditions(branch.base[,children[root,]], pars,
                             depth[root], root)
  branch.base[,root] <- y.in[1:4]
  lq[root] <- log(y.in[5])
  ret1 <- ode$run(c(0,depth[root]), c(y.in[1:2],1,0), pars)
  ret2 <- ode$run(c(0,depth[root]), c(y.in[1:2],0,1), pars)
  G_list[[root]] <- cbind(ret1, ret2) 
  branch.base[1:2,i] <- ret1[1:2]  
  branch.base[3:4,root] <- G_list[[root]][3:4,] %*% y.in[c(3,4)]
  
  list(vals = branch.base[,root], lq = sum(lq))
}

