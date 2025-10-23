fp <- file.path("~", "Downloads", "hiquassefit_updated")  

env <- new.env()
obj_names <- load(fp, envir = env)  
cat("Loaded objects:", paste(obj_names, collapse = ", "), "\n")

for (nm in obj_names) {
  cat("\n=== Object:", nm, "===\n")
  x <- env[[nm]]
  cat("class:", paste(class(x), collapse = " / "), "\n")
  if (is.data.frame(x) || is.matrix(x)) {
    print(dim(x)); print(utils::head(x, 5))
  } else if (is.list(x)) {
    utils::str(x, max.level = 1)
  } else {
    utils::str(x)
  }
}


#================== plot ==========================
library(ggplot2)

fits   <- env$fit
files  <- env$files
N_tag  <- as.integer(sub(".*_N([0-9]+)_rep.*", "\\1", files))  
xgrid  <- seq(-3, 3, length.out = 400)                         
lam_tr <- 0.1 + 0.2 / (1 + exp(-5 * (xgrid - 0)))              

for (N in c(50, 100, 200)) {
  idx <- which(N_tag == N)
  lam_mat <- sapply(idx, function(k) {
    pk <- tryCatch(coef(fits[[k]]), error = function(e) fits[[k]]$par)
    env$lambda(xgrid, pk["l.y0"], pk["l.y1"], pk["l.xmid"], pk["l.r"])
  })  
  
  mu  <- rowMeans(lam_mat, na.rm = TRUE)
  qs  <- apply(lam_mat, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  df_lines <- data.frame(
    x = rep(xgrid, times = ncol(lam_mat)),
    y = as.vector(lam_mat),
    rep = rep(seq_len(ncol(lam_mat)), each = length(xgrid))
  )
  df_mean <- data.frame(x = xgrid, y = mu)
  df_band <- data.frame(x = xgrid, lo = qs[1, ], hi = qs[2, ])
  df_true <- data.frame(x = xgrid, y = lam_tr)
  
  p <- ggplot() +
    geom_line(data = df_lines, aes(x, y, group = rep), alpha = 0.25, linewidth = 0.3) +
    geom_ribbon(data = df_band, aes(x, ymin = lo, ymax = hi), alpha = 0.15) +
    geom_line(data = df_mean, aes(x, y), linewidth = 1) +
    geom_line(data = df_true, aes(x, y), linewidth = 1, linetype = "dashed") +
    labs(title = paste0("HiQuaSSE: N = ", N), x = "Trait x", y = expression(lambda(x))) +
    theme_minimal()
  print(p)  
}


