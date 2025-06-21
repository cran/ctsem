# # Importance sampling for optimized ctstanfit
# importanceSamplingFunc <- function(
#     samples, mcovl, delta,
#   finishsamples, finishmultiply, isloopsize, tdf,
#   cl, cores, parlp, flexlapplytext, chancethreshold
# ) {
#   message('Importance sampling...')
#   log_sum_exp <- function(x) {
#     xmax <- which.max(x)
#     log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
#   }
#   targetsamples <- finishsamples * finishmultiply
#   j <- 0L
#   sample_prob <- numeric(0)
#   ess <- numeric(0)
#   qdiag <- numeric(0)
#   target_dens_list <- list()
# 
#   while (nrow(samples) < targetsamples) {
#     j <- j + 1L
#     if (j == 1L) {
#       samples    <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]], df = tdf)
#       newsamples <- samples
#     } else {
#       delta[[j]] <- colMeans(resamples)
#       mcovl[[j]] <- as.matrix(Matrix::nearPD(cov(resamples))$mat)
#       for (iteri in seq_len(j)) {
#         if (iteri == 1L) {
#           mcov <- mcovl[[j]]; mu <- delta[[j]]
#         } else {
#           mcov <- mcov * .5 + mcovl[[j]] * .5
#           mu   <- mu   * .5 + delta[[j]]   * .5
#         }
#       }
#       newsamples <- mvtnorm::rmvt(isloopsize, delta = mu, sigma = mcov, df = tdf)
#       samples    <- rbind(samples, newsamples)
#     }
# 
#     prop_dens <- mvtnorm::dmvt(tail(samples, isloopsize), delta[[j]], mcovl[[j]], df = tdf, log = TRUE)
#     if (cores > 1) parallel::clusterExport(cl, 'samples', envir = environment())
# 
#     target_dens_list[[j]] <- unlist(
#       flexlapplytext(
#         cl    = cl,
#         X     = seq_len(isloopsize),
#         fn    = "function(x){parlp(samples[x,])}",
#         cores = cores
#       )
#     )
#     target_dens_list[[j]][is.na(target_dens_list[[j]])] <- -1e200
#     if (all(target_dens_list[[j]] < -1e100)) stop('Could not sample from optimum! Try reparamaterizing?')
# 
#     targ <- target_dens_list[[j]]
#     targ[!is.finite(targ)] <- -1e30
#     weighted <- targ - prop_dens
# 
#     newsampleprob <- exp(weighted - log_sum_exp(weighted))
#     sample_prob   <- c(sample_prob, newsampleprob)
#     sample_prob[!is.finite(sample_prob)] <- 0
#     sample_prob[is.na(sample_prob)]     <- 0
# 
#     if (nrow(samples) >= targetsamples && max(sample_prob)/sum(sample_prob) < chancethreshold) {
#       message('Finishing importance sampling...')
#       nresamples <- finishsamples
#     } else {
#       nresamples <- max(5000, nrow(samples)/5)
#     }
# 
#     resample_i <- sample(
#       seq_len(nrow(samples)),
#       size    = nresamples,
#       replace = nrow(samples) <= targetsamples,
#       prob    = sample_prob / sum(sample_prob)
#     )
#     message(sprintf(
#       'Importance sample loop %d: %d unique samples from %d resamples of %d; prob sd=%.4f, max chance=%.4f',
#       j, length(unique(resample_i)), nresamples, nrow(samples),
#       sd(sample_prob), max(sample_prob)*isloopsize
#     ))
#     if (length(unique(resample_i)) < 100) {
#       message('Sampling ineffective (<100 unique) - consider increasing isloopsize or using HMC.')
#     }
# 
#     resamples <- samples[resample_i, , drop = FALSE]
#     ess[j]    <- sum(sample_prob[resample_i])^2 / sum(sample_prob[resample_i]^2)
#     qdiag[j]  <- mean(vapply(
#       seq_len(500),
#       function(i) {
#         idx <- sample(sample_prob, size = i, replace = TRUE)
#         max(idx) / sum(idx)
#       },
#       numeric(1)
#     ))
#   }
# 
#   combined_target <- unlist(target_dens_list)
#   lpsamples <- combined_target[resample_i]
# 
#   list(
#     samples     = samples,
#     resamples   = resamples,
#     resample_i  = resample_i,
#     lpsamples   = lpsamples,
#     sample_prob = sample_prob,
#     ess         = ess,
#     qdiag       = qdiag
#   )
# }


# # ---------------------------------------------------------------------------
# # AMIS version – fixed denominator – interrupt-safe
# # ---------------------------------------------------------------------------
# adaptive_is <- function(parlp,
#                         mu_hat,
#                         Sigma_hat,
#                         cl,
#                         n_batch     = 1000,
#                         target_ess  = 50,
#                         max_iter    = 10,
#                         scale_init  = 3,
#                         cov_inflate = 1.3,
#                         ridge       = 1e-8,
#                         finishsamples = 1000,
#                         use_t_df    = NULL,
#                         verbose     = TRUE) {
#   stopifnot(is.function(parlp),
#             is.numeric(mu_hat),
#             is.matrix(Sigma_hat))
#   cores=length(cl)
# 
#   if (!requireNamespace("mvtnorm", quietly = TRUE))
#     stop("Install the ‘mvtnorm’ package.")
#   if (!requireNamespace("diagis", quietly = TRUE))
#     stop("Install the ‘diagis’ package.")
# 
#   ## helpers -----------------------------------------------------------------
#   rmv <- if (is.null(use_t_df)) mvtnorm::rmvnorm else
#     function(n, mean, sigma) mvtnorm::rmvt(n, sigma, df = use_t_df, delta = mean)
#   dmv_log <- if (is.null(use_t_df))
#     function(x, m, S) mvtnorm::dmvnorm(x, m, S, log = TRUE) else
#       function(x, m, S) mvtnorm::dmvt(x, S, df = use_t_df, delta = m, log = TRUE)
#   ess  <- function(w) diagis::ess(w)
# 
#   logplus <- function(a, b) {                   # log(exp(a)+exp(b)), vectorised
#     idx <- a > b
#     res <- numeric(length(a))
#     res[idx]  <- a[idx] + log1p(exp(b[idx] - a[idx]))
#     res[!idx] <- b[!idx] + log1p(exp(a[!idx] - b[!idx]))
#     res
#   }
# 
#   ## storage -----------------------------------------------------------------
#   samples      <- matrix(0, 0, length(mu_hat))
#   log_p_target <- numeric(0)
#   log_qsum     <- numeric(0)
#   prop_mu_lst  <- list();  prop_cov_lst <- list()
#   T_props      <- 0L
#   ess_now      <- 0
#   interrupted  <- FALSE                      # <- NEW flag
# 
#   ## main loop ---------------------------------------------------------------
#   tryCatch({
#     for (iter in seq_len(max_iter)) {
# 
#       ## build proposal ------------------------------------------------------
#       if (iter == 1) {
#         prop_mu  <- mu_hat
#         prop_cov <- scale_init^2 * Sigma_hat
#       } else {
#         w_norm   <- w / sum(w)
#         prop_mu  <- as.numeric(diagis::weighted_mean(samples, w_norm))
#         prop_cov <- diagis::weighted_var(samples, w_norm) * cov_inflate
#       }
#       prop_cov <- prop_cov + diag(ridge, nrow(prop_cov))
# 
#       T_props <- T_props + 1L
#       prop_mu_lst[[T_props]]  <- prop_mu
#       prop_cov_lst[[T_props]] <- prop_cov
# 
#       ## update mixture for existing particles ------------------------------
#       if (nrow(samples) > 0L) {
#         lq_old  <- dmv_log(samples, prop_mu, prop_cov)
#         log_qsum <- logplus(log_qsum, lq_old)
#       }
# 
#       ## draw & score new batch ---------------------------------------------
#       x_new <- rmv(n_batch, prop_mu, prop_cov)
# 
#       if (cores > 1)
#         parallel::clusterExport(cl, 'x_new', envir = environment())
# 
#       lp_new <- unlist(
#         flexlapplytext(
#           cl    = cl,
#           X     = 1:nrow(x_new),
#           fn    = "function(i){parlp(x_new[i,])}",
#           cores = cores
#         ), use.names = FALSE
#       )
#       lp_new[is.na(lp_new)] <- -1e200
#       if (all(lp_new < -1e100))
#         stop("Could not sample from optimum! Try re-parameterising?")
# 
#       ## mixture denominator for new particles ------------------------------
#       log_mix_new <- rep(-Inf, n_batch)
#       if (T_props > 1L) {
#         for (k in seq_len(T_props - 1L))
#           log_mix_new <- logplus(
#             log_mix_new,
#             dmv_log(x_new, prop_mu_lst[[k]], prop_cov_lst[[k]])
#           )
#       }
#       lq_new      <- dmv_log(x_new, prop_mu, prop_cov)
#       log_mix_new <- logplus(log_mix_new, lq_new)
# 
#       ## append --------------------------------------------------------------
#       samples      <- rbind(samples, x_new)
#       log_p_target <- c(log_p_target, lp_new)
#       log_qsum     <- c(log_qsum,    log_mix_new)
# 
#       ## weights & ESS -------------------------------------------------------
#       log_mix <- log_qsum - log(T_props)
#       log_w   <- log_p_target - log_mix
#       log_w   <- log_w - max(log_w)
#       w       <- exp(log_w)
#       ess_now <- ess(w)
# 
#       if (verbose)
#         message(sprintf("Importance sample iter %2d | n = %6d | ESS = %.1f",
#                         iter, nrow(samples), ess_now))
# 
#       if (ess_now >= target_ess) break
#     }
#   }, interrupt = function(e) {               # <- NEW interrupt handler
#        interrupted <<- TRUE
#        if (verbose) message("⏹  interrupted – returning current state.")
#      })
# 
#   if (length(log_p_target) == 0)
#     stop("No samples collected before interrupt.")
# 
#   w_norm <- w / sum(w)
# 
#   ## resample equal-weight ---------------------------------------------------
#   idx_post <- sample.int(length(w_norm),
#                          size    = finishsamples,
#                          replace = TRUE,
#                          prob    = w_norm)
# 
#   list(theta       = samples[idx_post, , drop = FALSE],
#        ess         = ess_now,
#        interrupted = interrupted,           # <- flag for caller
#        mean        = as.numeric(diagis::weighted_mean(samples, w_norm)),
#        covariance  = diagis::weighted_var(samples, w_norm))
# }
# ======================================================================
#  adaptive_is()     –  AMIS with
#      * robust Student-t option (use_t_df = "auto")
#      * 2-component proposal (adaptive 70 %  +  base-tail 30 %)
#      * correct same-batch denominator update
#      * numerically stable mixture density
#      * correctly centred covariance update
#      * verbose ↑ / ↓ trace of proposal scale
# ======================================================================
# adaptive_is <- function(parlp,
#   mu_hat,
#   Sigma_hat,
#   cl,
#   n_batch       = 1000,
#   target_ess    = 100,
#   max_iter      = 50,
#   scale_init    = 3,
#   cov_inflate   = 1.3,
#   ridge         = 1e-8,
#   finishsamples = 1000,
#   use_t_df      = "auto",
#   power_upd     = 0.8,
#   window_size   = 3,     # batches kept in sliding window
#   ess_gate      = 3,     # ESS reference in multiples of d
#   lambda_ema    = 0.25,  # max EMA learning-rate
#   verbose       = TRUE) {
# 
#   stopifnot(is.function(parlp),
#             is.numeric(mu_hat),
#             is.matrix(Sigma_hat))
# 
#   cores <- length(cl)
#   d     <- length(mu_hat)
#   message(sprintf("Target ESS %d | max iter %d | Esc to interrupt.",
#                   target_ess, max_iter))
# 
#   for (pkg in c("mvtnorm", "diagis", "matrixStats"))
#     if (!requireNamespace(pkg, quietly = TRUE))
#       stop(sprintf("Install '%s' first.", pkg))
# 
#   ## ── helpers ------------------------------------------------------------
#   safe_pd <- function(S, eps = ridge) {
#     for (f in c(1,1e-4,1e-3,1e-2))
#       if (all(eigen(S + diag(eps*f, nrow(S)), symmetric=TRUE,
#                     only.values=TRUE)$values > 0))
#         return(S + diag(eps*f, nrow(S)))
#     diag(diag(S)) + diag(1e-2, nrow(S))
#   }
# 
#   make_k <- function(df) {
#     if (is.finite(df)) {
#       list(
#         rmv = function(n,m,S)
#           mvtnorm::rmvt(n, sigma = safe_pd(S), df = df, delta = m),
#         dmv = function(x,m,S,log=TRUE){
#           if (is.vector(x)) x <- matrix(x, ncol = length(m))
#           mvtnorm::dmvt(x, delta = m, sigma = safe_pd(S), df = df, log = log)
#         })
#     } else {
#       list(
#         rmv = function(n,m,S)
#           mvtnorm::rmvnorm(n, mean = m, sigma = safe_pd(S)),
#         dmv = function(x,m,S,log=TRUE){
#           if (is.vector(x)) x <- matrix(x, ncol = length(m))
#           mvtnorm::dmvnorm(x, mean = m, sigma = safe_pd(S), log = log)
#         })
#     }
#   }
# 
#   logplus <- function(a,b){
#     idx <- a > b
#     r   <- numeric(length(a))
#     r[idx]  <- a[idx] + log1p(exp(b[idx]-a[idx]))
#     r[!idx] <- b[!idx] + log1p(exp(a[!idx]-b[!idx]))
#     r
#   }
# 
#   stratified_rs <- function(w,N){
#     cs <- cumsum(w / sum(w))
#     findInterval((runif(1)+0:(N-1))/N, cs) + 1L
#   }
# 
#   ess <- function(w) diagis::ess(w)
# 
#   ## ── base component -----------------------------------------------------
#   base_mu  <- mu_hat
#   base_cov <- scale_init^2 * Sigma_hat
#   mix_w    <- 0.7
#   base_dmv <- make_k(Inf)$dmv
# 
#   ## ── state --------------------------------------------------------------
#   samples <- matrix(0,0,d)
#   log_p_target <- log_qsum <- numeric(0)
#   w <- numeric(0); ess_now <- 0; interrupted <- FALSE
# 
#   T_props <- 0L
#   current_df <- if (is.numeric(use_t_df)) use_t_df else Inf
#   kf <- make_k(current_df); rmv_fn <- kf$rmv; dmv_fn <- kf$dmv
#   prev_cov   <- base_cov
#   prev_trace <- sum(diag(prev_cov))
# 
#   win_x  <- vector("list", window_size)
#   win_w  <- vector("list", window_size)
#   win_ptr<- 0L
# 
#   ## ── loop ---------------------------------------------------------------
#   tryCatch({
#     for (iter in seq_len(max_iter)) {
# 
#       ## mean -------------------------------------------------------------
#       if (iter == 1 || length(w)==0) {
#         prop_mu <- mu_hat
#       } else {
#         w_up   <- (w^power_upd)/sum(w^power_upd)
#         prop_mu <- as.numeric(diagis::weighted_mean(samples, w_up))
#       }
# 
#       ## sliding-window covariance ---------------------------------------
#       if (iter == 1) {
#         new_cov <- base_cov
#       } else {
#         Xw <- do.call(rbind, win_x)
#         Ww <- unlist(win_w)
#         centred <- sweep(Xw,2,prop_mu)
#         new_cov <- crossprod(centred * sqrt(Ww)) / sum(Ww)
#         new_cov <- new_cov * cov_inflate
#       }
# 
#       ## ESS-weighted EMA -------------------------------------------------
#       ess_ratio <- ess_now / (ess_gate * d)
#       lambda_eff <- lambda_ema * min(1, ess_ratio)
#       prop_cov <- lambda_eff * new_cov + (1 - lambda_eff) * prev_cov
#       prev_cov <- prop_cov
# 
#       trace_now <- sum(diag(prop_cov))
#       arrow <- ifelse(trace_now > 1.001*prev_trace,"↑",
#                ifelse(trace_now < 0.999*prev_trace,"↓","–"))
#       prev_trace <- trace_now
# 
#       ## df auto-tune -----------------------------------------------------
#       if (identical(use_t_df,"auto") && nrow(samples) > 200) {
#         w_tmp <- (w^0.5)/sum(w^0.5)
#         mu_w  <- colSums(samples*w_tmp)
#         z <- sweep(samples,2,mu_w)
#         k_star <- max(colSums(w_tmp*z^4)/(colSums(w_tmp*z^2)^2))
#         if (k_star > 3.3) {
#           current_df <- max(6, 4 + 18/(k_star-3)); use_t_df <- current_df
#         } else current_df <- Inf
#         kf <- make_k(current_df); rmv_fn<-kf$rmv; dmv_fn<-kf$dmv
#       }
# 
#       mix_log <- function(x){
#         d1 <- dmv_fn(x, prop_mu, prop_cov)
#         d2 <- base_dmv(x, base_mu, base_cov)
#         m  <- pmax(d1,d2)
#         log(mix_w*exp(d1-m) + (1-mix_w)*exp(d2-m)) + m
#       }
# 
#       if (nrow(samples)>0)
#         log_qsum <- logplus(log_qsum, mix_log(samples))
# 
#       ## draw -------------------------------------------------------------
#       sel <- runif(n_batch) < mix_w
#       x_new <- rbind(rmv_fn(sum(sel),        prop_mu, prop_cov),
#                      rmv_fn(n_batch-sum(sel),base_mu, base_cov))[sample.int(n_batch),]
# 
#       if(cores>1) parallel::clusterExport(cl,"x_new",envir=environment())
# 
#       lp_new <- unlist(flexlapplytext(cl, seq_len(nrow(x_new)),
#                      fn="function(i){parlp(x_new[i,])}", cores = cores))
#       lp_new[is.na(lp_new)] <- -1e200
#       if(all(lp_new< -1e100)) stop("All lp = -Inf.")
# 
#       log_mix_new <- mix_log(x_new)
# 
#       ## store batch in window -------------------------------------------
#       win_ptr <- (win_ptr %% window_size) + 1L
#       log_w_batch <- lp_new - log_mix_new
#       w_batch <- exp(log_w_batch - max(log_w_batch))
#       win_x[[win_ptr]] <- x_new
#       win_w[[win_ptr]] <- w_batch / sum(w_batch)
# 
#       ## global bookkeeping ---------------------------------------------
#       samples      <- rbind(samples, x_new)
#       log_p_target <- c(log_p_target, lp_new)
#       log_qsum     <- c(log_qsum,    log_mix_new)
# 
#       T_props <- T_props + 1L
#       log_mix_all <- log_qsum - log(T_props)
#       log_w <- log_p_target - log_mix_all
#       log_w <- log_w - max(log_w)
#       w     <- pmax(exp(log_w), 0)
#       ess_now <- ess(w)
# 
#       if(verbose)
#         message(sprintf("iter %2d | Σ %s %.4f | λ %.3f | ESS %.1f | df %s",
#               iter, arrow, trace_now, lambda_eff, ess_now,
#               if(is.finite(current_df)) sprintf('%.1f',current_df) else "Inf"))
# 
#       if(ess_now >= target_ess) break
#     }
#   }, interrupt=function(e){
#        interrupted <<- TRUE
#        if(verbose) message("⏹  interrupted – returning current state.")
#   })
# 
#   if(!length(w)) stop("No samples collected before interrupt.")
#   w_norm <- w/sum(w); idx_post <- stratified_rs(w_norm,finishsamples)
# 
#   list(theta        = samples[idx_post,,drop=FALSE],
#        lpsamples    = log_p_target[idx_post],
#        weights      = rep(1/finishsamples,finishsamples),
#        full_theta   = samples,
#        full_weights = w_norm,
#        ess          = ess_now,
#        interrupted  = interrupted,
#        mean         = as.numeric(diagis::weighted_mean(samples,w_norm)),
#        covariance   = diagis::weighted_var(samples,w_norm),
#        df_used      = current_df)
# }



# # ======================================================================
# #  basic_is()  –  one-shot Laplace proposal with cluster + diagnostics
# # ======================================================================
# basic_is <- function(parlp,
#   mu_hat,
#   Sigma_hat,
#   cl,                       # parallel cluster for log-p
#   n_batch       = 5000,
#   target_ess    = 100,      # kept for API compatibility (unused)
#   max_iter      = 1,        # "
#   scale_init    = 3,        # "
#   cov_inflate   = 1.3,      # "
#   ridge         = 1e-8,
#   finishsamples = 1000,     # size of equal-weight resample
#   use_t_df      = NULL,     # "
#   power_upd     = 0.8,      # "
#   verbose       = TRUE,
#   diag_plots    = TRUE) {   # show weight_plot diagnostics
# 
#   ## ---- dependencies ----------------------------------------------------
#   for (pkg in c("mvtnorm", "diagis", "gridExtra", "ggplot2"))
#     if (!requireNamespace(pkg, quietly = TRUE))
#       stop(sprintf("Install the '%s' package.", pkg))
# 
#   ## ---- safe Σ ----------------------------------------------------------
#   safe_pd <- function(S, eps = ridge) {
#     S2 <- S + diag(eps, nrow(S))
#     while (any(eigen(S2, symmetric = TRUE,
#                      only.values = TRUE)$values <= 0))
#       S2 <- S2 + diag(eps, nrow(S))
#     S2
#   }
#   Sigma_hat <- safe_pd(Sigma_hat*scale_init, ridge)
# 
#   ## ---- draw ------------------------------------------------------------
#   theta <- mvtnorm::rmvnorm(n_batch, mu_hat, Sigma_hat)
# 
#   ## ---- log densities ---------------------------------------------------
#   if (!is.null(cl) && length(cl) > 1) {
#     parallel::clusterExport(cl, "theta", envir = environment())
#     log_p <- unlist(parallel::parLapply(
#       cl, seq_len(n_batch), function(i) parlp(theta[i, ])))
#   } else {
#     log_p <- vapply(seq_len(n_batch), \(i) parlp(theta[i, ]), numeric(1))
#   }
#   log_q <- mvtnorm::dmvnorm(theta, mu_hat, Sigma_hat, log = TRUE)
# 
#   ## ---- weights ---------------------------------------------------------
#   log_w <- log_p - log_q
#   log_w <- log_w - max(log_w)
#   w_raw  <- exp(log_w)
#   w_norm <- w_raw / sum(w_raw)
# 
#   ## ---- diagnostics -----------------------------------------------------
#   ess_val <- diagis::ess(w_raw)
#   if (verbose) {
#     cat("------------------------------------------------------------\n")
#     cat("Plain importance-sampling diagnostics\n")
#     cat("------------------------------------------------------------\n")
#     cat(sprintf("Samples drawn : %d\n", n_batch))
#     cat(sprintf("ESS           : %.1f (%.1f%% of n)\n",
#                 ess_val, 100 * ess_val / n_batch))
#     cat(sprintf("max weight    : %.3g\n", max(w_norm)))
#     cat(sprintf("min weight    : %.3g\n", min(w_norm)))
#     cat("------------------------------------------------------------\n\n")
#   }
# 
#   if (diag_plots) {
#     wg <- diagis::weight_plot(w_raw)
#     gridExtra::grid.arrange(
#       wg,
#       top = grid::textGrob(
#         sprintf("Weight diagnostics – ESS %.1f / %d",
#                 ess_val, n_batch),
#         gp = grid::gpar(fontface = "bold", fontsize = 14)))
#   }
# 
#   ## ---- summaries -------------------------------------------------------
#   mu_post <- as.numeric(diagis::weighted_mean(theta, w_norm))
#   Sigma_post <- diagis::weighted_var(theta, w_norm)
# 
#   ## equal-weight resample
#   stratified_rs <- function(w, N) {
#     cs <- cumsum(w / sum(w))
#     findInterval((runif(1) + 0:(N - 1)) / N, cs) + 1L
#   }
#   idx_eq <- stratified_rs(w_norm, finishsamples)
# 
#   list(
#     theta        = theta[idx_eq, , drop = FALSE],        # equal-weight draws
#     lpsamples    = log_p[idx_eq],
#     weights      = rep(1 / finishsamples, finishsamples),
#     full_theta   = theta,
#     full_weights = w_norm,
#     ess          = ess_val,
#     interrupted  = FALSE,
#     mean         = mu_post,
#     covariance   = Sigma_post,
#     df_used      = Inf                                   # Gaussian proposal
#   )
# }

# imis_is <- function(parlp,
#   mu_hat,
#   Sigma_hat,
#   cl,                          # parallel cluster for log-p
#   n_batch       = 1000,        # draws *per* component
#   target_ess    = 100,         # stop once ESS ≥ this
#   max_iter      = 10,          # max extra mixture components
#   scale_init    = 1.0,         # inflate Hessian for comp-0
#   tail_scale    = 4.0,         # SD of the added “bump” components
#   ridge         = 1e-8,
#   finishsamples = 1000,        # equal-weight resample
#   use_t_df      = NULL,        # kept for API compatibility
#   power_upd     = 0.8,         # unused (compat.)
#   verbose       = TRUE,
#   diag_plots    = TRUE) {
# 
#   for (pkg in c("mvtnorm", "diagis", "ggplot2", "gridExtra"))
#     if (!requireNamespace(pkg, quietly = TRUE))
#       stop(sprintf("Install '%s' first.", pkg))
# 
#   ## ---- helpers ----------------------------------------------------------
#   safe_pd <- function(S, eps = ridge) {
#     S2 <- S + diag(eps, nrow(S))
#     while (any(eigen(S2, symmetric = TRUE,
#                      only.values = TRUE)$values <= 0))
#       S2 <- S2 + diag(eps, nrow(S))
#     S2
#   }
#   Sigma_hat <- safe_pd(Sigma_hat)
# 
#   logplus <- function(a, b) {              # log(exp(a)+exp(b)) safely
#     idx <- a > b
#     r   <- numeric(length(a))
#     r[idx]  <- a[idx] + log1p(exp(b[idx]-a[idx]))
#     r[!idx] <- b[!idx] + log1p(exp(a[!idx]-b[!idx]))
#     r
#   }
# 
#   stratified_rs <- function(w, N) {
#     cs <- cumsum(w / sum(w))
#     findInterval((runif(1)+0:(N-1))/N, cs) + 1L
#   }
# 
#   ess <- function(w) diagis::ess(w)
# 
#   ## ---- component storage ------------------------------------------------
#   comp_mu  <- list(mu_hat)
#   comp_cov <- list(Sigma_hat * scale_init^2)
#   T_comp   <- 1L
# 
#   ## ---- global particle store -------------------------------------------
#   samples  <- matrix(0, 0, length(mu_hat))
#   log_p    <- log_q_mix <- numeric(0)
#   w_raw    <- numeric(0); ess_now <- 0
# 
#   repeat {
# 
#     ## ---- draw new batch from *current* mixture -------------------------
#     draw_from_mixture <- function(n) {
#       if (T_comp == 1L) {
#         mvtnorm::rmvnorm(n, comp_mu[[1]], comp_cov[[1]])
#       } else {
#         which_comp <- sample.int(T_comp, n, replace = TRUE)
#         do.call(rbind, lapply(seq_len(T_comp), function(k) {
#           m <- sum(which_comp == k)
#           if (m)
#             mvtnorm::rmvnorm(m, comp_mu[[k]], comp_cov[[k]])
#         }))
#       }
#     }
#     x_new <- draw_from_mixture(n_batch)
# 
#     ## ---- log target ----------------------------------------------------
#     if (!is.null(cl) && length(cl) > 1) {
#       parallel::clusterExport(cl, "x_new", envir = environment())
#       log_p_new <- unlist(parallel::parLapply(
#         cl, seq_len(nrow(x_new)), function(i) parlp(x_new[i, ])))
#     } else {
#       log_p_new <- vapply(seq_len(nrow(x_new)),
#                           \(i) parlp(x_new[i, ]), numeric(1))
#     }
# 
#     ## ---- log proposal for ALL comps -----------------------------------
#     log_q_new <- rep(-Inf, n_batch)
#     for (k in seq_len(T_comp)) {
#       log_q_new <- logplus(
#         log_q_new,
#         mvtnorm::dmvnorm(x_new, comp_mu[[k]], comp_cov[[k]], log = TRUE)
#       )
#     }
#     log_q_new <- log_q_new - log(T_comp)   # equal-weight mixture
# 
#     ## ---- append to global stores --------------------------------------
#     samples    <- rbind(samples, x_new)
#     log_p      <- c(log_p,      log_p_new)
#     log_q_mix  <- c(log_q_mix,  log_q_new)
# 
#     ## ---- weights / ESS -------------------------------------------------
#     log_w <- log_p - log_q_mix
#     log_w <- log_w - max(log_w)
#     w_raw <- exp(log_w)
#     ess_now <- ess(w_raw)
# 
#     if (verbose)
#       message(sprintf("iter %2d | comp %d | n = %d | ESS = %.1f",
#               T_comp - 1, T_comp, nrow(samples), ess_now))
# 
#   if (diag_plots) {
#     wg <- diagis::weight_plot(w_raw)
#     gridExtra::grid.arrange(
#       wg,
#       top = grid::textGrob(
#         sprintf("Weight diagnostics – ESS %.1f / %d",
#                 ess_now, n_batch),
#         gp = grid::gpar(fontface = "bold", fontsize = 14)))
#   }
# 
#     ## ---- stopping rule --------------------------------------------------
#     if (ess_now >= target_ess || T_comp > max_iter)
#       break
# 
#     ## ---- add new component at highest weight ---------------------------
#     best_idx <- which.max(w_raw)
#     comp_mu[[T_comp + 1L]]  <- samples[best_idx, ]
#     comp_cov[[T_comp + 1L]] <- diag(tail_scale^2, nrow(Sigma_hat))
#     T_comp <- T_comp + 1L
#   }
# 
#   ## ---- output summaries -------------------------------------------------
#   w_norm <- w_raw / sum(w_raw)
#   mu_post <- as.numeric(diagis::weighted_mean(samples, w_norm))
#   Sigma_post <- diagis::weighted_var(samples, w_norm)
# 
#   idx_eq <- stratified_rs(w_norm, finishsamples)
# 
#   list(
#     theta        = samples[idx_eq, , drop = FALSE],
#     lpsamples    = log_p[idx_eq],
#     weights      = rep(1 / finishsamples, finishsamples),
#     full_theta   = samples,
#     full_weights = w_norm,
#     ess          = ess_now,
#     interrupted  = FALSE,
#     mean         = mu_post,
#     covariance   = Sigma_post,
#     df_used      = Inf         # all components Gaussian here
#   )
# }
