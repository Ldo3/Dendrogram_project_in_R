# Required Libraries
library(MASS)
library(matrixStats)
library(mclust)
library(ggplot2)
library(cluster)

# Weak merging function
merging_of_atoms_weak_R <- function(pi, mus, Sigmas, l_cur) {
     k <- length(pi)
     l_new <- l_cur
     
     dist_mat <- outer(1:k, 1:k, Vectorize(function(i, j) {
          if (i == j) return(1e4)
          pi[i] * sum((mus[i, ] - mus[j, ])^2)
     }))
     
     min_idx <- which(dist_mat == min(dist_mat), arr.ind = TRUE)[1, ]
     merge_indices <- sort(min_idx)
     
     i <- merge_indices[1]
     j <- merge_indices[2]
     
     merge_set <- list(union(l_cur[[i]], l_cur[[j]]))
     
     pi_merge <- pi[i] + pi[j]
     mus_merge <- (pi[i] * mus[i, ] + pi[j] * mus[j, ]) / pi_merge
     
     diff_i <- mus[i, ] - mus_merge
     diff_j <- mus[j, ] - mus_merge
     Sigmas_merge <- (pi[i]/pi_merge)*(Sigmas[[i]] + tcrossprod(diff_i)) +
          (pi[j]/pi_merge)*(Sigmas[[j]] + tcrossprod(diff_j))
     
     pi_new <- pi
     mus_new <- mus
     Sigmas_new <- Sigmas
     
     pi_new[i] <- pi_merge
     mus_new[i, ] <- mus_merge
     Sigmas_new[[i]] <- Sigmas_merge
     
     pi_new <- pi_new[-j]
     mus_new <- mus_new[-j, , drop = FALSE]
     Sigmas_new <- Sigmas_new[-j]
     
     l_new[[i]] <- merge_set[[1]]
     l_new <- l_new[-j]
     
     list(pi = pi_new, mus = mus_new, Sigmas = Sigmas_new, k = k - 1, l = l_new, d_min = dist_mat[i, j])
}

# Strong merging function
merging_of_atoms_strong_R <- function(pi, mus, l_cur) {
     k <- length(pi)
     l_new <- l_cur
     
     har_mean <- outer(pi, pi, function(a, b) (a * b) / (a + b))
     dist_mat <- outer(1:k, 1:k, Vectorize(function(i, j) {
          if (i == j) return(1e4)
          sqrt(har_mean[i, j] * sum((mus[i, ] - mus[j, ])^2))
     }))
     
     min_idx <- which(dist_mat == min(dist_mat), arr.ind = TRUE)[1, ]
     merge_indices <- sort(min_idx)
     
     i <- merge_indices[1]
     j <- merge_indices[2]
     
     merge_set <- list(union(l_cur[[i]], l_cur[[j]]))
     
     pi_merge <- pi[i] + pi[j]
     mus_merge <- (pi[i] * mus[i, ] + pi[j] * mus[j, ]) / pi_merge
     
     pi_new <- pi
     mus_new <- mus
     pi_new[i] <- pi_merge
     mus_new[i, ] <- mus_merge
     
     pi_new <- pi_new[-j]
     mus_new <- mus_new[-j, , drop = FALSE]
     
     l_new[[i]] <- merge_set[[1]]
     l_new <- l_new[-j]
     
     list(pi = pi_new, mus = mus_new, k = k - 1, l = l_new, d_min = dist_mat[i, j])
}

# DSC
Dendrogram_Inferred_Clustering <- function(pi, mus) {
     
     n <- length(pi)   # number of atoms
     
     # Preallocate lists
     result_list_group_merge <- vector("list", n)
     result_pi_merge <- vector("list", n)
     result_mus_merge <- vector("list", n)
     K_list <- numeric(n)
     d_list <- numeric(n)
     
     # Initial groups: each atom is its own set
     missing <- seq_len(n)
     l_cur <- lapply(missing, function(k) { c(k) })  # ensures integers
     
     pi_cur <- pi
     mus_cur <- mus
     
     result_pi_merge[[1]] <- pi_cur
     result_list_group_merge[[1]] <- l_cur
     result_mus_merge[[1]] <- mus_cur
     
     K_list[1] <- n
     
     # Iterative merging
     for (i in 2:n) {
          res <- merging_of_atoms_strong(pi_cur, mus_cur, l_cur)
          pi_new <- res$pi_new
          mus_new <- res$mus_new
          k_new <- res$k_new
          l_new <- res$l_new
          d_min <- res$d_min
          
          result_pi_merge[[i]] <- pi_new
          result_list_group_merge[[i]] <- l_new
          result_mus_merge[[i]] <- mus_new
          K_list[i] <- k_new
          d_list[i - 1] <- d_min
          
          # update for next iteration
          l_cur <- l_new
          pi_cur <- pi_new
          mus_cur <- mus_new
     }
     
     return(list(
          result_pi_merge = result_pi_merge,
          result_list_group_merge = result_list_group_merge,
          result_mus_merge = result_mus_merge,
          K_list = K_list,
          d_list = d_list
     ))
}

# DSC Location-Scale Gaussian
dsc_location_scale_gaussian_R <- function(pi, mus, Sigmas) {
     n <- length(pi)
     result_list_group_merge <- vector("list", n)
     result_pi_merge <- vector("list", n)
     result_mus_merge <- vector("list", n)
     result_Sigmas_merge <- vector("list", n)
     K_list <- numeric(n)
     d_list <- numeric(n)
     
     l_cur <- lapply(1:n, function(i) i)
     
     pi_cur <- pi
     mus_cur <- mus
     Sigmas_cur <- Sigmas
     
     result_pi_merge[[1]] <- pi_cur
     result_list_group_merge[[1]] <- l_cur
     result_mus_merge[[1]] <- mus_cur
     result_Sigmas_merge[[1]] <- Sigmas_cur
     K_list[1] <- n
     
     for (i in 2:n) {
          res <- merging_of_atoms_weak_R(pi_cur, mus_cur, Sigmas_cur, l_cur)
          pi_cur <- res$pi
          mus_cur <- res$mus
          Sigmas_cur <- res$Sigmas
          l_cur <- res$l
          
          result_pi_merge[[i]] <- pi_cur
          result_list_group_merge[[i]] <- l_cur
          result_mus_merge[[i]] <- mus_cur
          result_Sigmas_merge[[i]] <- Sigmas_cur
          K_list[i] <- res$k
          d_list[i - 1] <- res$d_min
     }
     
     list(pi = result_pi_merge, l = result_list_group_merge, mus = result_mus_merge,
          Sigmas = result_Sigmas_merge, K = K_list, d = d_list)
}

# Example GMM fit and usage with mclust
fit_gmm <- function(X, G = 3) {
     gmm_fit <- Mclust(X, G = G)
     list(pi = gmm_fit$parameters$pro,
          mus = t(gmm_fit$parameters$mean),
          Sigmas = lapply(1:G, function(i) gmm_fit$parameters$variance$sigma[,,i]),
          model = gmm_fit)
}

# Get ICL
get_icl <- function(gmm_model) {
     gmm_model$icl
}

# Full dendrogram plot with hierarchical merges
plot_dendrogram_R <- function(mus_merge, list_group_merge, d_list, K_bar, size = 3) {
     d_list_cumsum <- c(0, cumsum(d_list))[-(K_bar + 1)]
     
     plot(NULL, xlim = c(0, max(d_list_cumsum)), ylim = range(sapply(mus_merge, function(x) range(x[,1]))),
          xlab = "D", ylab = "mu[1]", main = "Dendrogram")
     
     for (i in seq_len(K_bar - 1)) {
          for (j in seq_len(K_bar - i)) {
               for (jj in seq_len(K_bar - i - 1)) {
                    if (all(list_group_merge[[i]][[j]] %in% list_group_merge[[i + 1]][[jj]])) {
                         y1 <- mus_merge[[i]][j, 1]
                         y2 <- mus_merge[[i + 1]][jj, 1]
                         x1 <- d_list_cumsum[i]
                         x2 <- d_list_cumsum[i + 1]
                         segments(x1, y1, x2, y1, col = 'gray')
                         segments(x2, y1, x2, y2, col = 'gray')
                    }
               }
          }
     }
     for (i in seq_len(K_bar)) {
          points(rep(d_list_cumsum[i], K_bar - i + 1), mus_merge[[i]][, 1], pch = 19, cex = size)
          text(rep(d_list_cumsum[i], K_bar - i + 1), mus_merge[[i]][, 1],
               labels = seq_len(nrow(mus_merge[[i]])) - 1, col = 'white', cex = 0.6, pos = 3)
     }
}