library(dplyr)
library(RANN)

connect_trees <- function(
    inventur_alt_path,
    inventur_neu_path,
    lower_limit = 0.4,
    upper_limit = 0.7
) {
  
  # --------------------------------------------------------------
  # Hilfsfunktionen
  # --------------------------------------------------------------
  center_point <- function(mat) {
    c(mean(range(mat[, 1])), mean(range(mat[, 2])))
  }
  shift_mat <- function(mat, shift) {
    mat + matrix(shift, nrow = nrow(mat), ncol = 2, byrow = TRUE)
  }
  rotate_mat <- function(mat, center, angle_deg) {
    ang <- angle_deg * pi / 180
    R   <- matrix(c(cos(ang), -sin(ang),
                    sin(ang),  cos(ang)), 2, 2)
    t(R %*% t(mat - matrix(center, nrow = nrow(mat), ncol = 2, byrow = TRUE))) +
      matrix(center, nrow = nrow(mat), ncol = 2, byrow = TRUE)
  }
  mean_nn_dist <- function(query, data_mat) {
    mean(nn2(data = data_mat, query = query, k = 1)$nn.dists)
  }
  compute_weights <- function(mat) {
    ctr <- center_point(mat)
    d   <- sqrt(rowSums((mat - matrix(ctr, nrow(mat), 2, byrow = TRUE))^2))
    if (max(d) == 0) return(rep(1, length(d)))
    pmax(1 - d / max(d), 0)
  }
  search_absolute_angle <- function(mat, center, data_mat, angles) {
    errs <- vapply(
      angles,
      function(a) mean_nn_dist(rotate_mat(mat, center, a), data_mat),
      numeric(1)
    )
    best <- angles[which.min(errs)]
    list(angle = best,
         mat   = rotate_mat(mat, center, best),
         error = min(errs))
  }
  refine_rotation <- function(mat, center, data_mat, window, n = 181) {
    angles <- seq(window[1], window[2], length.out = n)
    errs <- vapply(
      angles,
      function(a) mean_nn_dist(rotate_mat(mat, center, a), data_mat),
      numeric(1)
    )
    best <- angles[which.min(errs)]
    list(angle = best,
         mat   = rotate_mat(mat, center, best),
         error = min(errs))
  }
  optimise_shift <- function(mat, data_mat, limit) {
    fn <- function(shift) {
      mean_nn_dist(shift_mat(mat, shift), data_mat)
    }
    opt <- optim(
      par     = c(0, 0),
      fn      = fn,
      method  = "L-BFGS-B",
      lower   = rep(-limit, 2),
      upper   = rep(limit, 2)
    )
    list(shift = opt$par,
         mat   = shift_mat(mat, opt$par),
         value = opt$value)
  }
  optimise_weighted_transform <- function(mat, data_mat, center, weights) {
    fn <- function(par) {
      rot <- rotate_mat(mat, center, par[1])
      trn <- shift_mat(rot, par[2:3])
      dists <- nn2(data = data_mat, query = trn, k = 1)$nn.dists
      mean(weights * dists)
    }
    optim(par = c(0, 0, 0), fn = fn, method = "BFGS")
  }
  rotation_matrix_rad <- function(angle_rad) {
    matrix(c(cos(angle_rad), -sin(angle_rad),
             sin(angle_rad),  cos(angle_rad)), 2, 2)
  }
  
  # --------------------------------------------------------------
  # 1. Daten einlesen und filtern
  # --------------------------------------------------------------
  Inventur_alt_raw <- read.csv(inventur_alt_path)
  Inventur_neu_raw <- read.csv(inventur_neu_path)
  
  Inventur_alt <- transmute(Inventur_alt_raw, X = X, Y = Y, BHD = BHD / 100)
  Inventur_neu <- transmute(Inventur_neu_raw, X = X, Y = Y, BHD = BHD / 100)
  
  Inventur_alt_f <- filter(Inventur_alt, between(BHD, lower_limit, upper_limit))
  Inventur_neu_f <- filter(Inventur_neu, between(BHD, lower_limit, upper_limit))
  
  X_mat <- as.matrix(select(Inventur_alt_f, X, Y))
  Y_mat <- as.matrix(select(Inventur_neu_f, X, Y))
  
  # --------------------------------------------------------------
  # 2. Grobe Translation und Rotation
  # --------------------------------------------------------------
  center_X <- center_point(X_mat)
  center_Y <- center_point(Y_mat)
  
  translation_vec0 <- center_X - center_Y
  Y_trans_mat <- shift_mat(Y_mat, translation_vec0)
  
  coarse_angle <- search_absolute_angle(Y_trans_mat, center_X, X_mat, seq(0, 360, by = 1))
  angle0 <- coarse_angle$angle
  coords_rot_mat <- coarse_angle$mat
  
  ref1 <- refine_rotation(coords_rot_mat, center_X, X_mat, window = c(-2, 2), n = 181)
  delta1 <- ref1$angle
  coords_rot_mat <- ref1$mat
  
  ref2 <- refine_rotation(coords_rot_mat, center_X, X_mat, window = c(-0.5, 0.5), n = 181)
  delta2 <- ref2$angle
  coords_rot_mat <- ref2$mat
  
  # --------------------------------------------------------------
  # 3. Verschiebung optimieren (ungewichtet)
  # --------------------------------------------------------------
  extent_x <- diff(range(Y_mat[, 1]))
  extent_y <- diff(range(Y_mat[, 2]))
  shift_limit <- max(extent_x, extent_y) / 10
  
  shift_opt <- optimise_shift(coords_rot_mat, X_mat, shift_limit)
  translation_vec1 <- shift_opt$shift
  coords_shifted   <- shift_opt$mat
  
  # --------------------------------------------------------------
  # 4. Gewichtete Feinoptimierung
  # --------------------------------------------------------------
  weights <- compute_weights(coords_shifted)
  opt_res <- optimise_weighted_transform(coords_shifted, X_mat, center_X, weights)
  
  # --------------------------------------------------------------
  # 5. Affine Transformation y' = A y + b
  # --------------------------------------------------------------
  theta0_rad  <- (angle0 + delta1 + delta2) * pi / 180
  theta_opt   <- opt_res$par[1] * pi / 180
  A           <- rotation_matrix_rad(theta_opt) %*% rotation_matrix_rad(theta0_rad)
  
  c_vec  <- matrix(center_X, ncol = 1)
  t0_vec <- matrix(translation_vec0, ncol = 1)
  t1_vec <- matrix(translation_vec1, ncol = 1)
  s_vec  <- matrix(opt_res$par[2:3], ncol = 1)
  
  b_vec <- (diag(2) - A) %*% c_vec + A %*% t0_vec + rotation_matrix_rad(theta_opt) %*% t1_vec + s_vec
  
  # --------------------------------------------------------------
  # 6. Auf komplette Inventur_neu anwenden
  # --------------------------------------------------------------
  Y_full_mat    <- as.matrix(select(Inventur_neu, X, Y))
  Y_full_trans  <- Y_full_mat %*% t(A) + matrix(rep(t(b_vec), each = nrow(Y_full_mat)), ncol = 2, byrow = FALSE)
  
  Inventur_neu_trans <- Inventur_neu_raw
  Inventur_neu_trans$X <- Y_full_trans[, 1]
  Inventur_neu_trans$Y <- Y_full_trans[, 2]
  
  list(
    Inventur_alt_filtered = Inventur_alt_f,
    Inventur_neu_filtered = Inventur_neu_f,
    Inventur_neu_transformed = Inventur_neu_trans,
    params = list(
      A = A,
      b = as.numeric(b_vec),
      angle_components_deg = list(angle0 = angle0,
                                  delta1 = delta1,
                                  delta2 = delta2,
                                  optim_angle = opt_res$par[1]),
      shifts = list(
        translation_vec0 = translation_vec0,
        translation_vec1 = translation_vec1,
        final_shift      = opt_res$par[2:3]
      )
    )
  )
}




# Beispielaufruf:
result <- connect_trees(
  inventur_alt_path = "path/to/Inventur_alt.csv",
  inventur_neu_path = "path/to/Inventur_neu.csv"
)

transformed <- result$Inventur_neu_transformed
head(transformed)

# Optional: plotten
plot(transformed$X, transformed$Y, pch = 19, col = "black",cex = 0.7)
points(result$Inventur_alt_filtered$X, result$Inventur_alt_filtered$Y, pch = 20, col = "red", cex = 0.5)
 