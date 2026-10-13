# 03-dependencia.R — autocorrelacion y pruebas de dependencia.
# Solo define funciones: no ejecuta nada al hacer source().

# ACF muestral r_1..r_m: r_k = sum_{t=1}^{n-k} (y_t - ybar)(y_{t+k} - ybar) / sum_t (y_t - ybar)^2.
# Se divide por n (no por n-k) en gamma(k), como hace acf(); asi la matriz de autocovarianzas es semidefinida positiva.
acf_muestral <- function(y, m) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) > 2, is.numeric(m), length(m) == 1, m >= 1, m < length(y))
  n <- length(y)
  yc <- y - mean(y)
  g0 <- sum(yc^2)
  vapply(seq_len(m), function(k) sum(yc[1:(n - k)] * yc[(1 + k):n]) / g0, numeric(1))
}
