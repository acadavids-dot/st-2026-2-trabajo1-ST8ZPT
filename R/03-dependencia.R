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

# PACF muestral por la recursion de Durbin-Levinson, a partir de r_1..r_m.
# phi_kk = (r_k - sum_{j<k} phi_{k-1,j} r_{k-j}) / (1 - sum_{j<k} phi_{k-1,j} r_j): es el ultimo coeficiente del
# mejor predictor lineal de y_t con k rezagos. Los demas coeficientes se actualizan con
# phi_kj = phi_{k-1,j} - phi_kk * phi_{k-1,k-j}.
pacf_dl <- function(r) {
  stopifnot(is.numeric(r), !anyNA(r), length(r) >= 1)
  m <- length(r)
  pacf_k <- numeric(m)
  phi <- r[1]
  pacf_k[1] <- r[1]
  if (m >= 2) {
    for (k in 2:m) {
      prev <- phi
      num <- r[k] - sum(prev * r[(k - 1):1])
      den <- 1 - sum(prev * r[1:(k - 1)])
      pkk <- num / den
      phi <- c(prev - pkk * rev(prev), pkk)
      pacf_k[k] <- pkk
    }
  }
  pacf_k
}

# Banda de Bartlett para la ACF en cada rezago h: +-1.96 * sqrt((1 + 2 sum_{i<h} r_i^2) / n).
# Responde "es r_h distinta de cero SABIENDO que ya hay autocorrelacion en rezagos menores?" (no asume ruido blanco).
# La banda de ruido blanco +-1.96/sqrt(n) responde "es esta serie ruido blanco?".
banda_bartlett <- function(r, n) {
  stopifnot(is.numeric(r), !anyNA(r), n > 1)
  m <- length(r)
  1.96 * sqrt((1 + 2 * cumsum(c(0, r^2))[1:m]) / n)
}

# Varianza de largo plazo con nucleo de Bartlett: sigma2_LP = g0 + 2 sum_{j=1}^{L} (1 - j/(L+1)) g_j.
# Los pesos decrecientes garantizan un resultado no negativo. g_j usa divisor n, como la ACF.
varianza_largo_plazo <- function(y, L) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) > 2, is.numeric(L), length(L) == 1, L >= 0, L < length(y))
  n <- length(y)
  yc <- y - mean(y)
  g <- vapply(0:L, function(j) sum(yc[1:(n - j)] * yc[(1 + j):n]) / n, numeric(1))
  g[1] + if (L >= 1) 2 * sum((1 - (1:L) / (L + 1)) * g[-1]) else 0
}

# Media y varianza de k bloques consecutivos de tamano (casi) igual (seccion 5c).
bloques_media_var <- function(y, k = 4) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) >= 2 * k)
  # `grupo` no se llama `bloque` para que tibble() no lo confunda con la columna `bloque` que define.
  grupo <- ceiling(seq_along(y) * k / length(y))
  tibble::tibble(bloque = seq_len(k),
                 n = as.integer(table(grupo)),
                 media = as.numeric(tapply(y, grupo, mean)),
                 varianza = as.numeric(tapply(y, grupo, stats::var)))
}

# Dos intervalos al 95 % para la media: error estandar ingenuo sigma/sqrt(n) (supone independencia)
# y error estandar de largo plazo sqrt(sigma2_LP / n), con L = floor(4 (n/100)^(2/9)).
ic_media <- function(y) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) > 10)
  n <- length(y)
  L <- floor(4 * (n / 100)^(2 / 9))
  mu <- mean(y)
  se_ing <- stats::sd(y) / sqrt(n)
  se_lp <- sqrt(varianza_largo_plazo(y, L) / n)
  tibble::tibble(n = n, L = L, media = mu, se_ingenuo = se_ing, se_largo_plazo = se_lp,
                 ic_ingenuo_inf = mu - 1.96 * se_ing, ic_ingenuo_sup = mu + 1.96 * se_ing,
                 ic_lp_inf = mu - 1.96 * se_lp, ic_lp_sup = mu + 1.96 * se_lp,
                 razon_se = se_lp / se_ing)
}

# Prueba de Ljung-Box con los seis elementos. Referida a chi2_m: si y son errores o residuos con
# regresores deterministas, los coeficientes NO restan grados de libertad (ver enunciado).
# que: descripcion de la serie analizada para redactar la lectura.
ljung_box <- function(y, m, que = "la serie", alpha = 0.05) {
  stopifnot(is.numeric(y), !anyNA(y), is.numeric(m), length(m) == 1, m >= 1, m < length(y))
  n <- length(y)
  r <- acf_muestral(y, m)
  Q <- n * (n + 2) * sum(r^2 / (n - seq_len(m)))
  crit <- stats::qchisq(1 - alpha, m)
  nueva_prueba(
    nombre = paste0("Prueba de Ljung-Box sobre ", que, " (m = ", m, ")"),
    h0 = "rho_1 = ... = rho_m = 0 (no hay autocorrelacion hasta el rezago m)",
    h1 = "algun rho_k distinto de 0 con k <= m",
    estadistico = "Q",
    formula = "n(n+2) * sum_{k=1}^{m} r_k^2/(n-k)",
    distribucion = paste0("chi-cuadrado con ", m, " grados de libertad"),
    region = paste0("Q > ", formatear_num(crit)),
    valor_critico = crit,
    observado = Q,
    valor_p = stats::pchisq(Q, m, lower.tail = FALSE),
    lectura_rechaza = paste0("Hay autocorrelacion en ", que, " hasta el rezago ", m, ": no se comporta como ruido blanco."),
    lectura_no_rechaza = paste0("No hay evidencia de autocorrelacion en ", que, " hasta el rezago ", m, ": es compatible con ruido blanco."),
    alpha = alpha
  )
}
