# 02-operadores.R — operador de rezago B, diferencias y raices de polinomios en B.
# Solo define funciones: no ejecuta nada al hacer source().

# Aplica (1-B)^d (1-B^s)^D a y, sin usar diff().
# (1-B) y_t = y_t - y_{t-1}; (1-B^s) y_t = y_t - y_{t-s}. Los operadores conmutan, asi que el orden no importa.
# Cada aplicacion de lag L pierde L observaciones al inicio, por eso el resultado tiene
# length(y) - d - D*s valores (no se rellena con NA).
diferenciar <- function(y, d = 0, D = 0, s = 1) {
  stopifnot(is.numeric(y), !anyNA(y), is.numeric(d), length(d) == 1, d >= 0, d == round(d),
            is.numeric(D), length(D) == 1, D >= 0, D == round(D),
            is.numeric(s), length(s) == 1, s >= 1, s == round(s))
  if (length(y) - d - D * s < 1) stop("La serie es demasiado corta para esas diferencias.")
  for (i in seq_len(d)) y <- y[-1] - y[-length(y)]
  for (i in seq_len(D)) y <- y[(s + 1):length(y)] - y[1:(length(y) - s)]
  y
}

# Multiplica dos polinomios en B dados por sus coeficientes (de menor a mayor grado): convolucion.
multiplicar_polinomios <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1)
  for (i in seq_along(a)) {
    idx <- i:(i + length(b) - 1)
    out[idx] <- out[idx] + a[i] * b
  }
  out
}

# Coeficientes c1..cp del polinomio c(B) = (1-B)^d (1-B^s)^D = 1 + c1 B + ... + cp B^p.
coef_diferencia <- function(d = 0, D = 0, s = 1) {
  stopifnot(is.numeric(d), d >= 0, is.numeric(D), D >= 0, is.numeric(s), s >= 1)
  p <- 1
  for (i in seq_len(d)) p <- multiplicar_polinomios(p, c(1, -1))
  for (i in seq_len(D)) p <- multiplicar_polinomios(p, c(1, rep(0, s - 1), -1))
  p[-1]
}

# Raices de c(B) = 1 + c1 B + ... + cp B^p y su modulo. coef = c(c1, ..., cp).
# polyroot() pide los coeficientes de menor a mayor grado, incluido el 1 del termino constante.
# Si el ultimo coeficiente es 0 el grado real es menor: se recortan los ceros finales.
# Modulo > 1: raiz fuera del circulo unitario; modulo = 1: raiz unitaria (el factor elimina un componente determinista).
polinomio_rezago <- function(coef) {
  stopifnot(is.numeric(coef), !anyNA(coef))
  while (length(coef) > 0 && coef[length(coef)] == 0) coef <- coef[-length(coef)]
  if (length(coef) == 0) return(tibble::tibble(re = numeric(0), im = numeric(0), modulo = numeric(0)))
  r <- polyroot(c(1, coef))
  tibble::tibble(re = Re(r), im = Im(r), modulo = Mod(r))
}

# Tabla de la seccion 4c: n, media, varianza muestral y r1 de la serie y de cada candidata.
# candidatas: lista de list(nombre, log, d, D, s); log = TRUE aplica logaritmo natural antes de diferenciar.
evaluar_candidatas <- function(y, candidatas) {
  stopifnot(is.numeric(y), !anyNA(y), is.list(candidatas))
  filas <- lapply(candidatas, function(cd) {
    z <- aplicar_candidata(y, cd)
    tibble::tibble(candidata = cd$nombre, n = length(z), media = mean(z),
                   varianza = stats::var(z), r1 = acf_muestral(z, 1))
  })
  dplyr::bind_rows(filas)
}

# Aplica una candidata list(nombre, log, d, D, s, tend) a y: primero logaritmo (si log), luego se quita
# una tendencia polinomial de grado `tend` (si no es NULL), luego las diferencias. Una sola definicion
# para que la tabla de la seccion 4 y el analisis de la seccion 5 usen exactamente la misma transformacion.
aplicar_candidata <- function(y, cd) {
  z <- y
  if (isTRUE(cd$log)) {
    if (any(y <= 0)) stop("log requiere y > 0.")
    z <- log(y)
  }
  if (!is.null(cd$tend)) z <- residuos_tendencia(z, cd$tend)$residuos
  diferenciar(z, cd$d, cd$D, cd$s)
}

# Residuos de ajustar por minimos cuadrados una tendencia polinomial de grado `grado` en t = 1..n.
# Se resuelven las ecuaciones normales (X'X) b = X'y; sin lm(). Con un intercepto en el modelo,
# la media de los residuos es 0 por construccion (la tendencia absorbe el nivel).
residuos_tendencia <- function(y, grado = 1) {
  stopifnot(is.numeric(y), !anyNA(y), is.numeric(grado), length(grado) == 1, grado >= 0, grado == round(grado),
            length(y) > grado + 1)
  t <- seq_along(y)
  X <- outer(t, 0:grado, "^")
  b <- solve(crossprod(X), crossprod(X, y))
  list(coeficientes = as.numeric(b), ajustado = as.numeric(X %*% b), residuos = as.numeric(y - X %*% b))
}

# Identidad Var(grad Y) = 2 gamma(0) - 2 gamma(1), con gamma muestrales (divisor n).
# Devuelve ambos lados: son casi iguales; la diferencia viene solo de los efectos de borde
# (la varianza de grad Y usa n-1 terminos y su propia media).
identidad_var_diferencia <- function(y) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) > 3)
  n <- length(y)
  yc <- y - mean(y)
  g0 <- sum(yc^2) / n
  g1 <- sum(yc[-1] * yc[-n]) / n
  tibble::tibble(var_diferencia = stats::var(diferenciar(y, 1, 0, 1)),
                 dos_g0_menos_dos_g1 = 2 * g0 - 2 * g1)
}
