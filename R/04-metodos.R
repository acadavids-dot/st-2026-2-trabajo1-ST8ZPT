# 04-metodos.R — catalogo de metodos de pronostico, un metodo por funcion.
# Solo define funciones: no ejecuta nada al hacer source().
#
# Contrato de cada metodo (ajustar_xxx(y, ...)):
#   yhat        vector de largo length(y); yhat[t] usa SOLO y[1:(t-1)]; el calentamiento queda en NA.
#   pronosticar function(h) -> Y_{T+1}, ..., Y_{T+h}; solo usa el estado del ajuste, no objetos globales.
#   parametros  lista con las constantes usadas y los valores iniciales.

# Valida el vector de entrada comun a todos los metodos.
validar_y <- function(y, n_min) {
  if (!is.numeric(y)) stop("y debe ser numerico.")
  if (anyNA(y)) stop("y no puede tener NA.")
  if (length(y) < n_min) stop("y debe tener al menos ", n_min, " observaciones.")
  invisible(TRUE)
}

# Valida un parametro entero escalar dentro de [minimo, maximo].
validar_entero <- function(x, nombre, minimo, maximo = Inf) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != round(x) || x < minimo || x > maximo) {
    stop(nombre, " debe ser un entero entre ", minimo, " y ", maximo, ".")
  }
  invisible(TRUE)
}

# Media simple: el pronostico de t es el promedio de todo lo observado antes de t.
ajustar_media <- function(y) {
  validar_y(y, 2)
  n <- length(y)
  mu <- mean(y)
  # cumsum(y)[t-1] / (t-1) = promedio de y[1:(t-1)]
  yhat <- c(NA_real_, cumsum(y)[1:(n - 1)] / (1:(n - 1)))
  list(yhat = yhat,
       pronosticar = function(h) rep(mu, h),
       parametros = list(media = mu))
}

# Ingenuo: el pronostico de t es el ultimo valor observado.
ajustar_ingenuo <- function(y) {
  validar_y(y, 2)
  n <- length(y)
  ultimo <- y[n]
  list(yhat = c(NA_real_, y[-n]),
       pronosticar = function(h) rep(ultimo, h),
       parametros = list(ultimo = ultimo))
}

# Ingenuo estacional: el pronostico de t es el valor de hace s periodos.
ajustar_ingenuo_estacional <- function(y, s) {
  validar_entero(s, "s", 2)
  validar_y(y, s + 1)
  n <- length(y)
  ciclo <- y[(n - s + 1):n]  # ultimo ciclo completo, que se repite hacia adelante
  list(yhat = c(rep(NA_real_, s), y[1:(n - s)]),
       pronosticar = function(h) ciclo[((seq_len(h) - 1) %% s) + 1],
       parametros = list(s = s, ultimo_ciclo = ciclo))
}

# Media movil de orden k: promedio de las ultimas k observaciones; pronostico plano.
ajustar_mm <- function(y, k) {
  validar_entero(k, "k", 2)
  validar_y(y, k + 1)
  n <- length(y)
  acum <- c(0, cumsum(y))  # acum[t] = sum(y[1:(t-1)])
  yhat <- rep(NA_real_, n)
  # suma de y[(t-k):(t-1)] = acum[t] - acum[t-k]
  for (t in (k + 1):n) yhat[t] <- (acum[t] - acum[t - k]) / k
  nivel <- mean(y[(n - k + 1):n])
  list(yhat = yhat,
       pronosticar = function(h) rep(nivel, h),
       parametros = list(k = k, nivel = nivel))
}

# Matriz de diseno de las tendencias, filas t = 1..n.
#   lineal: 1, t | cuadratica: 1, t, t^2 | exponencial: 1, t (se ajusta log y)
#   estacional: 1, t y s-1 ficticias de suma cero: D_j = 1 en la estacion j, -1 en la estacion s, 0 en las demas
#   (asi el efecto de la estacion s es menos la suma de los demas y los efectos suman 0).
# La posicion 1 es la primera observacion de la serie.
matriz_diseno <- function(n, tipo = c("lineal", "cuadratica", "exponencial", "estacional"), s = NULL) {
  tipo <- match.arg(tipo)
  t <- seq_len(n)
  X <- switch(tipo,
    lineal = cbind(Intercepto = 1, t = t),
    exponencial = cbind(Intercepto = 1, t = t),
    cuadratica = cbind(Intercepto = 1, t = t, t2 = t^2),
    estacional = {
      pos <- ((t - 1) %% s) + 1
      D <- vapply(seq_len(s - 1), function(j) ifelse(pos == j, 1, ifelse(pos == s, -1, 0)), numeric(n))
      colnames(D) <- paste0("D", seq_len(s - 1))
      cbind(Intercepto = 1, t = t, D)
    })
  X
}

# MCO por descomposicion QR (mas estable que invertir X'X cuando hay t y t^2).
coef_mco <- function(X, y) as.numeric(qr.coef(qr(X), y))

# Nucleo comun de las tendencias. Para que yhat[t] use solo y[1:(t-1)], en cada t se reajusta la regresion
# con las primeras t-1 observaciones y se predice la fila t (ventana creciente).
# El calentamiento (NA) es min_obs = numero de coeficientes + 2: lo minimo para estimar con algun grado de libertad.
ajustar_regresion <- function(y, tipo, s = NULL, usar_log = FALSE) {
  n <- length(y)
  X <- matriz_diseno(n, tipo, s)
  p <- ncol(X)
  min_obs <- p + 2
  if (n < min_obs + 1) stop("y es demasiado corta para esta tendencia (se necesitan al menos ", min_obs + 1, " observaciones).")
  z <- if (usar_log) log(y) else y
  yhat <- rep(NA_real_, n)
  for (t in (min_obs + 1):n) {
    b <- coef_mco(X[1:(t - 1), , drop = FALSE], z[1:(t - 1)])
    pred <- sum(X[t, ] * b)
    yhat[t] <- if (usar_log) exp(pred) else pred
  }
  b_final <- coef_mco(X, z)
  names(b_final) <- colnames(X)
  list(yhat = yhat,
       pronosticar = function(h) {
         Xf <- matriz_diseno(n + h, tipo, s)[(n + 1):(n + h), , drop = FALSE]
         pred <- as.numeric(Xf %*% b_final)
         if (usar_log) exp(pred) else pred
       },
       parametros = list(tipo = tipo, s = s, coeficientes = b_final, min_obs = min_obs, escala_log = usar_log))
}

# Tendencia lineal, cuadratica o exponencial. La exponencial ajusta log(y) = a + b t y devuelve exp(a + b t)
# (sin correccion de sesgo de la transformacion).
ajustar_tendencia <- function(y, tipo = c("lineal", "cuadratica", "exponencial")) {
  tipo <- match.arg(tipo)
  validar_y(y, 5)
  if (tipo == "exponencial" && any(y <= 0)) stop("La tendencia exponencial requiere y > 0.")
  ajustar_regresion(y, tipo, usar_log = (tipo == "exponencial"))
}

# Tendencia lineal mas ficticias estacionales de suma cero.
ajustar_tendencia_estacional <- function(y, s) {
  validar_entero(s, "s", 2)
  validar_y(y, s + 4)
  ajustar_regresion(y, "estacional", s = s)
}
