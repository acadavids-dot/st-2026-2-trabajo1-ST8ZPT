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
#   estacional: 1, t y las s-1 ficticias de suma cero de las notas (clase 5): Z_j = D_j - D_1 para j = 2..s,
#   es decir 1 en la estacion j, -1 en la estacion 1 y 0 en las demas. Los coeficientes son directamente gamma_2..gamma_s
#   y gamma_1 = -(gamma_2 + ... + gamma_s), asi que los efectos estacionales suman 0.
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
      D <- vapply(2:s, function(j) ifelse(pos == j, 1, ifelse(pos == 1, -1, 0)), numeric(n))
      colnames(D) <- paste0("D", 2:s, "-D1")
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
  # Correccion de retransformacion de Duan (notas clase 4, prop. 4): exp(a + b t) estima la MEDIANA condicional, y el MSE
  # se minimiza en la media. El factor "smearing" (1/n) sum exp(e_t) de los residuos en log corrige sin exigir normalidad.
  duan <- function(Xa, za, b) if (usar_log) mean(exp(za - as.numeric(Xa %*% b))) else 1
  yhat <- rep(NA_real_, n)
  for (t in (min_obs + 1):n) {
    Xa <- X[1:(t - 1), , drop = FALSE]
    b <- coef_mco(Xa, z[1:(t - 1)])
    pred <- sum(X[t, ] * b)
    yhat[t] <- if (usar_log) exp(pred) * duan(Xa, z[1:(t - 1)], b) else pred
  }
  b_final <- coef_mco(X, z)
  names(b_final) <- colnames(X)
  factor_duan <- duan(X, z, b_final)
  list(yhat = yhat,
       pronosticar = function(h) {
         Xf <- matriz_diseno(n + h, tipo, s)[(n + 1):(n + h), , drop = FALSE]
         pred <- as.numeric(Xf %*% b_final)
         if (usar_log) exp(pred) * factor_duan else pred
       },
       parametros = list(tipo = tipo, s = s, coeficientes = b_final, min_obs = min_obs, escala_log = usar_log,
                         factor_duan = factor_duan))
}

# Tendencia lineal, cuadratica o exponencial. La exponencial ajusta log(y) = a + b t y devuelve
# exp(a + b t) * factor de Duan (ver ajustar_regresion).
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

# Validacion de una constante de suavizamiento en el intervalo abierto (0, 1), como en las notas.
validar_constante <- function(x, nombre) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0 || x >= 1) stop(nombre, " debe estar en el intervalo abierto (0, 1).")
  invisible(TRUE)
}

# Suavizamiento exponencial simple (notas clase 4): Yhat_{t+1} = alpha Y_t + (1 - alpha) Yhat_t.
# Inicializacion: L_1 = Y_1, es decir Yhat_2 = Y_1 (con alpha -> 1 el metodo tiende al ingenuo).
# Esa convencion es la que da 11 errores con 12 datos en el ejemplo del bono de las notas (la 1.a posicion queda en NA).
ajustar_ses <- function(y, alpha) {
  validar_constante(alpha, "alpha")
  validar_y(y, 3)
  n <- length(y)
  yhat <- rep(NA_real_, n)
  yhat[2] <- y[1]
  for (t in 2:(n - 1)) yhat[t + 1] <- alpha * y[t] + (1 - alpha) * yhat[t]
  nivel <- alpha * y[n] + (1 - alpha) * yhat[n]  # Yhat_{T+1}
  list(yhat = yhat,
       pronosticar = function(h) rep(nivel, h),
       parametros = list(alpha = alpha, nivel_inicial = y[1], nivel = nivel))
}

# Dobles medias moviles (notas clase 4, def. 3). MM_t(k) = promedio de y[(t-k+1):t];
# DMM_t(k) = promedio de MM_{t-i}(k), i = 0..k-1, definido para t >= 2k-1.
#   E_t = 2 MM_t - DMM_t (nivel);  b_t = 2/(k-1) (MM_t - DMM_t) (pendiente).
#   Pronostico muestral Yhat_{t+1} = E_t + b_t; extramuestral Yhat_{T+L} = E_T + b_T L.
# Hacen falta 2k-1 observaciones antes del primer pronostico, asi que yhat[t] existe desde t = 2k.
ajustar_dmm <- function(y, k) {
  validar_entero(k, "k", 2)
  validar_y(y, 2 * k)
  n <- length(y)
  acum <- c(0, cumsum(y))
  mm <- rep(NA_real_, n)
  for (t in k:n) mm[t] <- (acum[t + 1] - acum[t - k + 1]) / k
  dmm <- rep(NA_real_, n)
  for (t in (2 * k - 1):n) dmm[t] <- mean(mm[(t - k + 1):t])
  nivel_t <- 2 * mm - dmm
  pend_t <- 2 / (k - 1) * (mm - dmm)
  yhat <- c(NA_real_, nivel_t[-n] + pend_t[-n])  # yhat[t+1] = E_t + b_t (NA mientras DMM no exista)
  nivel <- nivel_t[n]
  pendiente <- pend_t[n]
  list(yhat = yhat,
       pronosticar = function(h) nivel + pendiente * seq_len(h),
       parametros = list(k = k, nivel = nivel, pendiente = pendiente))
}

# Holt (notas clase 4, def. 4): L_t = alpha Y_t + (1 - alpha) (L_{t-1} + T_{t-1}),
# T_t = beta (L_t - L_{t-1}) + (1 - beta) T_{t-1}, Yhat_{t+1} = L_t + T_t.
# Inicializacion de las notas: L_1 = Y_1 y T_1 = 0 (yhat[1] = NA, yhat[2] = Y_1).
ajustar_holt <- function(y, alpha, beta) {
  validar_constante(alpha, "alpha")
  validar_constante(beta, "beta")
  validar_y(y, 3)
  n <- length(y)
  L <- numeric(n)
  Tt <- numeric(n)
  L[1] <- y[1]
  Tt[1] <- 0
  yhat <- rep(NA_real_, n)
  for (t in 2:n) {
    yhat[t] <- L[t - 1] + Tt[t - 1]
    L[t] <- alpha * y[t] + (1 - alpha) * yhat[t]
    Tt[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * Tt[t - 1]
  }
  nivel <- L[n]
  pendiente <- Tt[n]
  list(yhat = yhat,
       pronosticar = function(h) nivel + pendiente * seq_len(h),
       parametros = list(alpha = alpha, beta = beta, L1 = y[1], T1 = 0, nivel = nivel, pendiente = pendiente,
                         nivel_t = L, pendiente_t = Tt))
}

# Estados iniciales de Holt-Winters (notas clase 4, diap. 68). Con m ciclos completos y promedios por ciclo ybar_1..ybar_m:
#   b0 = (ybar_m - ybar_1) / ((m - 1) s);  L0 = ybar_1 - b0 (s - 1)/2   (nivel en t = 1);
#   S0_j = promedio sobre los ciclos de [Y_{j+is} - (L0 + b0 (j + is - 1))], normalizado a suma cero.
# Multiplicativo: las notas solo dicen que los indices "promedian uno"; aqui se usa el analogo directo de la formula
# aditiva con cocientes: S0_j = promedio de Y_{j+is} / (L0 + b0 (j + is - 1)), dividido por su media (media 1).
# Se usan TODOS los ciclos completos de la serie que se ajusta (convencion de las notas, que reproduce su ejemplo).
inicializar_hw <- function(y, s, multiplicativo = FALSE) {
  m <- floor(length(y) / s)
  if (m < 2) stop("Holt-Winters necesita al menos 2 ciclos completos para inicializar.")
  yc <- y[1:(m * s)]
  ybar <- vapply(seq_len(m), function(i) mean(yc[((i - 1) * s + 1):(i * s)]), numeric(1))
  b0 <- (ybar[m] - ybar[1]) / ((m - 1) * s)
  L0 <- ybar[1] - b0 * (s - 1) / 2
  S0 <- vapply(seq_len(s), function(j) {
    t <- j + (0:(m - 1)) * s
    nivel_t <- L0 + b0 * (t - 1)
    if (multiplicativo) mean(yc[t] / nivel_t) else mean(yc[t] - nivel_t)
  }, numeric(1))
  S0 <- if (multiplicativo) S0 / mean(S0) else S0 - mean(S0)
  list(L0 = L0, b0 = b0, S0 = S0, m = m)
}

# Nucleo comun de Holt-Winters aditivo y multiplicativo. La recursion arranca en t = 1 con (L0, b0) como estado previo
# y con S0 como indices del ciclo "0" (convencion de las notas). Como la inicializacion usa los ciclos completos de y,
# yhat[t] de las primeras posiciones depende de toda la serie ajustada: es lo que hacen las notas. La evaluacion de origen
# movil NO se contamina, porque en cada origen el metodo se ajusta solo con y[1:o].
ajustar_hw_nucleo <- function(y, s, alpha, beta, gamma, multiplicativo) {
  validar_constante(alpha, "alpha")
  validar_constante(beta, "beta")
  validar_constante(gamma, "gamma")
  validar_entero(s, "s", 2)
  validar_y(y, 2 * s)
  if (multiplicativo && any(y <= 0)) stop("Holt-Winters multiplicativo requiere y > 0.")
  n <- length(y)
  ini <- inicializar_hw(y, s, multiplicativo)
  L <- ini$L0
  b <- ini$b0
  S <- ini$S0  # S[j] = ultimo indice disponible de la estacion j
  yhat <- numeric(n)
  for (t in seq_len(n)) {
    j <- ((t - 1) %% s) + 1
    P <- L + b  # nivel proyectado desde t-1
    yhat[t] <- if (multiplicativo) P * S[j] else P + S[j]
    L_nuevo <- if (multiplicativo) alpha * (y[t] / S[j]) + (1 - alpha) * P else alpha * (y[t] - S[j]) + (1 - alpha) * P
    b <- beta * (L_nuevo - L) + (1 - beta) * b
    S[j] <- if (multiplicativo) gamma * (y[t] / L_nuevo) + (1 - gamma) * S[j] else gamma * (y[t] - L_nuevo) + (1 - gamma) * S[j]
    L <- L_nuevo
  }
  nivel <- L
  pendiente <- b
  indices <- S
  list(yhat = yhat,
       pronosticar = function(h) {
         pos <- ((n + seq_len(h) - 1) %% s) + 1  # estacion de T+m; usa el indice del ultimo ciclo disponible
         if (multiplicativo) (nivel + pendiente * seq_len(h)) * indices[pos] else nivel + pendiente * seq_len(h) + indices[pos]
       },
       parametros = list(s = s, alpha = alpha, beta = beta, gamma = gamma, multiplicativo = multiplicativo,
                         L0 = ini$L0, b0 = ini$b0, S0 = ini$S0, ciclos_init = ini$m,
                         nivel = nivel, pendiente = pendiente, indices = indices))
}

# Holt-Winters aditivo (notas clase 4, diap. 67 y clase 5, diap. 49).
ajustar_hw <- function(y, s, alpha, beta, gamma) ajustar_hw_nucleo(y, s, alpha, beta, gamma, multiplicativo = FALSE)

# Holt-Winters multiplicativo: L_t = alpha Y_t/S_{t-s} + (1 - alpha)(L_{t-1} + b_{t-1}), S_t = gamma Y_t/L_t + (1 - gamma) S_{t-s},
# Yhat_{t+m} = (L_t + b_t m) S_{t+m-s ceiling(m/s)}.
ajustar_hw_mult <- function(y, s, alpha, beta, gamma) ajustar_hw_nucleo(y, s, alpha, beta, gamma, multiplicativo = TRUE)
