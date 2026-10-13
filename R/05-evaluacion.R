# 05-evaluacion.R — eleccion de parametros, origen movil, medidas y validacion de errores.
# Solo define funciones: no ejecuta nada al hacer source().

# Elige los parametros de un metodo por MSE de un paso en el tramo de estimacion 1..n0.
#   ajustar: funcion del catalogo; rejilla: data.frame con una columna por parametro (nombres = argumentos de ajustar);
#   fijos: lista de argumentos que no se optimizan (p. ej. list(s = 12)).
# Como yhat[t] solo usa el pasado, un unico ajuste sobre y[1:n0] da todos los pronosticos de un paso del tramo.
# Los calentamientos (NA) difieren entre puntos de la rejilla (k = 2 pierde 2 datos, k = 12 pierde 12): para que el MSE
# compare lo mismo, se promedia solo sobre las posiciones donde TODOS los puntos de la rejilla tienen pronostico.
optimizar <- function(y, ajustar, rejilla, n0, fijos = list()) {
  stopifnot(is.numeric(y), !anyNA(y), is.function(ajustar), is.data.frame(rejilla), nrow(rejilla) >= 1,
            is.numeric(n0), length(n0) == 1, n0 >= 3, n0 < length(y), is.list(fijos))
  yy <- y[1:n0]
  yhats <- lapply(seq_len(nrow(rejilla)), function(i) {
    do.call(ajustar, c(list(yy), fijos, as.list(rejilla[i, , drop = FALSE])))$yhat
  })
  comun <- Reduce(`&`, lapply(yhats, function(v) !is.na(v)))
  if (!any(comun)) stop("No hay posiciones con pronostico en todos los puntos de la rejilla.")
  mse <- vapply(yhats, function(v) mean((yy[comun] - v[comun])^2), numeric(1))
  curva <- tibble::as_tibble(rejilla)
  curva$mse <- mse
  list(mejor = curva[which.min(mse), ], curva = curva, n_comun = sum(comun), posiciones = which(comun))
}

# Origen movil con parametros fijos. Para cada origen o = n0..T-1 se ajusta el metodo SOLO con y[1:o] y se pronostica
# hasta min(H, T-o) pasos. Devuelve matrices (filas = origenes, columnas = horizontes h = 1..H) con NA mas alla de T:
#   errores = real - pronostico; real; pronostico; origenes.
# Los errores de horizonte h existen para los origenes o <= T-h: h = 1 usa n0..T-1; h = H usa n0..T-H.
origen_movil <- function(y, ajustar, parametros = list(), n0, H) {
  stopifnot(is.numeric(y), !anyNA(y), is.function(ajustar), is.list(parametros),
            is.numeric(n0), n0 >= 3, n0 < length(y), is.numeric(H), H >= 1)
  T <- length(y)
  origenes <- n0:(T - 1)
  vacio <- matrix(NA_real_, nrow = length(origenes), ncol = H, dimnames = list(origenes, paste0("h", seq_len(H))))
  real <- pron <- vacio
  for (i in seq_along(origenes)) {
    o <- origenes[i]
    hmax <- min(H, T - o)
    fit <- do.call(ajustar, c(list(y[1:o]), parametros))
    real[i, seq_len(hmax)] <- y[(o + 1):(o + hmax)]
    pron[i, seq_len(hmax)] <- fit$pronosticar(hmax)
  }
  list(errores = real - pron, real = real, pronostico = pron, origenes = origenes, n0 = n0, H = H)
}

# Errores (y valores reales) del horizonte h de un resultado de origen_movil, sin los NA de los origenes que no alcanzan.
errores_horizonte <- function(om, h) {
  stopifnot(is.list(om), h >= 1, h <= om$H)
  ok <- !is.na(om$errores[, h])
  list(e = om$errores[ok, h], real = om$real[ok, h], origenes = om$origenes[ok])
}

# RMSE, MAD (promedio de |e|) y MAPE (promedio de |e/real| en %).
medidas <- function(e, real) {
  stopifnot(is.numeric(e), is.numeric(real), length(e) == length(real), length(e) >= 1, !anyNA(e), !anyNA(real))
  if (any(real == 0)) warning("MAPE indefinido: hay valores reales iguales a 0.")
  tibble::tibble(n = length(e), rmse = sqrt(mean(e^2)), mad = mean(abs(e)), mape = 100 * mean(abs(e / real)))
}

# Validacion de los errores de origen movil de un paso (NO de los residuos del ajuste):
# correlograma con la banda calculada sobre n = numero de errores, t de media cero, Ljung-Box (chi2_m),
# Jarque-Bera, y Durbin-Watson como diagnostico.
validar_errores <- function(e, m, que = "los errores de origen movil a un paso") {
  stopifnot(is.numeric(e), !anyNA(e), length(e) > m + 1)
  n <- length(e)
  r <- acf_muestral(e, m)
  list(n = n, r = r, banda_rb = 1.96 / sqrt(n),
       media = prueba_t_media(e, que),
       lb = ljung_box(e, m, que),
       jb = jarque_bera(e, que),
       dw = durbin_watson(e))
}
