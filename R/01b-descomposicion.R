# 01b-descomposicion.R — descomposicion clasica a mano y evidencia de estacionalidad.
# Solo define funciones: no ejecuta nada al hacer source().

# Media movil centrada de periodo s, sin stats::filter().
# s impar: promedio simple de s valores centrados en t.
# s par: media movil 2 x s, con pesos 1/(2s) en los extremos y 1/s en el medio,
# para que la ventana quede centrada en t (una ventana de tamano par no lo esta).
# Los s/2 primeros y ultimos valores quedan NA (no hay ventana completa).
media_movil_centrada <- function(y, s) {
  stopifnot(is.numeric(y), !anyNA(y), is.numeric(s), length(s) == 1, s >= 2, length(y) > s)
  n <- length(y)
  tendencia <- rep(NA_real_, n)
  if (s %% 2 == 1) {
    a <- (s - 1) / 2
    w <- rep(1 / s, s)
  } else {
    a <- s / 2
    w <- c(0.5, rep(1, s - 1), 0.5) / s
  }
  for (t in (a + 1):(n - a)) tendencia[t] <- sum(w * y[(t - a):(t + a)])
  tendencia
}

# Descomposicion clasica. tipo: "aditiva" (y = T + S + I) o "multiplicativa" (y = T * S * I).
# 1) tendencia = media movil centrada; 2) se quita la tendencia (resta o cociente);
# 3) indice crudo de cada estacion = promedio de lo que queda en esa posicion del ciclo;
# 4) se normalizan los indices (suma 0 en aditiva, media 1 en multiplicativa) para que
#    la estacionalidad no se lleve parte del nivel; 5) el residuo es lo que resta.
# La posicion 1 es la primera observacion de y (no necesariamente enero).
descomposicion_clasica <- function(y, s, tipo = c("aditiva", "multiplicativa")) {
  tipo <- match.arg(tipo)
  stopifnot(is.numeric(y), !anyNA(y), is.numeric(s), length(s) == 1, s >= 2, length(y) >= 2 * s)
  if (tipo == "multiplicativa" && any(y <= 0)) stop("La descomposicion multiplicativa requiere y > 0.")
  n <- length(y)
  tendencia <- media_movil_centrada(y, s)
  sin_tendencia <- if (tipo == "aditiva") y - tendencia else y / tendencia
  pos <- ((seq_len(n) - 1) %% s) + 1
  crudo <- vapply(seq_len(s), function(j) mean(sin_tendencia[pos == j], na.rm = TRUE), numeric(1))
  indices <- if (tipo == "aditiva") crudo - mean(crudo) else crudo / mean(crudo)
  estacional <- indices[pos]
  residuo <- if (tipo == "aditiva") y - tendencia - estacional else y / (tendencia * estacional)
  list(tipo = tipo, s = s, tendencia = tendencia, indices = indices,
       estacional = estacional, residuo = residuo, sin_tendencia = sin_tendencia)
}

# Evidencia numerica para (i) estacional o no y (ii) aditiva o multiplicativa.
# - prop_var_estacional: parte de la varianza de la serie sin tendencia que explican los indices
#   (cerca de 0 = no hay patron que se repita cada s periodos).
# - r_s_residuo: autocorrelacion en el rezago s del residuo de la descomposicion aditiva.
# - cor_amplitud_nivel: correlacion, por ciclo completo, entre el nivel (media) y la amplitud (sd);
#   cerca de 1 = la amplitud crece con el nivel = senal de modelo multiplicativo.
# - cv_residuo_*: desviacion relativa del residuo en cada version (menor = mejor ajuste relativo).
evidencia_estacional <- function(y, s) {
  stopifnot(is.numeric(y), !anyNA(y), length(y) >= 2 * s)
  da <- descomposicion_clasica(y, s, "aditiva")
  sa <- da$sin_tendencia[!is.na(da$sin_tendencia)]
  n_ciclos <- floor(length(y) / s)
  m <- matrix(y[seq_len(n_ciclos * s)], ncol = s, byrow = TRUE)
  res <- da$residuo[!is.na(da$residuo)]
  out <- tibble::tibble(
    s = s,
    prop_var_estacional = stats::var(da$estacional[!is.na(da$tendencia)]) / stats::var(sa),
    r_s_residuo = acf_muestral(res, s)[s],
    cor_amplitud_nivel = stats::cor(rowMeans(m), apply(m, 1, stats::sd)),
    cv_residuo_aditiva = stats::sd(res / da$tendencia[!is.na(da$residuo)])
  )
  if (all(y > 0)) {
    dm <- descomposicion_clasica(y, s, "multiplicativa")
    out$cv_residuo_multiplicativa <- stats::sd(dm$residuo[!is.na(dm$residuo)])
  }
  out
}
