# 00-lectura.R — lectura y validacion de las series; ficha y particion.
# Solo define funciones: no ejecuta nada al hacer source().

# Lee un csv (fecha, valor) y valida que cumpla lo prometido en el LEEME:
# columnas correctas, sin NA, sin fechas duplicadas, orden cronologico.
leer_serie <- function(ruta) {
  stopifnot(is.character(ruta), length(ruta) == 1, file.exists(ruta))
  d <- readr::read_csv(
    ruta,
    col_types = readr::cols(fecha = readr::col_date("%Y-%m-%d"), valor = readr::col_double()),
    progress = FALSE
  )
  if (!identical(names(d), c("fecha", "valor"))) {
    stop("El archivo debe tener exactamente las columnas 'fecha' y 'valor'.")
  }
  if (anyNA(d$fecha) || anyNA(d$valor)) stop("La serie tiene faltantes.")
  if (anyDuplicated(d$fecha) > 0) stop("La serie tiene fechas duplicadas.")
  if (is.unsorted(d$fecha, strictly = TRUE)) stop("Las fechas no estan en orden cronologico.")
  d
}

# Ficha de una serie (seccion 2a): todo se calcula desde los datos.
ficha_serie <- function(d, nombre, unidad, frecuencia) {
  stopifnot(is.data.frame(d), all(c("fecha", "valor") %in% names(d)))
  tibble::tibble(
    nombre = nombre,
    unidad = unidad,
    frecuencia = frecuencia,
    inicio = min(d$fecha),
    fin = max(d$fecha),
    T = nrow(d),
    minimo = min(d$valor),
    maximo = max(d$valor),
    media = mean(d$valor)
  )
}

# Particion fija del enunciado (seccion 2c).
# n0 = floor(0.8 T); origenes h=1: n0..T-1; origenes h=H: n0..T-H;
# m = min(floor(T/4), 24) rezagos para el correlograma.
particion_serie <- function(d, H) {
  stopifnot(is.data.frame(d), is.numeric(H), length(H) == 1, H >= 4)
  T <- nrow(d)
  n0 <- floor(0.8 * T)
  tibble::tibble(
    T = T,
    n0 = n0,
    fin_estimacion = d$fecha[n0],
    inicio_evaluacion = d$fecha[n0 + 1],
    H = H,
    origenes_h1 = length(n0:(T - 1)),
    origenes_H = length(n0:(T - H)),
    m = min(floor(T / 4), 24)
  )
}
