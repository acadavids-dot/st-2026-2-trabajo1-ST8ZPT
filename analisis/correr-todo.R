# correr-todo.R — ejecutar desde la raiz del repositorio:
#   source("analisis/correr-todo.R")

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, encoding = "UTF-8")
source("analisis/verificaciones.R", encoding = "UTF-8")

# Metadatos de las series. H de la serie 2 pendiente de decision (7 o 9).
# ks: rezagos del grafico de rezagos; s_cand: periodos estacionales a examinar.
series <- list(
  list(id = 1, nombre = "Llegadas mensuales de turistas (Bahia Serena)", unidad = "miles de personas",
       frecuencia = "mensual", H = 12, ks = c(1, 6, 12), s_cand = 12),
  list(id = 2, nombre = "Habitaciones ocupadas por noche (Hotel Portal de la Bahia)", unidad = "habitaciones",
       frecuencia = "diaria", H = 9, ks = c(1, 7), s_cand = 7, date_breaks = "1 month"),
  list(id = 3, nombre = "Estudiantes matriculados en cursos virtuales", unidad = "estudiantes",
       frecuencia = "mensual", H = 6, ks = c(1, 6, 12), s_cand = c(12, 6))
)

datos <- lapply(series, function(s) leer_serie(sprintf("datos/ST8ZPT-serie-%d.csv", s$id)))

# --- Seccion 2: fichas y particion ---
fichas <- dplyr::bind_rows(Map(function(d, s) ficha_serie(d, s$nombre, s$unidad, s$frecuencia), datos, series))
particion <- dplyr::bind_rows(Map(function(d, s) particion_serie(d, s$H), datos, series))
print(as.data.frame(fichas))
print(as.data.frame(particion))

# --- Seccion 3: descripcion, estacionalidad y correlograma (serie sin transformar) ---
dir.create("figs", showWarnings = FALSE)
guardar <- function(p, nombre, w = 8, h = 4.5) ggplot2::ggsave(file.path("figs", nombre), p, width = w, height = h, dpi = 150)

evidencia <- list()
for (i in seq_along(series)) {
  s <- series[[i]]
  d <- datos[[i]]
  y <- d$valor
  m <- particion$m[i]
  guardar(graficar_serie(d, s$nombre, s$unidad, s$date_breaks), sprintf("serie%d-tiempo.png", i))
  guardar(grafico_rezagos(y, s$ks, s$unidad), sprintf("serie%d-rezagos.png", i), h = 4)
  guardar(correlograma(acf_muestral(y, m), length(y), paste0("ACF de la serie ", i, " sin transformar (m = ", m, ")")),
          sprintf("serie%d-acf.png", i))
  for (sc in s$s_cand) {
    guardar(grafico_estacional(d, sc, s$unidad), sprintf("serie%d-estacional-s%d.png", i, sc))
    evidencia[[length(evidencia) + 1]] <- dplyr::mutate(evidencia_estacional(y, sc), serie = i, .before = 1)
  }
}
evidencia <- dplyr::bind_rows(evidencia)
print(as.data.frame(evidencia))

# Descomposicion clasica donde el diagnostico sostiene estacionalidad (se decide con la evidencia de arriba).
guardar(grafico_descomposicion(datos[[1]], descomposicion_clasica(datos[[1]]$valor, 12, "multiplicativa"), series[[1]]$unidad),
        "serie1-descomposicion.png", h = 7)

# --- Seccion 4: operador de rezago, candidatas y raices ---
cand <- function(nombre, log = FALSE, d = 0, D = 0, s = 1, tend = NULL) list(nombre = nombre, log = log, d = d, D = D, s = s, tend = tend)
candidatas <- list(
  list(cand("Y"), cand("log Y", log = TRUE), cand("(1-B) log Y", TRUE, 1), cand("(1-B^12) log Y", TRUE, 0, 1, 12),
       cand("(1-B)(1-B^12) log Y", TRUE, 1, 1, 12), cand("(1-B)(1-B^12) Y", FALSE, 1, 1, 12)),
  list(cand("Y"), cand("(1-B) Y", FALSE, 1), cand("(1-B)^2 Y", FALSE, 2), cand("(1-B^7) Y", FALSE, 0, 1, 7),
       cand("(1-B) log Y", TRUE, 1)),
  list(cand("Y"), cand("(1-B) Y", FALSE, 1), cand("(1-B)^2 Y", FALSE, 2), cand("(1-B^6) Y", FALSE, 0, 1, 6),
       cand("(1-B^12) Y", FALSE, 0, 1, 12), cand("log Y", TRUE), cand("Y - tendencia lineal", tend = 1))
)
tabla_cand <- lapply(seq_along(series), function(i) evaluar_candidatas(datos[[i]]$valor, candidatas[[i]]))
for (i in seq_along(series)) { cat("\nSerie", i, "- candidatas\n"); print(as.data.frame(tabla_cand[[i]])) }

# Identidad Var(grad Y) = 2 gamma(0) - 2 gamma(1) sobre cada serie (con log en la 1).
identidad <- dplyr::bind_rows(lapply(seq_along(series), function(i) {
  dplyr::mutate(identidad_var_diferencia(datos[[i]]$valor), serie = i, .before = 1)
}))
print(as.data.frame(identidad))

# Raices de los polinomios de diferencia de cada serie.
for (i in seq_along(series)) {
  ps <- unique(lapply(candidatas[[i]], function(cd) c(cd$d, cd$D, cd$s)))
  ps <- ps[vapply(ps, function(p) p[1] + p[2] > 0, logical(1))]
  nombres <- vapply(ps, function(p) {
    paste0(if (p[1] > 0) paste0("(1-B)", if (p[1] > 1) paste0("^", p[1]) else "") else "",
           if (p[2] > 0) paste0("(1-B^", p[3], ")") else "")
  }, character(1))
  rs <- setNames(lapply(ps, function(p) polinomio_rezago(coef_diferencia(p[1], p[2], p[3]))), nombres)
  guardar(grafico_raices(rs), sprintf("serie%d-raices.png", i), h = 5)
}

# --- Seccion 5: dependencia sobre la serie transformada ELEGIDA ---
# Elegidas: S1 = (1-B)(1-B^12) log Y; S2 = (1-B) Y; S3 = Y - tendencia lineal (MCO).
elegida <- list(
  cand("(1-B)(1-B^12) log Y", TRUE, 1, 1, 12),
  cand("(1-B) Y", FALSE, 1),
  cand("Y - tendencia lineal", tend = 1)
)
m5 <- 24  # m para Ljung-Box y correlogramas de la serie transformada (ver justificacion en el informe)
dep <- list()
for (i in seq_along(series)) {
  z <- aplicar_candidata(datos[[i]]$valor, elegida[[i]])
  n <- length(z)
  r <- acf_muestral(z, m5)
  phi <- pacf_dl(r)
  bb <- banda_bartlett(r, n)
  guardar(correlograma(r, n, paste0("ACF de la serie ", i, " transformada: ", elegida[[i]]$nombre), bartlett = bb),
          sprintf("serie%d-acf-transformada.png", i))
  guardar(correlograma(phi, n, paste0("PACF de la serie ", i, " transformada: ", elegida[[i]]$nombre), "Autocorrelacion parcial"),
          sprintf("serie%d-pacf-transformada.png", i))
  dep[[i]] <- list(z = z, n = n, r = r, phi = phi, bartlett = bb,
                   lb = ljung_box(z, m5, paste0("la serie ", i, " transformada")),
                   bloques = bloques_media_var(z, 4), ic = ic_media(z))
  cat("\n=== Serie", i, "transformada:", elegida[[i]]$nombre, "(n =", n, ") ===\n")
  cat(imprimir_prueba(dep[[i]]$lb), sep = "\n")
  print(as.data.frame(dep[[i]]$bloques))
  print(as.data.frame(dep[[i]]$ic))
  cat("Rezagos de la ACF fuera de la banda de ruido blanco:", which(abs(r) > 1.96 / sqrt(n)), "\n")
  cat("Rezagos de la ACF fuera de la banda de Bartlett:", which(abs(r) > bb), "\n")
  cat("Rezagos de la PACF fuera de la banda de ruido blanco:", which(abs(phi) > 1.96 / sqrt(n)), "\n")
}
