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
