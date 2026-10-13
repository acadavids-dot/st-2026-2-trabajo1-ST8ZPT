# correr-todo.R — ejecutar desde la raiz del repositorio:
#   source("analisis/correr-todo.R")

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, encoding = "UTF-8")

# Metadatos de las series (H de la serie 2 pendiente de decision: 7 o 9).
series <- list(
  list(id = 1, nombre = "Llegadas mensuales de turistas (Bahia Serena)", unidad = "miles de personas", frecuencia = "mensual", H = 12),
  list(id = 2, nombre = "Habitaciones ocupadas por noche (Hotel Portal de la Bahia)", unidad = "habitaciones", frecuencia = "diaria", H = 9),
  list(id = 3, nombre = "Estudiantes matriculados en cursos virtuales", unidad = "estudiantes", frecuencia = "mensual", H = 6)
)

datos <- lapply(series, function(s) leer_serie(sprintf("datos/ST8ZPT-serie-%d.csv", s$id)))

fichas <- dplyr::bind_rows(Map(function(d, s) ficha_serie(d, s$nombre, s$unidad, s$frecuencia), datos, series))
particion <- dplyr::bind_rows(Map(function(d, s) particion_serie(d, s$H), datos, series))
print(as.data.frame(fichas))
print(as.data.frame(particion))
