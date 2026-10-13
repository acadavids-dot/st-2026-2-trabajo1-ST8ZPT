# 09-config.R — configuracion de las tres series, compartida por correr-todo.R y por el informe.
# Solo define funciones: no ejecuta nada al hacer source().

# Metadatos de las series.
#   H: horizonte de la pregunta de pronostico (serie 2: 9, porque la serie termina un viernes y la semana
#      lunes-domingo siguiente queda a 9 pasos: sabado y domingo son los pasos 1 y 2).
#   ks: rezagos del grafico de rezagos; s_cand: periodos estacionales a examinar (el primero es el que se usa en la seccion 6).
definir_series <- function() {
  list(
    list(id = 1, nombre = "Llegadas mensuales de turistas (Bahia Serena)", unidad = "miles de personas",
         frecuencia = "mensual", H = 12, ks = c(1, 6, 12), s_cand = 12),
    list(id = 2, nombre = "Habitaciones ocupadas por noche (Hotel Portal de la Bahia)", unidad = "habitaciones",
         frecuencia = "diaria", H = 9, ks = c(1, 7), s_cand = 7, date_breaks = "1 month"),
    list(id = 3, nombre = "Estudiantes matriculados en cursos virtuales", unidad = "estudiantes",
         frecuencia = "mensual", H = 6, ks = c(1, 6, 12), s_cand = c(12, 6))
  )
}

# Una transformacion candidata: logaritmo (log), diferencias (1-B)^d (1-B^s)^D y/o tendencia polinomial quitada (tend = grado).
cand <- function(nombre, log = FALSE, d = 0, D = 0, s = 1, tend = NULL) {
  list(nombre = nombre, log = log, d = d, D = D, s = s, tend = tend)
}

# Transformaciones candidatas por serie (seccion 4b y 4c).
definir_candidatas <- function() {
  list(
    list(cand("Y"), cand("log Y", log = TRUE), cand("(1-B) log Y", TRUE, 1), cand("(1-B^12) log Y", TRUE, 0, 1, 12),
         cand("(1-B)(1-B^12) log Y", TRUE, 1, 1, 12), cand("(1-B)(1-B^12) Y", FALSE, 1, 1, 12)),
    list(cand("Y"), cand("(1-B) Y", FALSE, 1), cand("(1-B)^2 Y", FALSE, 2), cand("(1-B^7) Y", FALSE, 0, 1, 7),
         cand("(1-B) log Y", TRUE, 1)),
    list(cand("Y"), cand("(1-B) Y", FALSE, 1), cand("(1-B)^2 Y", FALSE, 2), cand("(1-B^6) Y", FALSE, 0, 1, 6),
         cand("(1-B^12) Y", FALSE, 0, 1, 12), cand("log Y", TRUE), cand("Y - tendencia lineal", tend = 1))
  )
}

# Transformacion ELEGIDA por serie (seccion 4c): sobre ella se hace todo el analisis de la seccion 5.
#   S1: (1-B)(1-B^12) log Y; S2: (1-B) Y; S3: Y menos su tendencia lineal (MCO), porque la primera diferencia
#   sobrediferencia (r1 = -0,51) y los residuos de la recta no rechazan ruido blanco.
definir_elegidas <- function() {
  list(cand("(1-B)(1-B^12) log Y", TRUE, 1, 1, 12), cand("(1-B) Y", FALSE, 1), cand("Y - tendencia lineal", tend = 1))
}
