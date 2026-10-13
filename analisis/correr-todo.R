# correr-todo.R — ejecutar desde la raiz del repositorio:
#   source("analisis/correr-todo.R")

for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f, encoding = "UTF-8")
source("analisis/verificaciones.R", encoding = "UTF-8")

# Metadatos de las series (R/09-config.R).
series <- definir_series()

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
candidatas <- definir_candidatas()
tabla_cand <- lapply(seq_along(series), function(i) evaluar_candidatas(datos[[i]]$valor, candidatas[[i]]))
for (i in seq_along(series)) { cat("\nSerie", i, "- candidatas\n"); print(as.data.frame(tabla_cand[[i]])) }

# Identidad: varianza de (grad Y) = 2 gamma(0) - 2 gamma(1) sobre cada serie (con log en la 1).
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
elegida <- definir_elegidas()
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

# --- Seccion 6: pronostico y evaluacion (origen movil, ganador, regresiones, validacion, pronostico final) ---
# Los resultados pesados se guardan en resultados/serieN.rds y el informe los lee (se regeneran corriendo este script).
dir.create("resultados", showWarnings = FALSE)
asignacion <- readr::read_csv("analisis/asignacion-metodos.csv", col_types = readr::cols(.default = "c"), progress = FALSE)
resultados <- list()
for (i in seq_along(series)) {
  s <- series[[i]]
  d <- datos[[i]]
  cfg <- list(s = s$s_cand[1], H = s$H, frecuencia = s$frecuencia, unidad = s$unidad,
              asignacion = as.data.frame(dplyr::filter(asignacion, serie == as.character(i))))
  res <- analizar_serie(d$valor, d$fecha, cfg)
  resultados[[i]] <- res
  cat("\n=== Serie", i, "- origen movil h = 1 (n0 =", res$n0, ") ===\n")
  print(as.data.frame(res$tabla_h1[, c("metodo", "rol", "parametros", "n", "rmse", "mad", "mape", "razon_rmse_ref")]), digits = 4)
  cat("\n=== Serie", i, "- origen movil h =", res$H, "===\n")
  print(as.data.frame(res$tabla_H[, c("metodo", "rol", "parametros", "n", "rmse", "mad", "mape", "razon_rmse_ref")]), digits = 4)
  cat("Ganador h = 1:", res$ganador_h1, "| ganador h =", res$H, ":", res$ganador_H, "\n")

  # figuras: curvas y mapas de MSE, residuos de las regresiones, validacion del ganador y pronostico final
  for (a in res$ajustes) {
    if (!is.null(a$opt)) guardar(grafico_mse(a$opt, a$etiqueta), sprintf("serie%d-mse-%s.png", i, a$id), h = if (length(a$parametros_elegidos) >= 3) 6 else 4.5)
  }
  for (id in names(res$regresiones)) {
    rg <- res$regresiones[[id]]
    guardar(correlograma(rg$r, rg$mco$n, paste0("ACF de los residuos: ", rg$etiqueta, " (serie ", i, ")")), sprintf("serie%d-residuos-%s.png", i, id))
  }
  v <- res$validacion
  g1 <- res$ajustes[[res$ganador_h1]]
  diaria <- s$frecuencia == "diaria"
  brk <- if (diaria) "1 week" else NULL
  lab <- if (diaria) "%m-%d" else "%Y-%m"
  guardar(patchwork::wrap_plots(
    grafico_errores(v$e, v$fechas_objetivo, s$unidad, paste0("Errores de origen movil a un paso: ", g1$etiqueta, " (serie ", i, ")"), brk, lab),
    correlograma(v$r, v$n, paste0("ACF de los errores (m = ", v$m, ")")), ncol = 1),
    sprintf("serie%d-validacion.png", i), h = 7)
  gH <- res$ajustes[[res$ganador_H]]
  guardar(grafico_pronostico(d$fecha, d$valor, res$pronostico, s$unidad, gH$etiqueta, ultimos = if (diaria) 40 else 48,
                             date_breaks = brk, date_labels = if (diaria) "%Y-%m-%d" else "%Y-%m"),
          sprintf("serie%d-pronostico.png", i))
  res$ajustes <- lapply(res$ajustes, function(a) { a$ajustar <- NULL; a })  # las funciones no se guardan en el .rds
  saveRDS(res, sprintf("resultados/serie%d.rds", i))
}

# --- Sesion: sessionInfo() de la maquina que corre el analisis ---
writeLines(capture.output(sessionInfo()), "sesion-info.txt")
