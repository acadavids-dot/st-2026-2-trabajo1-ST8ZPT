# 08-flujo.R — flujo completo de la seccion 6 para UNA serie (metodos, rejillas, origen movil, tablas, ganador,
# regresiones, validacion y pronostico final). Solo define funciones: no ejecuta nada al hacer source().
# Todo entra por argumentos: ninguna funcion lee objetos globales.

# Catalogo de metodos con sus rejillas (solo en 1..n0, por MSE de un paso). s = periodo estacional de la serie.
# Las constantes de suavizamiento viven en el intervalo abierto (0, 1), como en las notas.
catalogo_metodos <- function(s) {
  cat1 <- function(etiqueta, ajustar, rejilla = NULL, fijos = list()) list(etiqueta = etiqueta, ajustar = ajustar, rejilla = rejilla, fijos = fijos)
  list(
    media = cat1("Media simple", ajustar_media),
    ingenuo = cat1("Ingenuo", ajustar_ingenuo),
    ingenuo_est = cat1(paste0("Ingenuo estacional (s = ", s, ")"), ajustar_ingenuo_estacional, fijos = list(s = s)),
    mm = cat1("Media movil MM(k)", ajustar_mm, data.frame(k = 2:12)),
    ses = cat1("Suavizamiento exponencial simple", ajustar_ses, data.frame(alpha = round(seq(0.02, 0.98, by = 0.02), 2))),
    dmm = cat1("Doble media movil DMM(k)", ajustar_dmm, data.frame(k = 3:12)),
    tend_lineal = cat1("Tendencia lineal", function(y) ajustar_tendencia(y, "lineal")),
    tend_cuad = cat1("Tendencia cuadratica", function(y) ajustar_tendencia(y, "cuadratica")),
    tend_exp = cat1("Tendencia exponencial", function(y) ajustar_tendencia(y, "exponencial")),
    tend_est = cat1(paste0("Tendencia + ficticias (s = ", s, ")"), ajustar_tendencia_estacional, fijos = list(s = s)),
    holt = cat1("Holt", ajustar_holt,
                expand.grid(alpha = round(seq(0.05, 0.95, by = 0.05), 2), beta = round(seq(0.05, 0.95, by = 0.05), 2))),
    hw = cat1(paste0("Holt-Winters aditivo (s = ", s, ")"), ajustar_hw,
              expand.grid(alpha = seq(0.1, 0.9, by = 0.1), beta = seq(0.1, 0.9, by = 0.1), gamma = seq(0.1, 0.9, by = 0.1)),
              list(s = s)),
    hw_mult = cat1(paste0("Holt-Winters multiplicativo (s = ", s, ")"), ajustar_hw_mult,
                   expand.grid(alpha = seq(0.1, 0.9, by = 0.1), beta = seq(0.1, 0.9, by = 0.1), gamma = seq(0.1, 0.9, by = 0.1)),
                   list(s = s))
  )
}

# Texto con los parametros de un metodo (para las tablas).
texto_parametros <- function(params) {
  if (length(params) == 0) return("-")
  paste(paste0(names(params), " = ", vapply(params, function(v) format(round(v, 3), nsmall = 0), character(1))), collapse = ", ")
}

# Fechas de los H periodos siguientes a `ultima`.
fechas_futuras <- function(ultima, H, frecuencia) {
  stopifnot(inherits(ultima, "Date"), H >= 1)
  if (frecuencia == "mensual") seq(ultima, by = "month", length.out = H + 1)[-1]
  else if (frecuencia == "diaria") ultima + seq_len(H)
  else stop("frecuencia debe ser 'mensual' o 'diaria'.")
}

# Tabla de origen movil para un horizonte h: RMSE, MAD y MAPE por metodo, ordenada por RMSE, con el referente como una fila.
# Exige que todos los metodos tengan exactamente los mismos origenes (condicion de comparabilidad).
tabla_horizonte <- function(ajustes, h) {
  err <- lapply(ajustes, function(a) errores_horizonte(a$om, h))
  origenes <- lapply(err, function(e) e$origenes)
  if (length(unique(origenes)) != 1) stop("Los metodos no comparten los mismos origenes en h = ", h, ".")
  if (any(vapply(err, function(e) anyNA(e$e), logical(1)))) stop("Hay errores NA en h = ", h, ".")
  filas <- lapply(seq_along(ajustes), function(i) {
    a <- ajustes[[i]]
    dplyr::bind_cols(tibble::tibble(id = a$id, metodo = a$etiqueta, rol = a$rol, parametros = texto_parametros(a$parametros_elegidos)),
                     medidas(err[[i]]$e, err[[i]]$real))
  })
  tabla <- dplyr::bind_rows(filas)
  tabla <- tabla[order(tabla$rmse), ]
  rmse_ref <- tabla$rmse[tabla$rol == "referente"]
  tabla$razon_rmse_ref <- tabla$rmse / rmse_ref
  tabla$supera_referente <- tabla$rmse < rmse_ref
  tabla$h <- h
  tabla
}

# Diagnostico de una regresion de tendencia sobre toda la serie: coeficientes con EE ordinario y HAC, pruebas t robustas y,
# sobre los residuos, correlograma, Ljung-Box (chi2_m, sin restar grados de libertad), Jarque-Bera y Durbin-Watson.
diagnostico_regresion <- function(y, tipo, s, etiqueta, m) {
  X <- matriz_diseno(length(y), tipo, s)
  z <- if (tipo == "exponencial") log(y) else y
  mh <- mco_hac(z, X)
  res <- mh$residuos
  list(etiqueta = etiqueta, tipo = tipo, escala_log = (tipo == "exponencial"), mco = mh,
       r = acf_muestral(res, m), m = m,
       lb = ljung_box(res, m, paste0("los residuos de la ", tolower(etiqueta))),
       jb = jarque_bera(res, paste0("los residuos de la ", tolower(etiqueta))),
       dw = durbin_watson(res))
}

# Flujo completo de una serie.
#   y, fechas: la serie y sus fechas; cfg: list(s, H, frecuencia, unidad, asignacion), donde asignacion es un data.frame con
#   columnas id (nombre en catalogo_metodos), rol ("sostenido", "contraste" o "referente"; exactamente un referente) y evidencia.
analizar_serie <- function(y, fechas, cfg) {
  stopifnot(is.numeric(y), !anyNA(y), length(fechas) == length(y), is.list(cfg),
            all(c("s", "H", "frecuencia", "asignacion") %in% names(cfg)))
  asig <- cfg$asignacion
  stopifnot(sum(asig$rol == "referente") == 1, all(asig$rol %in% c("sostenido", "contraste", "referente")))
  T <- length(y)
  n0 <- floor(0.8 * T)
  H <- cfg$H
  catalogo <- catalogo_metodos(cfg$s)
  stopifnot(all(asig$id %in% names(catalogo)))

  # 1) parametros por MSE de un paso en 1..n0 y 2) origen movil con parametros fijos
  ajustes <- lapply(seq_len(nrow(asig)), function(i) {
    sp <- catalogo[[asig$id[i]]]
    opt <- NULL
    params <- sp$fijos
    elegidos <- list()
    if (!is.null(sp$rejilla)) {
      opt <- optimizar(y, sp$ajustar, sp$rejilla, n0, sp$fijos)
      elegidos <- as.list(opt$mejor[, names(opt$mejor) != "mse", drop = FALSE])
      params <- c(sp$fijos, elegidos)
    }
    list(id = asig$id[i], etiqueta = sp$etiqueta, rol = asig$rol[i], evidencia = asig$evidencia[i],
         parametros = params, parametros_elegidos = c(sp$fijos, elegidos), opt = opt,
         om = origen_movil(y, sp$ajustar, params, n0, H), ajustar = sp$ajustar)
  })
  names(ajustes) <- asig$id

  # 3) tablas de origen movil en h = 1 y h = H
  tabla_h1 <- tabla_horizonte(ajustes, 1)
  tabla_H <- tabla_horizonte(ajustes, H)
  n_h1 <- T - 1 - n0 + 1
  n_H <- T - H - n0 + 1
  stopifnot(all(tabla_h1$n == n_h1), all(tabla_H$n == n_H))

  # 4) ganadores: el mejor RMSE entre los metodos que no son el referente (el referente es una fila de comparacion)
  ganador_h1 <- tabla_h1$id[tabla_h1$rol != "referente"][1]
  ganador_H <- tabla_H$id[tabla_H$rol != "referente"][1]

  # 5) regresiones de tendencia (MCO con EE ordinario y HAC) sobre toda la serie
  tipos_reg <- c(tend_lineal = "lineal", tend_cuad = "cuadratica", tend_exp = "exponencial", tend_est = "estacional")
  m_reg <- min(floor(T / 4), 24)
  regresiones <- lapply(intersect(names(tipos_reg), asig$id), function(id) {
    diagnostico_regresion(y, tipos_reg[[id]], cfg$s, catalogo[[id]]$etiqueta, m_reg)
  })
  names(regresiones) <- intersect(names(tipos_reg), asig$id)

  # 6) validacion del ganador de h = 1 sobre sus errores de origen movil (no sobre los residuos del ajuste)
  e1 <- errores_horizonte(ajustes[[ganador_h1]]$om, 1)
  m_val <- max(1, floor(length(e1$e) / 4))
  validacion <- c(list(id = ganador_h1, e = e1$e, origenes = e1$origenes, fechas_objetivo = fechas[e1$origenes + 1], m = m_val),
                  validar_errores(e1$e, m_val, paste0("los errores de origen movil a un paso de ", ajustes[[ganador_h1]]$etiqueta)))

  # 7) pronostico final: ganador de h = H reajustado sobre toda la serie; banda heuristica con la sd de los errores
  #    de origen movil a cada horizonte (no es un intervalo de confianza)
  g <- ajustes[[ganador_H]]
  fc <- do.call(g$ajustar, c(list(y), g$parametros))$pronosticar(H)
  sd_h <- vapply(seq_len(H), function(h) stats::sd(errores_horizonte(g$om, h)$e), numeric(1))
  pronostico <- tibble::tibble(h = seq_len(H), fecha = fechas_futuras(fechas[T], H, cfg$frecuencia), pronostico = fc,
                               sd_error = sd_h, banda_inf = fc - 1.96 * sd_h, banda_sup = fc + 1.96 * sd_h)

  list(T = T, n0 = n0, H = H, s = cfg$s, ajustes = ajustes, tabla_h1 = tabla_h1, tabla_H = tabla_H,
       ganador_h1 = ganador_h1, ganador_H = ganador_H, regresiones = regresiones, validacion = validacion,
       pronostico = pronostico, asignacion = asig)
}

# --- Figuras de la seccion 6 ---

# Curva (1 parametro) o mapa (2 o 3 parametros, paneles por el tercero) de MSE de un paso en 1..n0.
grafico_mse <- function(opt, etiqueta) {
  stopifnot(is.list(opt), is.data.frame(opt$curva))
  curva <- as.data.frame(opt$curva)
  pars <- names(curva)[names(curva) != "mse"]
  mejor <- as.data.frame(opt$mejor)
  cap <- paste0("n = ", opt$n_comun, " errores de un paso en el tramo de estimacion; optimo (rojo): ",
                paste(pars, "=", vapply(pars, function(p) format(round(mejor[[p]], 3)), character(1)), collapse = ", "),
                "; MSE = ", format(round(mejor$mse, 4)))
  base <- ggplot2::theme_minimal(base_size = 11)
  if (length(pars) == 1) {
    p <- ggplot2::ggplot(curva, ggplot2::aes(.data[[pars[1]]], mse)) +
      ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
      ggplot2::geom_point(data = mejor, colour = "firebrick", size = 3) +
      ggplot2::labs(x = pars[1], y = "MSE de un paso")
  } else {
    p <- ggplot2::ggplot(curva, ggplot2::aes(.data[[pars[1]]], .data[[pars[2]]], fill = mse)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(name = "MSE (log)", trans = "log10") +
      ggplot2::labs(x = pars[1], y = pars[2])
    if (length(pars) == 3) {
      p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~", pars[3])), labeller = ggplot2::label_both)
      mejor_panel <- mejor
    } else {
      mejor_panel <- mejor
    }
    p <- p + ggplot2::geom_point(data = mejor_panel, colour = "firebrick", size = 2.5, shape = 4, stroke = 1.2,
                                 ggplot2::aes(.data[[pars[1]]], .data[[pars[2]]]), inherit.aes = FALSE)
  }
  p + ggplot2::labs(title = paste0("MSE de un paso: ", etiqueta), caption = cap) + base
}

# Errores de origen movil a un paso en el tiempo (fecha del dato pronosticado).
grafico_errores <- function(e, fechas_objetivo, unidad, titulo, date_breaks = NULL, date_labels = "%Y-%m") {
  stopifnot(length(e) == length(fechas_objetivo))
  d <- data.frame(fecha = fechas_objetivo, e = e)
  ggplot2::ggplot(d, ggplot2::aes(fecha, e)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_line(linewidth = 0.5) + ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_x_date(date_labels = date_labels, date_breaks = if (is.null(date_breaks)) ggplot2::waiver() else date_breaks) +
    ggplot2::labs(title = titulo, x = "Fecha del dato pronosticado", y = paste0("Error (", unidad, ")"),
                  caption = paste0("n = ", length(e), " errores de origen movil a un paso")) +
    ggplot2::theme_minimal(base_size = 11)
}

# Ultimo tramo de la serie, los H pronosticos y la banda heuristica.
grafico_pronostico <- function(fechas, y, pronostico, unidad, metodo, ultimos = 48, date_breaks = NULL, date_labels = "%Y-%m") {
  stopifnot(length(fechas) == length(y), is.data.frame(pronostico))
  idx <- seq(max(1, length(y) - ultimos + 1), length(y))
  obs <- data.frame(fecha = fechas[idx], valor = y[idx])
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = pronostico, ggplot2::aes(fecha, ymin = banda_inf, ymax = banda_sup), fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(data = obs, ggplot2::aes(fecha, valor)) +
    ggplot2::geom_line(data = pronostico, ggplot2::aes(fecha, pronostico), colour = "firebrick") +
    ggplot2::geom_point(data = pronostico, ggplot2::aes(fecha, pronostico), colour = "firebrick", size = 1.2) +
    ggplot2::scale_x_date(date_labels = date_labels, date_breaks = if (is.null(date_breaks)) ggplot2::waiver() else date_breaks) +
    ggplot2::labs(title = paste0("Pronostico a ", nrow(pronostico), " pasos: ", metodo),
                  x = "Fecha", y = unidad,
                  caption = paste0("Ultimas ", length(idx), " observaciones (n total = ", length(y), ").\n",
                                   "Banda HEURISTICA: pronostico +-1,96 sd de los errores de origen movil a cada horizonte (no es un intervalo de confianza).")) +
    ggplot2::theme_minimal(base_size = 11)
}
