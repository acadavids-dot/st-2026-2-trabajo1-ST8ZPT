# 01-graficos.R — graficos descriptivos (ggplot2 + patchwork).
# Solo define funciones: no ejecuta nada al hacer source().
# Todos llevan titulo, ejes con unidad y caption con el numero de observaciones.

# Serie en el tiempo. d: data.frame con fecha (Date) y valor.
# date_breaks: separacion de las marcas del eje (p. ej. "1 month", "1 year"); NULL = automatico.
graficar_serie <- function(d, titulo, unidad, date_breaks = NULL) {
  stopifnot(is.data.frame(d), inherits(d$fecha, "Date"), is.numeric(d$valor))
  ggplot2::ggplot(d, ggplot2::aes(fecha, valor)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::scale_x_date(date_labels = "%Y-%m", date_breaks = if (is.null(date_breaks)) ggplot2::waiver() else date_breaks) +
    ggplot2::labs(title = titulo, x = "Fecha", y = unidad, caption = paste0("n = ", nrow(d), " observaciones")) +
    ggplot2::theme_minimal(base_size = 11)
}

# Grafico de rezagos: y_t (eje y) contra y_{t-k} (eje x), un panel por k, con r_k en el titulo del panel.
grafico_rezagos <- function(y, ks, unidad) {
  stopifnot(is.numeric(y), !anyNA(y), all(ks >= 1), all(ks < length(y)))
  n <- length(y)
  r <- acf_muestral(y, max(ks))
  d <- do.call(rbind, lapply(ks, function(k) {
    data.frame(x = y[1:(n - k)], y = y[(1 + k):n],
               panel = sprintf("k = %d  (r = %.2f)", k, r[k]))
  }))
  d$panel <- factor(d$panel, levels = unique(d$panel))
  ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_point(size = 1, alpha = 0.7) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::labs(title = "Grafico de rezagos", x = paste0("Y(t-k) (", unidad, ")"),
                  y = paste0("Y(t) (", unidad, ")"), caption = paste0("n = ", n, " observaciones")) +
    ggplot2::theme_minimal(base_size = 11)
}

# Grafico estacional: un trazo por ciclo sobre el indice de estacion (posicion dentro del ciclo).
# La posicion 1 es la primera observacion de la serie.
grafico_estacional <- function(d, s, unidad, etiquetas = NULL) {
  stopifnot(is.data.frame(d), is.numeric(s), s >= 2)
  n <- nrow(d)
  g <- data.frame(pos = ((seq_len(n) - 1) %% s) + 1, ciclo = ((seq_len(n) - 1) %/% s) + 1, valor = d$valor)
  p <- ggplot2::ggplot(g, ggplot2::aes(pos, valor, group = ciclo, colour = ciclo)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::scale_colour_viridis_c(name = "Ciclo") +
    ggplot2::labs(title = paste0("Grafico estacional (s = ", s, ")"), x = "Posicion en el ciclo",
                  y = unidad, caption = paste0("n = ", n, " observaciones; ", ceiling(n / s), " ciclos")) +
    ggplot2::theme_minimal(base_size = 11)
  if (is.null(etiquetas)) {
    p + ggplot2::scale_x_continuous(breaks = seq_len(s))
  } else {
    p + ggplot2::scale_x_continuous(breaks = seq_len(s), labels = etiquetas)
  }
}

# Correlograma: barras de r_1..r_m con la banda de ruido blanco +-1.96/sqrt(n).
# n es el tamano de la serie sobre la que se calculo r (no m).
correlograma <- function(r, n, titulo = "Correlograma (ACF)", ylab = "Autocorrelacion", bartlett = NULL) {
  stopifnot(is.numeric(r), is.numeric(n), length(n) == 1, n > 1)
  d <- data.frame(k = seq_along(r), r = r)
  b <- 1.96 / sqrt(n)
  p <- ggplot2::ggplot(d, ggplot2::aes(k, r)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(yintercept = c(-b, b), linetype = "dashed", colour = "steelblue") +
    ggplot2::geom_segment(ggplot2::aes(xend = k, yend = 0)) +
    ggplot2::geom_point(size = 1.2)
  cap <- paste0("n = ", n, " observaciones; banda de ruido blanco +-1,96/sqrt(n) = +-", round(b, 3), " (azul)")
  if (!is.null(bartlett)) {
    stopifnot(length(bartlett) == length(r))
    db <- data.frame(k = rep(d$k, 2), b = c(bartlett, -bartlett), lado = rep(c("sup", "inf"), each = length(r)))
    p <- p + ggplot2::geom_step(data = db, ggplot2::aes(k, b, group = lado), colour = "firebrick",
                                linetype = "dotted", inherit.aes = FALSE)
    cap <- paste0(cap, "; banda de Bartlett (rojo)")
  }
  p + ggplot2::labs(title = titulo, x = "Rezago k", y = ylab, caption = cap) +
    ggplot2::theme_minimal(base_size = 11)
}

# Raices de polinomios en B en el plano complejo, con el circulo unitario.
# raices: lista con nombre de tibbles (re, im, modulo) como los de polinomio_rezago().
grafico_raices <- function(raices) {
  stopifnot(is.list(raices), !is.null(names(raices)))
  d <- do.call(rbind, lapply(names(raices), function(nm) {
    r <- raices[[nm]]
    data.frame(polinomio = nm, re = r$re, im = r$im)
  }))
  d$polinomio <- factor(d$polinomio, levels = names(raices))
  th <- seq(0, 2 * pi, length.out = 200)
  circ <- data.frame(re = cos(th), im = sin(th))
  ggplot2::ggplot(d, ggplot2::aes(re, im)) +
    ggplot2::geom_path(data = circ, colour = "grey50") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_point(colour = "firebrick", size = 2, alpha = 0.8) +
    ggplot2::facet_wrap(~polinomio) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Raices de los polinomios en B", x = "Parte real", y = "Parte imaginaria",
                  caption = "Circulo unitario en gris; todas las raices de un factor de diferencia tienen modulo 1") +
    ggplot2::theme_minimal(base_size = 11)
}

# Panel de la descomposicion: observada, tendencia, componente estacional y residuo.
grafico_descomposicion <- function(d, desc, unidad) {
  stopifnot(is.data.frame(d), is.list(desc), length(desc$tendencia) == nrow(d))
  x <- d$fecha
  panel <- function(v, ylab, ref = NULL) {
    p <- ggplot2::ggplot(data.frame(x = x, v = v), ggplot2::aes(x, v)) +
      ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE) +
      ggplot2::labs(x = NULL, y = ylab) +
      ggplot2::theme_minimal(base_size = 10)
    if (!is.null(ref)) p <- p + ggplot2::geom_hline(yintercept = ref, linetype = "dashed")
    p
  }
  ref <- if (desc$tipo == "aditiva") 0 else 1
  patchwork::wrap_plots(
    panel(d$valor, "Observada"),
    panel(desc$tendencia, "Tendencia"),
    panel(desc$estacional, "Estacional", ref),
    panel(desc$residuo, "Residuo", ref),
    ncol = 1
  ) +
    patchwork::plot_annotation(
      title = paste0("Descomposicion clasica ", desc$tipo, " (s = ", desc$s, ")"),
      subtitle = paste0("Observada y tendencia en ", unidad, "; tendencia = media movil centrada"),
      caption = paste0("n = ", nrow(d), " observaciones"))
}
