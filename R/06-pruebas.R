# 06-pruebas.R — plantilla de las pruebas de hipotesis con los seis elementos.
# Solo define funciones: no ejecuta nada al hacer source().

# Construye una prueba con los seis elementos que exige el enunciado, en este orden:
# 1) H0 y H1; 2) estadistico, formula y distribucion bajo H0 con grados de libertad;
# 3) region de rechazo con el valor critico calculado; 4) valor observado y valor p;
# 5) decision; 6) lectura en terminos de la serie o del metodo.
# La decision se calcula aqui (p < alpha), no se escribe a mano.
# region: texto de la region de rechazo con el critico ya calculado (lo arma cada prueba).
nueva_prueba <- function(nombre, h0, h1, estadistico, formula, distribucion, region, valor_critico,
                         observado, valor_p, lectura_rechaza, lectura_no_rechaza, alpha = 0.05) {
  stopifnot(is.numeric(observado), is.numeric(valor_p), valor_p >= 0, valor_p <= 1)
  rechaza <- valor_p < alpha
  structure(list(
    nombre = nombre, h0 = h0, h1 = h1,
    estadistico = estadistico, formula = formula, distribucion = distribucion,
    alpha = alpha, region = region, valor_critico = valor_critico,
    observado = observado, valor_p = valor_p,
    rechaza = rechaza,
    decision = if (rechaza) "Se rechaza H0" else "No se rechaza H0",
    lectura = if (rechaza) lectura_rechaza else lectura_no_rechaza
  ), class = "prueba6")
}

# Numero con coma decimal (prosa del informe).
formatear_num <- function(x, cifras = 4) format(signif(x, cifras), decimal.mark = ",", scientific = FALSE, trim = TRUE)

# Valor p legible: "< 0,0001" si es muy pequeno.
formatear_p <- function(p) if (p < 1e-4) "< 0,0001" else formatear_num(p, 3)

# Devuelve los seis elementos como lineas de markdown, para imprimirlas en el .qmd con cat(sep = "\n").
imprimir_prueba <- function(p) {
  stopifnot(inherits(p, "prueba6"))
  c(
    paste0("**", p$nombre, "**"),
    "",
    paste0("1. **Hipótesis.** $H_0$: ", p$h0, ". $H_1$: ", p$h1, "."),
    paste0("2. **Estadístico.** $", p$estadistico, " = ", p$formula, "$. Bajo $H_0$: ", p$distribucion, "."),
    paste0("3. **Región de rechazo al ", formatear_num(100 * p$alpha, 2), " %.** ", p$region, "."),
    paste0("4. **Valor observado.** $", p$estadistico, " = ", formatear_num(p$observado), "$; valor p = ", formatear_p(p$valor_p), "."),
    paste0("5. **Decisión.** ", p$decision, "."),
    paste0("6. **Lectura.** ", p$lectura)
  )
}

# t de media cero: H0: mu = 0 contra H1: mu != 0. Bajo H0 e independencia, t = ybar / (s / sqrt(n)) ~ t_{n-1}.
# Para errores de pronostico de un paso, mu != 0 significa sesgo sistematico.
prueba_t_media <- function(e, que = "los errores", alpha = 0.05) {
  stopifnot(is.numeric(e), !anyNA(e), length(e) > 2)
  n <- length(e)
  t_obs <- mean(e) / (stats::sd(e) / sqrt(n))
  crit <- stats::qt(1 - alpha / 2, n - 1)
  nueva_prueba(
    nombre = paste0("Prueba t de media cero sobre ", que),
    h0 = "$\\mu = 0$", h1 = "$\\mu \\neq 0$",
    estadistico = "t", formula = "\\frac{\\bar{x}}{s/\\sqrt{n}}",
    distribucion = paste0("$t$ de Student con ", n - 1, " grados de libertad"),
    region = paste0("$|t| > ", formatear_num(crit), "$"), valor_critico = crit,
    observado = t_obs, valor_p = 2 * stats::pt(abs(t_obs), n - 1, lower.tail = FALSE),
    lectura_rechaza = paste0("La media de ", que, " (", formatear_num(mean(e)), ") es distinta de 0: el método se equivoca en una dirección (sesgo)."),
    lectura_no_rechaza = paste0("No hay evidencia de que la media de ", que, " (", formatear_num(mean(e)), ") sea distinta de 0: no se detecta sesgo."),
    alpha = alpha
  )
}

# Jarque-Bera: H0: normalidad (asimetria 0 y curtosis 3). JB = n/6 * (S^2 + (K-3)^2/4) ~ chi2_2.
# S y K son la asimetria y la curtosis muestrales con momentos de divisor n.
jarque_bera <- function(x, que = "los errores", alpha = 0.05) {
  stopifnot(is.numeric(x), !anyNA(x), length(x) > 3)
  n <- length(x)
  xc <- x - mean(x)
  m2 <- mean(xc^2)
  S <- mean(xc^3) / m2^1.5
  K <- mean(xc^4) / m2^2
  JB <- n / 6 * (S^2 + (K - 3)^2 / 4)
  crit <- stats::qchisq(1 - alpha, 2)
  res <- nueva_prueba(
    nombre = paste0("Prueba de Jarque-Bera sobre ", que),
    h0 = "asimetría $= 0$ y curtosis $= 3$ (normalidad)", h1 = "asimetría $\\neq 0$ o curtosis $\\neq 3$",
    estadistico = "JB", formula = "\\frac{n}{6}\\left(S^2 + \\frac{(K-3)^2}{4}\\right)",
    distribucion = "$\\chi^2$ con 2 grados de libertad",
    region = paste0("$JB > ", formatear_num(crit), "$"), valor_critico = crit,
    observado = JB, valor_p = stats::pchisq(JB, 2, lower.tail = FALSE),
    lectura_rechaza = paste0(que, " no son normales ($S =$ ", formatear_num(S), ", $K =$ ", formatear_num(K), "): la banda de ±1,96 desviaciones debe leerse con cautela."),
    lectura_no_rechaza = paste0("No hay evidencia contra la normalidad de ", que, " ($S =$ ", formatear_num(S), ", $K =$ ", formatear_num(K), ")."),
    alpha = alpha
  )
  res$asimetria <- S
  res$curtosis <- K
  res
}

# Durbin-Watson como diagnostico (no como prueba con tabla): d = sum_{t>=2} (e_t - e_{t-1})^2 / sum e_t^2.
# Se lee con d ~ 2 (1 - r1): d cerca de 2 = sin autocorrelacion de orden 1; d < 2 = positiva; d > 2 = negativa.
durbin_watson <- function(e) {
  stopifnot(is.numeric(e), !anyNA(e), length(e) > 2)
  n <- length(e)
  d <- sum((e[-1] - e[-n])^2) / sum(e^2)
  r1 <- acf_muestral(e, 1)
  tibble::tibble(d = d, r1 = r1, aprox_2_1_menos_r1 = 2 * (1 - r1),
                 lectura = if (d < 1.5) "d claramente menor que 2: autocorrelación positiva de orden 1"
                           else if (d > 2.5) "d claramente mayor que 2: autocorrelación negativa de orden 1"
                           else "d cerca de 2: sin autocorrelación de orden 1 apreciable")
}
