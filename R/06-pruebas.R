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
    paste0("1. **Hipotesis.** H0: ", p$h0, ". H1: ", p$h1, "."),
    paste0("2. **Estadistico.** ", p$estadistico, " = ", p$formula, ". Bajo H0: ", p$distribucion, "."),
    paste0("3. **Region de rechazo al ", formatear_num(100 * p$alpha, 2), " %.** ", p$region, "."),
    paste0("4. **Valor observado.** ", p$estadistico, " = ", formatear_num(p$observado), "; valor p = ", formatear_p(p$valor_p), "."),
    paste0("5. **Decision.** ", p$decision, "."),
    paste0("6. **Lectura.** ", p$lectura)
  )
}
