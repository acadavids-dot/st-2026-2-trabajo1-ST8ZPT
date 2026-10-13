# 10-informe.R — formato de numeros y tablas para el informe (coma decimal en la prosa).
# Solo define funciones: no ejecuta nada al hacer source().

# Numero con `d` decimales y coma decimal.
fnum <- function(x, d = 2) format(round(x, d), nsmall = d, decimal.mark = ",", big.mark = "", trim = TRUE, scientific = FALSE)

# Porcentaje (x en proporcion) con `d` decimales.
fpct <- function(x, d = 1) paste0(fnum(100 * x, d), " %")

# Valor p legible para tablas.
fp <- function(p) ifelse(p < 1e-4, "< 0,0001", fnum(p, 3))

# Tabla html/markdown con coma decimal.
tabla <- function(df, digits = 2, ...) {
  knitr::kable(df, digits = digits, format.args = list(decimal.mark = ","), row.names = FALSE, ...)
}

# Pesos colombianos: "450.000 COP" y millones "19,2 M COP".
fcop <- function(x) paste0(format(round(x), big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE), " COP")
fmm <- function(x, d = 1) paste0(fnum(x / 1e6, d), " M COP")

# Fecha en formato AAAA-MM-DD.
ffecha <- function(x) format(x, "%Y-%m-%d")

# Polinomio c(B) = 1 + c1 B + ... + cp B^p en notacion LaTeX para el informe. coef = c(c1, ..., cp).
texto_polinomio <- function(coef) {
  out <- "1"
  for (k in seq_along(coef)) {
    if (coef[k] == 0) next
    a <- abs(coef[k])
    mono <- if (k == 1) "B" else paste0("B^{", k, "}")
    out <- paste0(out, if (coef[k] < 0) " - " else " + ", if (a == 1) "" else format(a), mono)
  }
  paste0("$", out, "$")
}

# Texto de parametros para el informe: letras griegas y coma decimal ("alpha = 0.5" -> "α = 0,5").
bonito <- function(x) {
  x <- gsub("alpha", "α", x, fixed = TRUE)
  x <- gsub("beta", "β", x, fixed = TRUE)
  x <- gsub("gamma", "γ", x, fixed = TRUE)
  gsub("([0-9])\\.([0-9])", "\\1,\\2", x)
}

# Marca los parametros elegidos que quedan en un limite de la rejilla (lectura de diagnostico de las notas:
# un valor pegado al borde sugiere que la familia del metodo no es la adecuada o que la rejilla es corta).
en_borde <- function(opt) {
  if (is.null(opt)) return("")
  pars <- names(opt$curva)[names(opt$curva) != "mse"]
  marcas <- vapply(pars, function(p) {
    v <- opt$mejor[[p]]
    rng <- range(opt$curva[[p]])
    nombre <- bonito(p)
    if (isTRUE(all.equal(v, rng[1]))) paste0(nombre, ": límite inferior") else if (isTRUE(all.equal(v, rng[2]))) paste0(nombre, ": límite superior") else ""
  }, character(1))
  paste(marcas[marcas != ""], collapse = "; ")
}

# Nombres abreviados de los meses (la serie 1 empieza en enero, la serie 3 en marzo).
meses_abrev <- function() c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")
