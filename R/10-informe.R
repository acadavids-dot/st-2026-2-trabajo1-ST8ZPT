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

# Nombres abreviados de los meses (la serie 1 empieza en enero, la serie 3 en marzo).
meses_abrev <- function() c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")
