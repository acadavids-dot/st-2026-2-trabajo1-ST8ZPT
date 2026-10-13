# 07-regresion.R — MCO con error estandar ordinario y robusto a autocorrelacion (HAC, nucleo de Bartlett).
# Solo define funciones: no ejecuta nada al hacer source().

# MCO de y sobre X (X ya incluye el intercepto). Devuelve la tabla de coeficientes con:
#  - error estandar ordinario: sqrt(diag(s2 (X'X)^-1)), s2 = SCE/(n-k);
#  - error estandar HAC (Newey-West): V = (X'X)^-1 S (X'X)^-1, con
#    S = G0 + sum_{j=1}^{L} (1 - j/(L+1)) (Gj + Gj'), Gj = sum_{t>j} e_t e_{t-j} x_t x_{t-j}'.
#    Los pesos de Bartlett (1 - j/(L+1)) hacen que S sea semidefinida positiva. Sin correccion de muestra finita.
#  - prueba t con el error robusto (seis elementos) para H0: beta_j = 0, con n - k grados de libertad.
# L por defecto: floor(4 (n/100)^(2/9)), el mismo de la varianza de largo plazo.
mco_hac <- function(y, X, L = NULL, nombres = NULL) {
  stopifnot(is.numeric(y), !anyNA(y), is.matrix(X), nrow(X) == length(y), !anyNA(X))
  n <- length(y)
  k <- ncol(X)
  if (n <= k) stop("Se necesitan mas observaciones que coeficientes.")
  if (is.null(L)) L <- floor(4 * (n / 100)^(2 / 9))
  if (is.null(nombres)) nombres <- if (is.null(colnames(X))) paste0("b", seq_len(k) - 1) else colnames(X)
  b <- coef_mco(X, y)
  e <- as.numeric(y - X %*% b)
  XtXi <- solve(crossprod(X))
  se_ord <- sqrt(diag(sum(e^2) / (n - k) * XtXi))
  Z <- X * e  # fila t = e_t x_t'
  S <- crossprod(Z)
  for (j in seq_len(L)) {
    Gj <- crossprod(Z[(j + 1):n, , drop = FALSE], Z[1:(n - j), , drop = FALSE])
    S <- S + (1 - j / (L + 1)) * (Gj + t(Gj))
  }
  V <- XtXi %*% S %*% XtXi
  se_hac <- sqrt(diag(V))
  gl <- n - k
  crit <- stats::qt(0.975, gl)
  t_hac <- b / se_hac
  tabla <- tibble::tibble(termino = nombres, estimacion = b, ee_ordinario = se_ord, ee_hac = se_hac,
                          t_hac = t_hac, valor_p = 2 * stats::pt(abs(t_hac), gl, lower.tail = FALSE))
  pruebas <- lapply(seq_len(k), function(j) {
    nueva_prueba(
      nombre = paste0("Prueba t (error robusto HAC) para el coeficiente ", nombres[j]),
      h0 = paste0("beta_", nombres[j], " = 0"), h1 = paste0("beta_", nombres[j], " distinto de 0"),
      estadistico = "t", formula = "beta_hat / EE_HAC(beta_hat)",
      distribucion = paste0("t de Student con ", gl, " grados de libertad (aproximada)"),
      region = paste0("|t| > ", formatear_num(crit)), valor_critico = crit,
      observado = t_hac[j], valor_p = tabla$valor_p[j],
      lectura_rechaza = paste0("El termino ", nombres[j], " aporta a la explicacion de la serie."),
      lectura_no_rechaza = paste0("No hay evidencia de que el termino ", nombres[j], " sea distinto de 0.")
    )
  })
  names(pruebas) <- nombres
  list(tabla = tabla, pruebas = pruebas, coeficientes = b, residuos = e, ajustado = as.numeric(X %*% b),
       vcov_hac = V, L = L, n = n, k = k, gl = gl)
}
