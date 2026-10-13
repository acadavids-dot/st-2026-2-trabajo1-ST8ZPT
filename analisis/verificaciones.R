# verificaciones.R — contraste de las implementaciones propias contra R y contra ejemplos a mano.
# Unico lugar (junto a bloques marcados) donde se permiten acf(), pacf(), lm(), Box.test() y diff().
# Se carga desde correr-todo.R despues de source() de R/. Si algo falla, se detiene con stop().

verificar <- function(cond, msg) {
  if (!isTRUE(cond)) stop("VERIFICACION FALLIDA: ", msg, call. = FALSE)
  cat("OK:", msg, "\n")
}

# --- VERIFICACION: acf_muestral() contra stats::acf() ---
set.seed(1)
y_prueba <- cumsum(rnorm(200)) + sin(seq_len(200) / 5)
dif_acf <- max(abs(acf_muestral(y_prueba, 24) - as.numeric(acf(y_prueba, lag.max = 24, plot = FALSE)$acf)[-1]))
verificar(dif_acf < 1e-10, sprintf("acf_muestral coincide con acf() (diferencia maxima %.2e)", dif_acf))

# --- VERIFICACION: media_movil_centrada() contra un calculo a mano ---
# y = 1..10, s = 4 (2x4): en t = 3 -> (0.5*1 + 2 + 3 + 4 + 0.5*5)/4 = 3. Una recta se conserva: tendencia[t] = t.
mm <- media_movil_centrada(1:10, 4)
verificar(isTRUE(all.equal(mm[3:8], 3:8)) && all(is.na(mm[c(1, 2, 9, 10)])),
          "media_movil_centrada (s par) reproduce una recta y deja NA en los extremos")
# s impar: y = 1..9, s = 3 -> promedio de (t-1, t, t+1) = t
mm3 <- media_movil_centrada(1:9, 3)
verificar(isTRUE(all.equal(mm3[2:8], 2:8)) && all(is.na(mm3[c(1, 9)])),
          "media_movil_centrada (s impar) reproduce una recta y deja NA en los extremos")

# --- VERIFICACION: descomposicion_clasica() recupera los indices de una serie construida ---
S_real <- c(-3, 1, 5, -2, 4, -5)            # suma 0
y_aditiva <- 100 + 0.5 * seq_len(60) + rep(S_real, 10)
da <- descomposicion_clasica(y_aditiva, 6, "aditiva")
verificar(max(abs(da$indices - S_real)) < 1e-10 && abs(sum(da$indices)) < 1e-10,
          "descomposicion aditiva recupera los indices (suma 0) de una serie sin ruido")
M_real <- c(0.8, 1.1, 1.3, 0.9, 1.2, 0.7)   # media 1
# Con nivel constante la recuperacion es exacta. Con tendencia lineal el metodo clasico tiene un sesgo
# pequeno (la media movil de T*S no es exactamente T), por eso ahi se usa una tolerancia de 0.01.
dm0 <- descomposicion_clasica(100 * rep(M_real, 10), 6, "multiplicativa")
verificar(max(abs(dm0$indices - M_real)) < 1e-10 && abs(mean(dm0$indices) - 1) < 1e-10,
          "descomposicion multiplicativa recupera los indices (media 1) con nivel constante")
dm <- descomposicion_clasica((50 + 0.4 * seq_len(60)) * rep(M_real, 10), 6, "multiplicativa")
verificar(max(abs(dm$indices - M_real)) < 0.01 && abs(mean(dm$indices) - 1) < 1e-10,
          "descomposicion multiplicativa con tendencia: error de indices < 0.01 y media 1")

# --- VERIFICACION: diferenciar() contra diff() ---
combos <- list(c(1, 0, 1), c(2, 0, 1), c(0, 1, 12), c(1, 1, 12), c(0, 1, 7), c(1, 1, 7), c(2, 1, 6), c(0, 2, 4))
ok_dif <- vapply(combos, function(p) {
  d <- p[1]; D <- p[2]; s <- p[3]
  ref <- y_prueba
  if (D > 0) for (i in seq_len(D)) ref <- diff(ref, lag = s)
  if (d > 0) ref <- diff(ref, differences = d)
  length(ref) == length(y_prueba) - d - D * s && max(abs(diferenciar(y_prueba, d, D, s) - ref)) < 1e-10
}, logical(1))
verificar(all(ok_dif), sprintf("diferenciar coincide con diff() en %d combinaciones de (d, D, s)", length(combos)))

# --- VERIFICACION: polinomio_rezago() y coef_diferencia() contra casos conocidos ---
r1 <- polinomio_rezago(coef_diferencia(1, 0, 1))
verificar(nrow(r1) == 1 && abs(r1$re - 1) < 1e-10 && abs(r1$modulo - 1) < 1e-10, "1 - B tiene una raiz en B = 1")
r4 <- polinomio_rezago(coef_diferencia(0, 1, 4))
verificar(nrow(r4) == 4 && all(abs(r4$modulo - 1) < 1e-10), "1 - B^4 tiene 4 raices de modulo 1 (+-1, +-i)")
verificar(identical(coef_diferencia(2, 0, 1), c(-2, 1)) && identical(coef_diferencia(1, 1, 2), c(-1, -1, 1)),
          "coef_diferencia: (1-B)^2 = 1 - 2B + B^2 y (1-B)(1-B^2) = 1 - B - B^2 + B^3")
# Aplicar el polinomio con sus coeficientes (y_t + c1 y_{t-1} + ...) debe dar lo mismo que diferenciar().
aplicar_coef <- function(y, cf) {
  p <- length(cf)
  n <- length(y)
  sapply((p + 1):n, function(t) y[t] + sum(cf * y[(t - 1):(t - p)]))
}
ok_cf <- all(vapply(combos, function(p) {
  max(abs(aplicar_coef(y_prueba, coef_diferencia(p[1], p[2], p[3])) - diferenciar(y_prueba, p[1], p[2], p[3]))) < 1e-8
}, logical(1)))
verificar(ok_cf, "coef_diferencia aplicado como polinomio coincide con diferenciar()")

# --- VERIFICACION: media de la serie diferenciada y identidad Var(grad Y) = 2 gamma(0) - 2 gamma(1) ---
# La suma de (1-B)y telescopa: mean(grad y) = (y_n - y_1)/(n - 1) exactamente.
verificar(abs(mean(diferenciar(y_prueba, 1, 0, 1)) - (y_prueba[200] - y_prueba[1]) / 199) < 1e-12,
          "media de grad y = (y_n - y_1)/(n-1) (telescopa)")
set.seed(2)
y_ar <- numeric(3000)
for (t in 2:3000) y_ar[t] <- 0.5 * y_ar[t - 1] + rnorm(1)
id <- identidad_var_diferencia(y_ar)
verificar(abs(id$var_diferencia / id$dos_g0_menos_dos_g1 - 1) < 0.01,
          sprintf("identidad de Var(grad Y): %.4f vs %.4f (AR(1) estacionario, n = 3000)", id$var_diferencia, id$dos_g0_menos_dos_g1))

# --- VERIFICACION: pacf_dl() contra pacf() (maxima diferencia absoluta < 1e-10) ---
# Series: AR(2) simulado, y las tres series transformadas elegidas (se calculan aqui directamente).
set.seed(3)
y_ar2 <- numeric(500)
for (t in 3:500) y_ar2[t] <- 0.6 * y_ar2[t - 1] - 0.3 * y_ar2[t - 2] + rnorm(1)
d_ver <- lapply(1:3, function(i) leer_serie(sprintf("datos/ST8ZPT-serie-%d.csv", i))$valor)
z_ver <- list(AR2 = y_ar2,
              S1 = diferenciar(log(d_ver[[1]]), 1, 1, 12),
              S2 = diferenciar(d_ver[[2]], 1, 0, 1),
              S3 = residuos_tendencia(d_ver[[3]], 1)$residuos)
dif_pacf <- vapply(z_ver, function(z) {
  max(abs(pacf_dl(acf_muestral(z, 24)) - as.numeric(pacf(z, lag.max = 24, plot = FALSE)$acf)))
}, numeric(1))
verificar(all(dif_pacf < 1e-10), sprintf("pacf_dl coincide con pacf() (diferencia maxima %.2e en AR2, S1, S2, S3)", max(dif_pacf)))
# A mano: con r1 = 0.5, r2 = 0.2 -> phi_22 = (0.2 - 0.5*0.5)/(1 - 0.5*0.5) = -0.05/0.75
verificar(abs(pacf_dl(c(0.5, 0.2))[2] - (-0.05 / 0.75)) < 1e-12 && abs(pacf_dl(c(0.5, 0.2))[1] - 0.5) < 1e-12,
          "pacf_dl: ejemplo a mano r1 = 0.5, r2 = 0.2 da phi_11 = 0.5 y phi_22 = -1/15")
# Teoria: r_k = rho^k (AR(1)) -> pacf = (rho, 0, 0, ...)
verificar(max(abs(pacf_dl(0.7^(1:6))[-1])) < 1e-12 && abs(pacf_dl(0.7^(1:6))[1] - 0.7) < 1e-12,
          "pacf_dl: r_k = 0.7^k da phi_11 = 0.7 y phi_kk = 0 para k >= 2")

# --- VERIFICACION: ljung_box() contra Box.test(type = 'Ljung-Box') ---
lb_propia <- ljung_box(z_ver$S2, 24)
lb_r <- Box.test(z_ver$S2, lag = 24, type = "Ljung-Box")
verificar(abs(lb_propia$observado - as.numeric(lb_r$statistic)) < 1e-8 && abs(lb_propia$valor_p - lb_r$p.value) < 1e-10,
          sprintf("ljung_box coincide con Box.test: Q = %.4f, p = %.4f", lb_propia$observado, lb_propia$valor_p))

# --- VERIFICACION: banda de Bartlett y varianza de largo plazo contra calculos a mano ---
bb <- banda_bartlett(c(0.5, 0.3, 0.1), 100)
verificar(abs(bb[1] - 1.96 / 10) < 1e-12 && abs(bb[2] - 1.96 * sqrt((1 + 2 * 0.25) / 100)) < 1e-12 &&
            abs(bb[3] - 1.96 * sqrt((1 + 2 * (0.25 + 0.09)) / 100)) < 1e-12,
          "banda_bartlett: h = 1 es 1.96/sqrt(n); h = 2 y 3 coinciden con el calculo a mano")
# y = (1, 3, 2, 4): media 2.5, g0 = 1.25, g1 = -0.4375; con L = 1: g0 + 2 (1 - 1/2) g1 = 0.8125
verificar(abs(varianza_largo_plazo(c(1, 3, 2, 4), 1) - 0.8125) < 1e-12 && abs(varianza_largo_plazo(c(1, 3, 2, 4), 0) - 1.25) < 1e-12,
          "varianza_largo_plazo: ejemplo a mano y = (1,3,2,4) da 1.25 (L = 0) y 0.8125 (L = 1)")

# --- VERIFICACION: residuos_tendencia() contra lm() y bloques ---
tt <- seq_along(d_ver[[3]])
ajuste_lm <- lm(d_ver[[3]] ~ tt)
verificar(max(abs(residuos_tendencia(d_ver[[3]], 1)$coeficientes - as.numeric(coef(ajuste_lm)))) < 1e-8 &&
            max(abs(residuos_tendencia(d_ver[[3]], 1)$residuos - as.numeric(resid(ajuste_lm)))) < 1e-8 &&
            abs(mean(residuos_tendencia(d_ver[[3]], 1)$residuos)) < 1e-10,
          "residuos_tendencia coincide con lm() y su media es 0")
verificar(identical(sum(bloques_media_var(1:130, 4)$n), 130L) && all(bloques_media_var(1:130, 4)$n >= 32),
          "bloques_media_var reparte las 130 observaciones en 4 bloques consecutivos")
