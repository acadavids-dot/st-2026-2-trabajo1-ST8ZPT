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

# --- VERIFICACION: metodos simples contra ejemplos hechos a mano ---
igual_na <- function(a, b) isTRUE(all.equal(a, b)) && identical(is.na(a), is.na(b))
verificar(igual_na(ajustar_media(1:6)$yhat, c(NA, 1, 1.5, 2, 2.5, 3)) && all(ajustar_media(1:6)$pronosticar(3) == 3.5),
          "ajustar_media: y = 1..6 da yhat = (NA, 1, 1.5, 2, 2.5, 3) y pronostico 3.5")
verificar(igual_na(ajustar_ingenuo(1:6)$yhat, c(NA, 1:5)) && all(ajustar_ingenuo(1:6)$pronosticar(2) == 6),
          "ajustar_ingenuo: yhat[t] = y[t-1] y pronostico = ultimo valor")
ie <- ajustar_ingenuo_estacional(c(5, 7, 9, 6, 8, 10), 3)
verificar(igual_na(ie$yhat, c(NA, NA, NA, 5, 7, 9)) && all(ie$pronosticar(4) == c(6, 8, 10, 6)),
          "ajustar_ingenuo_estacional (s = 3): yhat[t] = y[t-3] y repite el ultimo ciclo")
mm2 <- ajustar_mm(c(2, 4, 6, 8, 10), 2)
verificar(igual_na(mm2$yhat, c(NA, NA, 3, 5, 7)) && all(mm2$pronosticar(2) == 9),
          "ajustar_mm (k = 2): yhat = (NA, NA, 3, 5, 7) y pronostico = media de los ultimos 2 = 9")

# --- VERIFICACION: tendencias recuperan exactamente series construidas sin ruido ---
tt <- 1:40
fl <- ajustar_tendencia(3 + 2 * tt, "lineal")
verificar(max(abs(fl$yhat - (3 + 2 * tt)), na.rm = TRUE) < 1e-8 && max(abs(fl$pronosticar(3) - (3 + 2 * 41:43))) < 1e-8,
          "ajustar_tendencia lineal recupera y = 3 + 2t (yhat y pronostico exactos)")
fc <- ajustar_tendencia(1 + tt + 0.5 * tt^2, "cuadratica")
verificar(max(abs(fc$pronosticar(2) - (1 + 41:42 + 0.5 * (41:42)^2))) < 1e-6,
          "ajustar_tendencia cuadratica recupera y = 1 + t + 0.5 t^2")
fe <- ajustar_tendencia(2 * exp(0.05 * tt), "exponencial")
verificar(max(abs(fe$pronosticar(2) / (2 * exp(0.05 * 41:42)) - 1)) < 1e-8,
          "ajustar_tendencia exponencial recupera y = 2 exp(0.05 t)")
S_ver <- c(-2, 1, 1)
fs <- ajustar_tendencia_estacional(10 + 0.5 * tt + rep(S_ver, length.out = 40)[1:40], 3)
esperado <- 10 + 0.5 * 41:43 + rep(S_ver, length.out = 43)[41:43]
verificar(max(abs(fs$pronosticar(3) - esperado)) < 1e-8,
          "ajustar_tendencia_estacional recupera tendencia + efectos de suma cero (s = 3)")

# --- VERIFICACION: contrato (yhat, pronosticar, parametros), sin mirar el futuro y validacion de entradas ---
set.seed(5)
y_c <- 100 + cumsum(rnorm(60))
metodos_c <- list(
  media = function(y) ajustar_media(y), ingenuo = function(y) ajustar_ingenuo(y),
  ing_est = function(y) ajustar_ingenuo_estacional(y, 4), mm = function(y) ajustar_mm(y, 5),
  lineal = function(y) ajustar_tendencia(y, "lineal"), cuad = function(y) ajustar_tendencia(y, "cuadratica"),
  expo = function(y) ajustar_tendencia(y, "exponencial"), est = function(y) ajustar_tendencia_estacional(y, 4))
t0 <- 40
y_mod <- y_c
y_mod[t0] <- y_mod[t0] + 50  # cambiar y[t0] no debe alterar yhat[1..t0]
ok_contrato <- vapply(metodos_c, function(f) {
  a <- f(y_c); b <- f(y_mod)
  identical(names(a), c("yhat", "pronosticar", "parametros")) && length(a$yhat) == length(y_c) &&
    is.function(a$pronosticar) && is.list(a$parametros) && length(a$pronosticar(5)) == 5 &&
    igual_na(a$yhat[1:t0], b$yhat[1:t0]) && !isTRUE(all.equal(a$yhat[(t0 + 1):(t0 + 5)], b$yhat[(t0 + 1):(t0 + 5)]))
}, logical(1))
verificar(all(ok_contrato), "los 8 metodos cumplen el contrato y yhat[1..t] no depende de y[t] ni del futuro")
ok_valid <- all(vapply(list(
  function() ajustar_mm(c(1, NA, 3, 4), 2), function() ajustar_mm(1:10, 1), function() ajustar_media("a"),
  function() ajustar_ingenuo_estacional(1:3, 3), function() ajustar_tendencia(c(1, -1, 2, 3, 4, 5), "exponencial")),
  function(f) inherits(try(f(), silent = TRUE), "try-error"), logical(1)))
verificar(ok_valid, "los metodos rechazan entradas invalidas (NA, k fuera de rango, no numerico, y <= 0 en exponencial)")

# --- VERIFICACION: optimizar() y origen_movil() con casos de resultado conocido ---
y_o <- c(1, 3, 2, 4, 3, 5, 4, 6, 99)  # el 99 queda fuera del tramo de estimacion (n0 = 8): no debe influir
solo_k2 <- optimizar(y_o, ajustar_mm, data.frame(k = 2), 8)
k2_k3 <- optimizar(y_o, ajustar_mm, data.frame(k = 2:3), 8)
# k = 2 solo: errores en t = 3..8 son (0, 1.5, 0, 1.5, 0, 1.5) -> MSE = 6.75/6 = 1.125
# con k = 2:3 la muestra comun es t = 4..8 -> errores (1.5, 0, 1.5, 0, 1.5) -> MSE(k = 2) = 6.75/5 = 1.35
verificar(abs(solo_k2$mejor$mse - 1.125) < 1e-12 && abs(k2_k3$curva$mse[1] - 1.35) < 1e-12 && k2_k3$n_comun == 5,
          "optimizar: MSE a mano (1.125 con k = 2; 1.35 sobre la muestra comun t = 4..8 con k = 2:3)")
om <- origen_movil(1:10, ajustar_ingenuo, list(), n0 = 6, H = 3)
verificar(all(errores_horizonte(om, 1)$e == 1) && length(errores_horizonte(om, 1)$e) == 4 &&
            all(errores_horizonte(om, 3)$e == 3) && length(errores_horizonte(om, 3)$e) == 2 &&
            identical(errores_horizonte(om, 3)$origenes, 6:7),
          "origen_movil: ingenuo sobre y = 1..10 da error h en el horizonte h; 4 origenes en h = 1 y 2 en h = 3 (T - H - n0 + 1)")
verificar(isTRUE(all.equal(unname(unlist(medidas(c(1, -2, 3), c(10, 20, 30))[, c("rmse", "mad", "mape")])),
                           c(sqrt(14 / 3), 2, 10))),
          "medidas: e = (1, -2, 3), real = (10, 20, 30) da RMSE = sqrt(14/3), MAD = 2, MAPE = 10 %")

# --- VERIFICACION: SES, DMM, Holt y Holt-Winters contra cifras de las notas de clase (clases 4 y 5) ---
# SES a mano: y = (10, 20, 30), alpha = 0.5 -> yhat = (NA, 10, 15) y Yhat_4 = 0.5*30 + 0.5*15 = 22.5
ses_h <- ajustar_ses(c(10, 20, 30), 0.5)
verificar(igual_na(ses_h$yhat, c(NA, 10, 15)) && abs(ses_h$pronosticar(2)[1] - 22.5) < 1e-12,
          "ajustar_ses: ejemplo a mano (yhat = NA, 10, 15; pronostico 22.5)")
# Creditos al consumo (notas, clase 4 diap. 66 y clase 5 diap. 47), 11 meses
cred <- c(133, 155, 165, 171, 194, 231, 274, 312, 313, 333, 343)
ho <- ajustar_holt(cred, 0.3, 0.4)
L_notas <- c(133, 139.600, 149.068, 159.407, 174.936, 199.194, 232.891, 272.812, 305.759, 335.691, 359.321)
T_notas <- c(0, 2.640, 5.371, 7.358, 10.627, 16.079, 23.126, 29.844, 31.085, 30.624, 27.826)
yh_notas <- c(NA, 133, 142.240, 154.439, 166.766, 185.563, 215.273, 256.017, 302.656, 336.845, 366.315)
e_ho <- cred - ho$yhat
verificar(max(abs(ho$parametros$nivel_t - L_notas)) < 1e-3 && max(abs(ho$parametros$pendiente_t - T_notas)) < 1e-3 &&
            max(abs(ho$yhat - yh_notas), na.rm = TRUE) < 1e-3,
          "ajustar_holt(0.3, 0.4) reproduce las columnas L, T y yhat de la tabla de las notas (creditos al consumo)")
verificar(sum(!is.na(e_ho)) == 10 && abs(mean(e_ho^2, na.rm = TRUE) - 1133.08) < 0.01 && abs(mean(abs(e_ho), na.rm = TRUE) - 28.62) < 0.01 &&
            max(abs(ho$pronosticar(3) - c(387.15, 414.97, 442.80))) < 0.01,
          "ajustar_holt: N = 10, MSE = 1133.08, MAD = 28.62 y pronosticos 387.15, 414.97, 442.80 como en las notas")
dm3 <- ajustar_dmm(cred, 3)
e_dm <- cred - dm3$yhat
verificar(sum(!is.na(e_dm)) == 6 && abs(mean(e_dm^2, na.rm = TRUE) - 841.24) < 0.01 && abs(mean(abs(e_dm), na.rm = TRUE) - 28.20) < 0.01 &&
            max(abs(dm3$pronosticar(3) - c(356.56, 370.00, 383.44))) < 0.01,
          "ajustar_dmm(3): N = 6, MSE = 841.24, MAD = 28.20 y pronosticos 356.56, 370.00, 383.44 como en las notas")
# Ventas trimestrales (notas, clase 5 diap. 54): Holt-Winters aditivo (0.3, 0.05, 0.3), paso t = 5
ventas <- c(1.5, 4.3, 8.9, 12.5, 1.8, 4.5, 9.0, 12.1, 1.2, 4.6, 8.1, 12.9, 1.9, 5.1, 8.6, 13.2)
hw_v <- ajustar_hw(ventas, 4, 0.3, 0.05, 0.3)
verificar(abs(hw_v$parametros$b0 - 0.4 / 12) < 1e-12 && abs(hw_v$parametros$L0 - 6.75) < 1e-12 && abs(sum(hw_v$parametros$S0)) < 1e-12,
          "ajustar_hw: inicializacion a mano (b0 = 0.4/12, L0 = 6.8 - 1.5 b0 = 6.75, indices con suma cero hasta 1e-12)")
verificar(abs(hw_v$yhat[5] - 1.6505) < 5e-5 && abs((ventas[5] - hw_v$yhat[5]) - 0.1495) < 5e-5,
          "ajustar_hw: el pronostico del paso t = 5 es 1.6505 y su error 0.1495, como en 'Un periodo al detalle' de las notas")
# Serie aditiva exacta sin ruido: la inicializacion recupera exactamente (L0, b0, S0). Las notas usan (L0, b0) como estado
# PREVIO a t = 1, asi que yhat[1] = L0 + b0 (el nivel de t = 2) y el primer error es -b0; la recursion lo va absorbiendo.
y_hw_exacta <- 100 + 0.5 * (1:48) + rep(c(-3, 1, 5, -2, 4, -5), 8)
hw_e <- ajustar_hw(y_hw_exacta, 6, 0.4, 0.2, 0.3)
verificar(abs(hw_e$parametros$L0 - 100.5) < 1e-9 && abs(hw_e$parametros$b0 - 0.5) < 1e-9 &&
            max(abs(hw_e$parametros$S0 - c(-3, 1, 5, -2, 4, -5))) < 1e-9 && abs((y_hw_exacta[1] - hw_e$yhat[1]) + 0.5) < 1e-9,
          "ajustar_hw: sobre una serie aditiva exacta, la inicializacion es exacta y el primer error es -b0 (convencion de las notas)")
verificar(abs(y_hw_exacta[48] - hw_e$yhat[48]) < 0.05 &&
            max(abs(hw_e$pronosticar(6) - (100 + 0.5 * 49:54 + c(-3, 1, 5, -2, 4, -5)))) < 0.05,
          "ajustar_hw: la recursion absorbe el desfase inicial (error final < 0.05 en yhat[48] y en el pronostico a 6 pasos)")
y_mult_ver <- (50 + 0.4 * (1:60)) * rep(c(0.8, 1.1, 1.3, 0.9, 1.2, 0.7), 10)
hm <- ajustar_hw_mult(y_mult_ver, 6, 0.2, 0.1, 0.2)
verif_mult <- (50 + 0.4 * 61:66) * c(0.8, 1.1, 1.3, 0.9, 1.2, 0.7)
verificar(abs(mean(hm$parametros$S0) - 1) < 1e-12 && max(abs(hm$pronosticar(6) / verif_mult - 1)) < 0.05,
          "ajustar_hw_mult: indices iniciales con media 1 y pronostico con error relativo < 5 % sobre una serie multiplicativa sin ruido")
ok_valid2 <- all(vapply(list(
  function() ajustar_ses(1:5, 1), function() ajustar_holt(1:5, 0.5, 0), function() ajustar_dmm(1:5, 3),
  function() ajustar_hw(1:10, 6, 0.3, 0.1, 0.1), function() ajustar_hw_mult(c(1:11, -1), 6, 0.3, 0.1, 0.1)),
  function(f) inherits(try(f(), silent = TRUE), "try-error"), logical(1)))
verificar(ok_valid2, "SES, Holt, DMM y HW rechazan constantes fuera de (0,1), series cortas (< 2k o < 2 ciclos) y y <= 0 en el multiplicativo")
metodos_c2 <- list(ses = function(y) ajustar_ses(y, 0.3), dmm = function(y) ajustar_dmm(y, 4), holt = function(y) ajustar_holt(y, 0.3, 0.2))
ok_contrato2 <- vapply(metodos_c2, function(f) {
  a <- f(y_c); b <- f(y_mod)
  identical(names(a), c("yhat", "pronosticar", "parametros")) && length(a$yhat) == length(y_c) && length(a$pronosticar(5)) == 5 &&
    igual_na(a$yhat[1:t0], b$yhat[1:t0]) && !isTRUE(all.equal(a$yhat[(t0 + 1):(t0 + 5)], b$yhat[(t0 + 1):(t0 + 5)]))
}, logical(1))
verificar(all(ok_contrato2), "SES, DMM y Holt cumplen el contrato y su yhat[1..t] no depende de y[t] ni del futuro")
hw_c <- ajustar_hw(y_c, 4, 0.3, 0.1, 0.2); hwm_c <- ajustar_hw_mult(y_c, 4, 0.3, 0.1, 0.2)
verificar(all(vapply(list(hw_c, hwm_c), function(a) identical(names(a), c("yhat", "pronosticar", "parametros")) &&
                       length(a$yhat) == length(y_c) && length(a$pronosticar(9)) == 9 && !anyNA(a$yhat), logical(1))),
          "Holt-Winters (aditivo y multiplicativo) cumplen el contrato; su yhat inicial usa los ciclos completos, como las notas")

# --- VERIFICACION: numero de origenes sobre las series reales contra la tabla de particion ---
H_ver <- c(12, 9, 6)
esperado_orig <- list(c(36, 25), c(36, 28), c(26, 21))
ok_orig <- vapply(1:3, function(i) {
  om_i <- origen_movil(d_ver[[i]], ajustar_ingenuo, list(), n0 = floor(0.8 * length(d_ver[[i]])), H = H_ver[i])
  identical(c(length(errores_horizonte(om_i, 1)$e), length(errores_horizonte(om_i, H_ver[i])$e)), as.integer(esperado_orig[[i]]))
}, logical(1))
verificar(all(ok_orig), "origen_movil sobre las series reales: origenes (h = 1, h = H) = (36, 25), (36, 28), (26, 21) como la tabla de particion")
# validar_errores() corre de punta a punta y devuelve los cinco componentes con los seis elementos
ve <- validar_errores(errores_horizonte(origen_movil(d_ver[[2]], ajustar_ingenuo, list(), 144, 9), 1)$e, 10)
verificar(all(c("media", "lb", "jb", "dw") %in% names(ve)) && ve$n == 36 && length(imprimir_prueba(ve$lb)) == 8 &&
            abs(ve$banda_rb - 1.96 / 6) < 1e-12,
          "validar_errores: devuelve t de media cero, Ljung-Box, Jarque-Bera y DW con la banda sobre n = 36 errores")

# --- VERIFICACION: mco_hac() contra lm() y contra el calculo explicito de Newey-West ---
yy <- d_ver[[3]]
Xv <- matriz_diseno(length(yy), "lineal")
mh <- mco_hac(yy, Xv)
lm_v <- lm(yy ~ Xv[, "t"])
verificar(max(abs(mh$coeficientes - as.numeric(coef(lm_v)))) < 1e-8 &&
            max(abs(mh$tabla$ee_ordinario - summary(lm_v)$coefficients[, "Std. Error"])) < 1e-8,
          "mco_hac: coeficientes y error estandar ordinario coinciden con lm()")
# Newey-West con bucles explicitos sobre productos x_t x_s' (otro camino que el matricial)
nw_bucle <- function(y, X, L) {
  n <- length(y); k <- ncol(X)
  b <- solve(t(X) %*% X, t(X) %*% y); e <- as.numeric(y - X %*% b)
  S <- matrix(0, k, k)
  for (t in 1:n) S <- S + e[t]^2 * X[t, ] %*% t(X[t, ])
  for (j in seq_len(L)) for (t in (j + 1):n) {
    Gt <- e[t] * e[t - j] * (X[t, ] %*% t(X[t - j, ]))
    S <- S + (1 - j / (L + 1)) * (Gt + t(Gt))
  }
  XtXi <- solve(t(X) %*% X)
  sqrt(diag(XtXi %*% S %*% XtXi))
}
verificar(max(abs(mh$tabla$ee_hac - nw_bucle(yy, Xv, mh$L))) < 1e-10 &&
            max(abs(mco_hac(yy, Xv, L = 0)$tabla$ee_hac - nw_bucle(yy, Xv, 0))) < 1e-10,
          sprintf("mco_hac: error estandar HAC coincide con Newey-West a bucles (L = %d y L = 0, este ultimo es White)", mh$L))

# --- VERIFICACION: pruebas t de media cero, Jarque-Bera y Durbin-Watson ---
e_v <- z_ver$S2
tm <- prueba_t_media(e_v)
verificar(abs(tm$observado - as.numeric(t.test(e_v)$statistic)) < 1e-10 && abs(tm$valor_p - t.test(e_v)$p.value) < 1e-10,
          "prueba_t_media coincide con t.test()")
# x = (-2, -1, 0, 1, 2): m2 = 2, m3 = 0, m4 = 6.8, K = 1.7 -> JB = (5/6)(0 + (1.7 - 3)^2/4) = 0.3520833
jb_h <- jarque_bera(c(-2, -1, 0, 1, 2))
verificar(abs(jb_h$observado - 5 / 6 * (1.3^2 / 4)) < 1e-12 && abs(jb_h$asimetria) < 1e-12 && abs(jb_h$curtosis - 1.7) < 1e-12,
          "jarque_bera: ejemplo a mano x = (-2,...,2) da S = 0, K = 1.7 y JB = 0.3521")
e_dw <- z_ver$S3
dw <- durbin_watson(e_dw)
verificar(abs(dw$d - sum(diff(e_dw)^2) / sum(e_dw^2)) < 1e-12 && abs(dw$d - dw$aprox_2_1_menos_r1) < 0.05,
          sprintf("durbin_watson: coincide con su definicion via diff() y con 2(1 - r1) (d = %.3f, 2(1-r1) = %.3f)", dw$d, dw$aprox_2_1_menos_r1))
