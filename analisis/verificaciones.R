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
