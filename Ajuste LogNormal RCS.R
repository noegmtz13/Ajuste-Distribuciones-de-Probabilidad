# Ajuste del RCS Trimestral mediante la distribucion LOGNORMAL

# install.packages(c("MASS","goftest"))
# MASS:
#   Contiene la función fitdistr(), utilizada para estimar
#   parámetros mediante Máxima Verosimilitud (MLE).

# goftest:
#   Contiene la prueba Anderson-Darling (ad.test),
#   utilizada para validar bondad de ajuste.

library(MASS)  
library(goftest)

# Evitar notación científica global. Esto permite visualizar cantidades grandes completas.
options(scipen = 999)

# Datos
# Vector con valores históricos trimestrales del RCS. Cada observación representa un trimestre histórico.
rcs <- c(
  1728451934.27, 1184927638.51, 1612375402.89, 1015832741.66,
  846928510.43, 879214365.77, 1087563210.92, 982145637.58,
  857339421.04, 944681275.39, 913472856.28, 1124983751.64,
  774195283.55, 835742918.36, 882164730.91, 841309552.18,
  668492375.82, 721843690.44, 903174856.37, 1265983471.25,
  1159243780.83, 1193764208.59, 1048756213.94, 1312475986.72,
  1819645328.41, 1715083476.66, 1542398705.18
)

# Estadísticos básicos
# Media aritmética de los datos históricos.
m    <- mean(rcs)

# Desviación estándar muestral.
# Mide la dispersión de los datos.
desv <- sd(rcs)

# Transformación logarítmica natural.
# Si X es Lognormal, entonces log(X) es Normal.
lr   <- log(rcs)

# Ajuste LOGNORMAL (MLE)
# Ajuste de distribución lognormal utilizando
# Máxima Verosimilitud (MLE).

# fitdistr estima:
#   meanlog -> media de log(X)
#   sdlog   -> desviación estándar de log(X)
ajus <- fitdistr(rcs, "log-normal")

# Mostrar resultados del ajuste.
print(ajus)

# Extraer parámetros estimados.
meanlog <- ajus$estimate["meanlog"]
sdlog   <- ajus$estimate["sdlog"]


# KS y AD 

# Prueba Kolmogorov-Smirnov (KS)
# Hipótesis nula: Los datos siguen una distribución Lognormal.
KS <- ks.test(rcs, "plnorm", meanlog = meanlog, sdlog = sdlog)

# Mostrar resultados KS.
print(KS)


# Prueba Anderson-Darling (AD)
# Se aplica sobre log(RCS).
# Si X ~ Lognormal: log(X) ~ Normal
AD <- ad.test(lr, null = "pnorm", mean = mean(lr), sd = sd(lr))

# Mostrar resultados AD.
print(AD)

# CDF de valores específicos y cuantiles/VaR
# Valores específicos para evaluar la función acumulada.
# Se calcula: F(x) = P(X <= x)
val_desv <- c(  1800000000, 1900000000, 2000000000,
                2100000000, 2200000000, 2300000000,
                2400000000, 2800000000, 3800000000,
                4800000000)

# Evaluar probabilidades acumuladas teóricas.
probas <- plnorm(val_desv, meanlog = meanlog, sdlog = sdlog)

# Exportar resultados a CSV.
write.csv(data.frame(x = val_desv, F = probas), "probas_lognormal.csv", row.names = FALSE)

# Vector de probabilidades para cuantiles/VaR.
# Ejemplo:
#   0.99  -> VaR 99%
#   0.995 -> VaR 99.5%
probas_2 <- c(seq(0.5,0.95,by=0.05), 0.97, 0.99, 0.995, 0.9997, 0.9999)

# Cálculo de cuantiles teóricos.
tabla_2  <- qlnorm(probas_2, meanlog = meanlog, sdlog = sdlog)

# Exportar tabla VaR.
write.csv(data.frame(p = probas_2, VaR = tabla_2), "tabla_VaR_lognormal.csv", row.names = FALSE)

# Gráfica: densidad empírica + teórica
# - Sin ejes visibles, solo números del eje X
# - 5 líneas verticales, con 1089 resaltada
# - Sin fondo coloreado

# Densidad empírica mediante Kernel Density Estimation.
d_emp  <- density(rcs, bw = "nrd0")                        # densidad empírica (gris)

# Malla de valores X para la curva teórica. Se evita iniciar en 0 porque la lognormal solo está definida para x > 0.
x_grid <- seq(1, 3000000000, length.out = 2000)            

# Evaluar densidad teórica lognormal.
f_theo <- dlnorm(x_grid, meanlog = meanlog, sdlog = sdlog) # densidad teórica (azul)

# Límite superior del eje Y.Se toma el máximo entre ambas curvas y se deja un margen adicional.
y_max <- max(c(d_emp$y, f_theo), na.rm = TRUE) * 1.10

# Gráfica inicial usando la densidad empírica.
# type = "l": Gráfica tipo línea.
# axes = FALSE: Oculta los ejes.
plot(d_emp$x, d_emp$y,
     type = "l", lwd = 1.5, col = "gray40",
     xlim = c(0, 3000000000),
     ylim = c(0, y_max),
     xlab = "", ylab = "", main = "",
     axes = FALSE)

# Añadir curva teórica sobre el gráfico.
lines(x_grid, f_theo, lwd = 2.2, col = "#2E75B6")

# Agregar únicamente etiquetas del eje X. Los valores se muestran en millones.
axis(1, at = c(114,439,764,1089,1414,1739,2063,2388)*1e6,
     labels = c("114","439","764","1,089","1,414","1,739","2,063","2,388"))

# Cinco líneas verticales: 1089 destacada
# Otras cuatro (tenues, punteadas)
vlines_tenues <- c(439, 764, 1414, 1739) * 1e6
abline(v = vlines_tenues, col = "#E07B7B", lty = 3, lwd = 1)

# Línea destacada en 1,089 millones
abline(v = 1089e6, col = "#C00000", lwd = 2.5, lty = 1)

# Marca corta para la media empírica.
segments(x0 = m, y0 = 0, y1 = y_max*0.03, lty = 1)  # media empírica

# Marca corta para la mediana teórica.
segments(x0 = qlnorm(0.5, meanlog = meanlog, sdlog = sdlog),
         y0 = 0, y1 = y_max*0.03, lty = 2, col = "darkgreen")

# Observaciones trimestrales
# Evaluar la densidad teórica en cada observación. Esto permite posicionar los puntos sobre la curva.
y_obs <- dlnorm(rcs, meanlog = meanlog, sdlog = sdlog)
points(rcs[seq(1,26,by=4)],     y_obs[seq(1,26,by=4)],     col="#C00000", pch=4, lwd=2)  # 1T (Marzo)
points(rcs[seq(1,26,by=4)+1],   y_obs[seq(1,26,by=4)+1],   col="#2E75B6", pch=4, lwd=2)  # 2T (Junio)
points(rcs[seq(1,26,by=4)+2],   y_obs[seq(1,26,by=4)+2],   col="#595959", pch=4, lwd=2)  # 3T (Sept)
points(rcs[seq(1,26,by=4)+3],   y_obs[seq(1,26,by=4)+3],   col="orange",  pch=4, lwd=2)  # 4T (Dic)

# (H) Punto proyectado Dic-25 sobre la densidad teórica
proy   <- 1425000000

# Evaluar densidad teórica del punto proyectado.
y_proy <- dlnorm(proy, meanlog = meanlog, sdlog = sdlog)

# Graficar punto proyectado.
points(proy, y_proy, col="orange", pch=18)

# Líneas “media ± k·sd” 
segments(x0 = m + desv,     y0 = 0, y1 = y_max*0.017, lty = 3, col = "red")
segments(x0 = m + 2*desv,   y0 = 0, y1 = y_max*0.017, lty = 3, col = "red")
segments(x0 = m - desv,     y0 = 0, y1 = y_max*0.017, lty = 3, col = "red")
segments(x0 = m - 2*desv,   y0 = 0, y1 = y_max*0.017, lty = 3, col = "red")

# (J) Leyenda
legend(x = 2.3e9, y = y_max*0.80,
       legend = c("1T (Marzo)","2T (Junio)","3T (Septiembre)","4T (Diciembre)",
                  "Densidad empírica","Densidad lognormal","Mediana"),
       col    = c("#C00000","#2E75B6","#595959","orange","gray40","#2E75B6","darkgreen"),
       pch    = c(4,4,4,4,NA,NA,NA),
       lty    = c(0,0,0,0,1,1,2),
       lwd    = c(2,2,2,2,1.5,2.2,2),
       bty    = "n")


# Métricas (loglik/AIC/BIC) -

# Número de observaciones.
n <- length(rcs)

# Número de parámetros estimados.La lognormal tiene: meanlog y sdlog
p <- 2  # parámetros de lognormal: meanlog, sdlog


# Función auxiliar para calcular: loglik, AIC y BIC. Esto se hace manualmente en caso de que la versión de MASS no devuelva dichas métricas.
get_loglik_aic_bic_ln <- function() {
  ll  <- tryCatch(ajus$loglik, error = function(e) NULL)
  aic <- tryCatch(ajus$aic,    error = function(e) NULL)
  # Si MASS no devuelve métricas:
  if (is.null(ll) || is.null(aic)) {
    # Densidad individual de cada observación.
    dens_i <- dlnorm(rcs, meanlog = meanlog, sdlog = sdlog)
    # Evitar log(0).
    eps    <- .Machine$double.xmin
    dens_i <- pmax(dens_i, eps)
    # Log-likelihood manual.
    ll_man <- sum(log(dens_i))
    # AIC manual.
    aic_man <- -2 * ll_man + 2 * p
    # BIC manual.
    bic_man <- -2 * ll_man + log(n) * p
    return(list(loglik = ll_man, AIC = aic_man, BIC = bic_man))
  } else {
    # Si MASS sí devolvió loglik y AIC, solo calcular BIC.
    bic <- -2 * ll + log(n) * p
    return(list(loglik = ll, AIC = aic, BIC = bic))
  }
}

# Ejecutar función de métricas.
met <- get_loglik_aic_bic_ln()

# Crear tabla resumen.
met_ln <- data.frame(
  modelo = "Lognormal",
  loglik = met$loglik,
  AIC    = met$AIC,
  BIC    = met$BIC
)

# Mostrar métricas.
cat("\n== Métricas de ajuste LOGNORMAL ==\n")
print(met_ln)


# VaR/TVaR 99.5% por simulación 

# Fijar semilla para reproducibilidad.
set.seed(12345)

# Simular 200,000 observaciones lognormales.
sim_ln <- rlnorm(2e5, meanlog = meanlog, sdlog = sdlog)

# Nivel de confianza.
alpha       <- 0.995

# VaR 99.5%.
VaR_ln_995  <- as.numeric(quantile(sim_ln, alpha))

# TVaR 99.5%.
# Promedio de pérdidas por encima del VaR.
TVaR_ln_995 <- mean(sim_ln[sim_ln > VaR_ln_995])

# Mostrar resultados.
cat("\n== VaR/TVaR 99.5% (LOGNORMAL, simulación) ==\n")
print(data.frame(modelo = "Lognormal", VaR_995 = VaR_ln_995, TVaR_995 = TVaR_ln_995))

# Exportar resumen a CSV

# Tabla consolidada con:
#   - parámetros
#   - métricas
#   - pruebas de ajuste
#   - VaR
#   - TVaR
resumen_ln <- data.frame(
  meanlog = meanlog, sdlog = sdlog,
  loglik  = met$loglik, AIC = met$AIC, BIC = met$BIC,
  KS_p    = KS$p.value,
  AD_p    = AD$p.value,
  VaR_99  = qlnorm(0.99,  meanlog = meanlog, sdlog = sdlog),
  VaR_995 = qlnorm(0.995, meanlog = meanlog, sdlog = sdlog),
  TVaR_995= TVaR_ln_995
)

# Exportar archivo CSV.
write.csv(resumen_ln, "resumen_lognormal_metrics.csv", row.names = FALSE)
cat("\n>> Guardado 'resumen_lognormal_metrics.csv'\n")

# Guardar imagen a PNG del Ajuste de Distribucion
png("grafica_lognormal_rcs.png", width = 1600, height = 900, res = 144)
dev.off()

