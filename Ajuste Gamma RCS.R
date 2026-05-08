# Ajuste del RCS Trimestral mediante la distribucion GAMMA

# install.packages(c("fitdistrplus","goftest"))

# fitdistrplus:
#   Contiene la función fitdist(), utilizada para estimar
#   parámetros mediante Máxima Verosimilitud (MLE).

# goftest:
#   Contiene la prueba Anderson-Darling (ad.test),
#   utilizada para validar bondad de ajuste.

library(fitdistrplus)
library(goftest)

# Evitar notación científica global.
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


# Reescalamiento
#
# Gamma suele presentar problemas numéricos cuando
# los datos son extremadamente grandes.
#
# Por estabilidad computacional, los datos se
# reescalan temporalmente dividiendo entre 1e9.

rcs_scaled <- rcs / 1e9


# Estadísticos básicos

# Media aritmética.
m <- mean(rcs_scaled)

# Desviación estándar muestral.
desv <- sd(rcs_scaled)

# Varianza muestral.
varianza <- var(rcs_scaled)


# Ajuste GAMMA (MLE)

# Primero se obtienen estimadores iniciales mediante Método de Momentos.
# Para Gamma:
#   shape = mean^2 / var
#   rate  = mean / var

:contentReference[oaicite:0]{index=0}

shape_ini <- m^2 / varianza
rate_ini  <- m / varianza


# Ajuste Gamma mediante Máxima Verosimilitud (MLE).
#
# Los parámetros obtenidos por Método de Momentos
# se utilizan como valores iniciales del optimizador.

ajus <- fitdist(
  rcs_scaled,
  "gamma",
  start = list(
    shape = shape_ini,
    rate  = rate_ini
  )
)

# Mostrar resultados del ajuste.
print(ajus)

# Extraer parámetros estimados.
shape_gamma <- ajus$estimate["shape"]
rate_gamma  <- ajus$estimate["rate"]


# KS y AD

# Prueba Kolmogorov-Smirnov (KS)
# Hipótesis nula:
# Los datos siguen una distribución Gamma.

KS <- ks.test(rcs_scaled,"pgamma", shape = shape_gamma,rate = rate_gamma)

# Mostrar resultados KS.
print(KS)


# Prueba Anderson-Darling (AD)
# Se utiliza la transformación:
#   U = F(X)
# Si el ajuste Gamma es correcto:
#   U ~ Uniforme(0,1)

u_gamma <- pgamma(rcs_scaled, shape = shape_gamma, rate = rate_gamma)
AD <- ad.test(u_gamma, null = "punif")

# Mostrar resultados AD.
print(AD)


# CDF de valores específicos y cuantiles/VaR

# Valores específicos para evaluar la función acumulada.

val_desv <- c(
  1800000000, 1900000000, 2000000000, 2100000000, 2200000000, 2300000000,
  2400000000, 2800000000, 3800000000,4800000000)

# Reescalar valores.
val_desv_scaled <- val_desv / 1e9

# Evaluar probabilidades acumuladas Gamma.
probas <- pgamma(val_desv_scaled, shape = shape_gamma, rate = rate_gamma)

# Exportar resultados.
write.csv(data.frame(x = val_desv, F = probas), "probas_gamma.csv", row.names = FALSE)

# Vector de probabilidades para cuantiles/VaR.
probas_2 <- c(seq(0.5,0.95,by=0.05), 0.97, 0.99, 0.995, 0.9997, 0.9999)

# Cuantiles Gamma reescalados.
tabla_2_scaled <- qgamma(probas_2, shape = shape_gamma, rate = rate_gamma)

# Regresar a escala original.
tabla_2 <- tabla_2_scaled * 1e9

# Exportar tabla VaR.
write.csv(data.frame(p = probas_2, VaR = tabla_2), "tabla_VaR_gamma.csv", row.names = FALSE)


# Gráfica: densidad empírica + teórica

# Densidad empírica reescalada.
d_emp <- density(rcs_scaled, bw = "nrd0")

# Regresar eje X a escala original.
d_emp$x <- d_emp$x * 1e9

# Malla de valores X reescalados.
x_grid_scaled <- seq(0.001, 3, length.out = 2000)

# Densidad teórica Gamma.
f_theo <- dgamma(x_grid_scaled, shape = shape_gamma, rate = rate_gamma)

# Regresar X a escala original.
x_grid <- x_grid_scaled * 1e9

# Límite superior eje Y.
y_max <- max(c(d_emp$y, f_theo), na.rm = TRUE) * 1.10


# Gráfica inicial usando densidad empírica.

plot(
  d_emp$x, d_emp$y,
  type = "l", lwd = 1.5, col = "gray40",
  xlim = c(0, 3000000000),
  ylim = c(0, y_max),
  xlab = "", ylab = "", main = "",
  axes = FALSE)

# Añadir curva Gamma teórica.
lines(x_grid, f_theo, lwd = 2.2, col = "#2E75B6")

# Etiquetas eje X.
axis(1, at = c(114,439,764,1089,1414,1739,2063,2388)*1e6,
  labels = c("114","439","764","1,089","1,414","1,739","2,063","2,388"))

# Líneas verticales auxiliares.
vlines_tenues <- c(439,764,1414,1739) * 1e6
abline(v = vlines_tenues, col = "#E07B7B", lty = 3, lwd = 1)

# Línea destacada.
abline(v = 1089e6, col = "#C00000", lwd = 2.5)

# Marca media empírica.
segments(x0 = mean(rcs), y0 = 0, y1 = y_max*0.03)

# Marca mediana Gamma.

segments(x0 = qgamma(0.5,shape = shape_gamma, rate = rate_gamma) * 1e9,
        y0 = 0, y1 = y_max*0.03, lty = 2, col = "darkgreen")


# Observaciones trimestrales

# Evaluar densidad Gamma en datos reescalados.
y_obs <- dgamma(rcs_scaled, shape = shape_gamma, rate = rate_gamma)

# 1T (Marzo)
points(rcs[seq(1,26,by=4)], y_obs[seq(1,26,by=4)],  col="#C00000",  pch=4,  lwd=2)

# 2T (Junio)
points(rcs[seq(1,26,by=4)+1],  y_obs[seq(1,26,by=4)+1], col="#2E75B6",  pch=4, lwd=2)

# 3T (Septiembre)
points(rcs[seq(1,26,by=4)+2],  y_obs[seq(1,26,by=4)+2],  col="#595959",  pch=4,  lwd=2)

# 4T (Diciembre)
points(rcs[seq(1,26,by=4)+3],  y_obs[seq(1,26,by=4)+3],  col="orange",  pch=4,  lwd=2)


# Punto proyectado ilustrativo.
proy <- 1425000000

# Evaluar densidad del punto proyectado.
y_proy <- dgamma(proy / 1e9, shape = shape_gamma, rate = rate_gamma)

# Graficar punto proyectado.
points(proy, y_proy, col="orange", pch=18)


# Líneas “media ± k·sd”
segments(x0 = mean(rcs) + sd(rcs), y0 = 0, y1 = y_max*0.017, lty = 3, col = "red")

segments(x0 = mean(rcs) + 2*sd(rcs), y0 = 0, y1 = y_max*0.017, lty = 3, col = "red")

segments(x0 = mean(rcs) - sd(rcs), y0 = 0, y1 = y_max*0.017, lty = 3, col = "red")

segments(x0 = mean(rcs) - 2*sd(rcs), y0 = 0, y1 = y_max*0.017, lty = 3, col = "red")


# Leyenda
legend(x = 2.3e9, y = y_max*0.80,
  
  legend = c("1T (Marzo)","2T (Junio)", "3T (Septiembre)", "4T (Diciembre)",
          "Densidad empírica", "Densidad Gamma", "Mediana"),
  col = c("#C00000","#2E75B6","#595959","orange","gray40","#2E75B6","darkgreen"),
  pch = c(4,4,4,4,NA,NA,NA),
  lty = c(0,0,0,0,1,1,2),
  lwd = c(2,2,2,2,1.5,2.2,2),
  bty = "n")


# Métricas (loglik/AIC/BIC)

# Número de observaciones.
n <- length(rcs_scaled)

# Número de parámetros.
p <- 2

# Log-likelihood Gamma.

loglik_gamma <- sum(
  log(
    dgamma(
      rcs_scaled,
      shape = shape_gamma,
      rate = rate_gamma
    )
  )
)

# AIC.
AIC_gamma <- -2 * loglik_gamma + 2 * p

# BIC.
BIC_gamma <- -2 * loglik_gamma + log(n) * p

# Tabla resumen.
met_gamma <- data.frame(
  modelo = "Gamma",
  loglik = loglik_gamma,
  AIC    = AIC_gamma,
  BIC    = BIC_gamma
)

# Mostrar métricas.
cat("\n== Métricas de ajuste GAMMA ==\n")
print(met_gamma)


# VaR/TVaR 99.5% por simulación

# Fijar semilla.
set.seed(12345)

# Simular distribución Gamma reescalada.

sim_gamma_scaled <- rgamma(2e5,shape = shape_gamma,rate = rate_gamma)

# Regresar a escala original.
sim_gamma <- sim_gamma_scaled * 1e9

# Nivel de confianza.
alpha <- 0.995

# VaR 99.5%.
VaR_gamma_995 <- as.numeric(quantile(sim_gamma, alpha))

# TVaR 99.5%.
TVaR_gamma_995 <- mean(sim_gamma[sim_gamma > VaR_gamma_995])

# Mostrar resultados.

cat("\n== VaR/TVaR 99.5% (GAMMA, simulación) ==\n")

print(
  data.frame(
    modelo = "Gamma",
    VaR_995 = VaR_gamma_995,
    TVaR_995 = TVaR_gamma_995
  )
)


# Exportar resumen a CSV

resumen_gamma <- data.frame(
  shape = shape_gamma,
  rate  = rate_gamma,
  
  loglik = loglik_gamma,
  AIC    = AIC_gamma,
  BIC    = BIC_gamma,
  
  KS_p = KS$p.value,
  AD_p = AD$p.value,
  
  VaR_99 = qgamma(
    0.99,
    shape = shape_gamma,
    rate = rate_gamma
  ) * 1e9,
  
  VaR_995 = qgamma(
    0.995,
    shape = shape_gamma,
    rate = rate_gamma
  ) * 1e9,
  
  TVaR_995 = TVaR_gamma_995
)

# Exportar archivo CSV.

write.csv(
  resumen_gamma,
  "resumen_gamma_metrics.csv",
  row.names = FALSE
)

cat("\n>> Guardado 'resumen_gamma_metrics.csv'\n")

