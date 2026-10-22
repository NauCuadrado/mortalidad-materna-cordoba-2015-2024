# Mortalidad materna en Cordoba, 2015 a 2024
# Naudith C, agosto de 2026
#
# Los mismos analisis del pipeline, escritos otra vez con las funciones de R,
# para comprobar que las cifras no dependen del programa. Se corre desde la raiz
# del paquete:  Rscript "Codigo/Analisis en R.R"

datos <- "Datos derivados"

R  <- read.csv(file.path(datos, "Serie anual numeradores y denominadores.csv"))
C  <- read.csv(file.path(datos, "Causas de las 162 muertes.csv"))
Td <- read.csv(file.path(datos, "Muertes tempranas y tardias por anio.csv"))
M  <- read.csv(file.path(datos, "Base municipal con suavizado y LISA.csv"))

R$d_resto <- R$d_nac33 - R$d_cor
R$b_resto <- R$b_nac33 - R$b_cor
C$tipo <- ifelse(C$grupo == "Indirect causes", "Indirect", "Direct")

anios <- R$anio
nv    <- R$b_cor
z     <- qnorm(0.975)


# razones e intervalos ---------------------------------------------------

rmm <- function(k, n, m = 1e5) {
  lo <- if (k > 0) qchisq(0.025, 2 * k) / 2 else 0
  hi <- qchisq(0.975, 2 * (k + 1)) / 2
  c(k / n, lo / n, hi / n) * m
}

rr <- function(k1, n1, k2, n2) {
  r  <- (k1 / n1) / (k2 / n2)
  se <- sqrt(1 / k1 + 1 / k2)
  c(r, r * exp(-z * se), r * exp(z * se))
}

rr_exacto <- function(k1, n1, k2, n2) {
  t <- n1 / n2
  p_lo <- if (k1 > 0) qbeta(0.025, k1, k2 + 1) else 0
  p_hi <- if (k2 > 0) qbeta(0.975, k1 + 1, k2) else 1
  c((k1 / n1) / (k2 / n2), (p_lo / (1 - p_lo)) / t, (p_hi / (1 - p_hi)) / t)
}

ror <- function(a1, m1, b1, n1, a2, m2, b2, n2) {
  q  <- ((a2 / m2) / (b2 / n2)) / ((a1 / m1) / (b1 / n1))
  se <- sqrt(1 / a1 + 1 / b1 + 1 / a2 + 1 / b2)
  c(q, q * exp(-z * se), q * exp(z * se), 2 * pnorm(-abs(log(q) / se)))
}

cat("\n== 1. MAGNITUD Y COMPARACION CON EL RESTO DEL PAIS ==\n")
a <- rmm(sum(R$d_cor), sum(R$b_cor))
b <- rmm(sum(R$d_resto), sum(R$b_resto))
q <- rr(sum(R$d_cor), sum(R$b_cor), sum(R$d_resto), sum(R$b_resto))
cat(sprintf("  departamento %.1f (%.1f a %.1f) | resto del pais %.1f (%.1f a %.1f) | RR %.2f (%.2f a %.2f)\n",
            a[1], a[2], a[3], b[1], b[2], b[3], q[1], q[2], q[3]))

anual <- do.call(rbind, lapply(seq_len(nrow(R)), function(i) {
  x <- R[i, ]
  q <- rr(x$d_cor, x$b_cor, x$d_resto, x$b_resto)
  e <- rr_exacto(x$d_cor, x$b_cor, x$d_resto, x$b_resto)
  data.frame(anio = x$anio, rr = q[1], lo = q[2], hi = q[3],
             lo_exacto = e[2], hi_exacto = e[3])
}))
cat("  anios distinguibles:",
    paste(anual$anio[anual$lo > 1 | anual$hi < 1], collapse = ", "), "(log-normal) |",
    paste(anual$anio[anual$lo_exacto > 1 | anual$hi_exacto < 1], collapse = ", "),
    "(exacto condicional)\n")

per <- list("2015-2021" = R[R$anio <= 2021, ], "2022-2024" = R[R$anio >= 2022, ])
for (k in names(per)) {
  s <- per[[k]]
  t <- rmm(sum(s$d_cor), sum(s$b_cor))
  q <- rr(sum(s$d_cor), sum(s$b_cor), sum(s$d_resto), sum(s$b_resto))
  cat(sprintf("  %s  RMM %.1f (%.1f a %.1f)  RR %.2f (%.2f a %.2f)\n",
              k, t[1], t[2], t[3], q[1], q[2], q[3]))
}
g <- with(per, ror(sum(`2015-2021`$d_cor), sum(`2015-2021`$b_cor),
                   sum(`2015-2021`$d_resto), sum(`2015-2021`$b_resto),
                   sum(`2022-2024`$d_cor), sum(`2022-2024`$b_cor),
                   sum(`2022-2024`$d_resto), sum(`2022-2024`$b_resto)))
cat(sprintf("  razon de razones 2022-2024 frente a 2015-2021: %.2f (%.2f a %.2f) p=%.3f\n",
            g[1], g[2], g[3], g[4]))


# tendencias -------------------------------------------------------------

series <- list(
  "Todas las tempranas" = R$d_cor,
  "Causas directas"     = sapply(anios, function(y) sum(C$tipo == "Direct"   & C$anio == y)),
  "Causas indirectas"   = sapply(anios, function(y) sum(C$tipo == "Indirect" & C$anio == y)),
  "Muertes tardias"     = Td$late)

tendencia <- function(y, off, yy) {
  d <- data.frame(y = y, yc = yy - min(yy), covid = as.integer(yy %in% c(2020, 2021)))
  m <- glm(y ~ yc + covid, family = poisson, offset = log(off), data = d)
  s <- summary(m)$coefficients
  disp <- sum(residuals(m, "pearson")^2) / df.residual(m)
  wald <- function(i, sc = 1) {
    e <- s[i, 1]; se <- s[i, 2] * sqrt(sc)
    c(e, e - z * se, e + z * se, 2 * pnorm(-abs(e / se)))
  }
  list(apc = wald("yc"), covid = wald("covid"), disp = disp,
       apc_q = wald("yc", disp), covid_q = wald("covid", disp))
}

pct <- function(v) (exp(v) - 1) * 100

cat("\n== 3. TENDENCIAS ANUALES ==\n")
for (k in names(series)) {
  r <- tendencia(series[[k]], nv, anios)
  cat(sprintf("  %-20s APC %+.1f%% (%+.1f a %+.1f) p=%.3f | pandemia RR %.2f (%.2f-%.2f) p=%.3f | dispersion %.2f\n",
              k, pct(r$apc[1]), pct(r$apc[2]), pct(r$apc[3]), r$apc[4],
              exp(r$covid[1]), exp(r$covid[2]), exp(r$covid[3]), r$covid[4], r$disp))
  cat(sprintf("  %-20s quasi-Poisson APC %+.1f%% (%+.1f a %+.1f)\n",
              "", pct(r$apc_q[1]), pct(r$apc_q[2]), pct(r$apc_q[3])))
  p <- sapply(anios, function(o) {
    i <- anios != o
    tendencia(series[[k]][i], nv[i], anios[i])$apc[4]
  })
  cat(sprintf("  %-20s pierde significacion en %d de %d omisiones\n", "", sum(p >= 0.05), length(p)))
}


# modelo segmentado ------------------------------------------------------

cat("\n== 4. MODELO SEGMENTADO CON NODOS EN 2019 Y 2021 ==\n")
seg <- data.frame(y  = R$d_cor,
                  x0 = anios - 2015,
                  x1 = pmax(anios - 2019, 0),
                  x2 = pmax(anios - 2021, 0))
ms <- glm(y ~ x0 + x1 + x2, family = poisson, offset = log(nv), data = seg)
V  <- vcov(ms)
tramos <- list("Pre-pandemia 2015-2019"  = list(c(0, 1, 0, 0), anios <= 2019),
               "Pandemia 2020-2021"      = list(c(0, 1, 1, 0), anios %in% c(2020, 2021)),
               "Pos-pandemia 2022-2024"  = list(c(0, 1, 1, 1), anios >= 2022))
for (k in names(tramos)) {
  cc <- tramos[[k]][[1]]
  e  <- sum(cc * coef(ms))
  se <- sqrt(as.numeric(t(cc) %*% V %*% cc))
  cat(sprintf("  %-24s n=%3d  APC %+.1f%% (%+.1f a %+.1f) p=%.3f\n",
              k, sum(R$d_cor[tramos[[k]][[2]]]), pct(e), pct(e - z * se), pct(e + z * se),
              2 * pnorm(-abs(e / se))))
}


# modelos conjuntos ------------------------------------------------------

conjunto <- function(y1, y2, etiqueta) {
  d <- data.frame(y = c(y1, y2),
                  g = rep(0:1, each = length(anios)),
                  yc = rep(anios - 2015, 2),
                  covid = rep(as.integer(anios %in% c(2020, 2021)), 2),
                  off = rep(nv, 2))
  m <- glm(y ~ g + yc + covid + g:yc + g:covid, family = poisson,
           offset = log(off), data = d)
  s <- summary(m)$coefficients
  disp <- sum(residuals(m, "pearson")^2) / df.residual(m)
  for (i in c("g:yc", "g:covid")) {
    e <- s[i, 1]; se <- s[i, 2]
    cat(sprintf("  %-8s x %-9s RR %.2f (%.2f a %.2f) p=%.4f | quasi-Poisson (%.2f a %.2f)\n",
                etiqueta, c("g:yc" = "anio", "g:covid" = "pandemia")[i],
                exp(e), exp(e - z * se), exp(e + z * se), 2 * pnorm(-abs(e / se)),
                exp(e - z * se * sqrt(disp)), exp(e + z * se * sqrt(disp))))
  }
  invisible(m)
}

cat("\n== 5. MODELOS CONJUNTOS CON INTERACCION ==\n")
conjunto(series[["Causas directas"]], series[["Causas indirectas"]], "causa")
m2 <- conjunto(series[["Todas las tempranas"]], series[["Muertes tardias"]], "ventana")
b <- coef(m2)
cat(sprintf("     pendiente tempranas %+.1f%% | tardias %+.1f%%\n",
            pct(b["yc"]), pct(b["yc"] + b["g:yc"])))
cat(sprintf("     pandemia tempranas RR %.2f | tardias RR %.2f\n",
            exp(b["covid"]), exp(b["covid"] + b["g:covid"])))


# causas por periodo -----------------------------------------------------

cat("\n== 6. CAUSAS ==\n")
n1 <- sum(R$b_cor[R$anio <= 2021])
n2 <- sum(R$b_cor[R$anio >= 2022])
for (g in sort(unique(C$grupo))) {
  k1 <- sum(C$grupo == g & C$anio <= 2021)
  k2 <- sum(C$grupo == g & C$anio >= 2022)
  a <- rmm(k1, n1); b <- rmm(k2, n2)
  cat(sprintf("  %-30s n=%3d  %.1f (%.1f-%.1f) -> %.1f (%.1f-%.1f)\n",
              g, k1 + k2, a[1], a[2], a[3], b[1], b[2], b[3]))
}


# espacial y ecologico ---------------------------------------------------

cat("\n== 7. ESPACIAL Y ECOLOGICO ==\n")
if (requireNamespace("sf", quietly = TRUE) && requireNamespace("spdep", quietly = TRUE)) {
  sf::sf_use_s2(FALSE)
  gj <- sf::st_read(file.path(datos, "cordoba_municipios.geojson"), quiet = TRUE)
  gj$cod <- as.character(gj$cod)
  gdf <- merge(gj, transform(M, cod_municipio = as.character(cod_municipio)),
               by.x = "cod", by.y = "cod_municipio")
  gdf <- gdf[order(gdf$cod), ]
  w <- spdep::nb2listw(spdep::poly2nb(gdf, queen = TRUE), style = "W")
  set.seed(42)
  for (col in c("RMM_cruda", "RMM_EBglobal")) {
    mi <- spdep::moran.mc(gdf[[col]], w, nsim = 999)
    cat(sprintf("  %-14s I=%+.3f E=%+.3f p=%.3f\n", col, mi$statistic,
                -1 / (nrow(gdf) - 1), mi$p.value))
  }
  mr <- spdep::EBImoran.mc(gdf$muertes_tempranas, gdf$nacidos_vivos, w, nsim = 999)
  cat(sprintf("  %-14s I=%+.3f E=%+.3f p=%.3f\n", "EB ajustado", mr$statistic,
              -1 / (nrow(gdf) - 1), mr$p.value))
  li <- spdep::localmoran_perm(gdf$RMM_cruda, w, nsim = 999)
  cat(sprintf("  LISA sobre la tasa cruda: municipios con p<0.05: %d de %d\n",
              sum(li[, "Pr(z != E(Ii))"] < 0.05), nrow(gdf)))
} else {
  cat("  sf o spdep no instalados; se omite la prueba espacial\n")
}

ipm10 <- M$IPM / 10
me    <- glm(muertes_tempranas ~ I(IPM / 10), family = poisson,
             offset = log(nacidos_vivos), data = M)
s <- summary(me)$coefficients[2, ]
cat(sprintf("  ecologico IRR %.2f (%.2f a %.2f) p=%.3f dispersion %.2f\n",
            exp(s[1]), exp(s[1] - z * s[2]), exp(s[1] + z * s[2]), 2 * pnorm(-abs(s[1] / s[2])),
            sum(residuals(me, "pearson")^2) / df.residual(me)))

sin_capital <- M[M$cod_municipio != 23001, ]
s <- summary(glm(muertes_tempranas ~ I(IPM / 10), family = poisson,
                 offset = log(nacidos_vivos), data = sin_capital))$coefficients[2, ]
cat(sprintf("  sin la capital IRR %.2f (%.2f a %.2f) p=%.3f\n",
            exp(s[1]), exp(s[1] - z * s[2]), exp(s[1] + z * s[2]), 2 * pnorm(-abs(s[1] / s[2]))))

obs <- sum(M$muertes_tempranas)
for (irr in c(1.05, 1.10, 1.15, 1.20)) {
  set.seed(7)
  lam <- M$nacidos_vivos * exp(log(irr) * (ipm10 - mean(ipm10)))
  lam <- lam * obs / sum(lam)
  h <- sum(replicate(3000, {
    y <- rpois(length(lam), lam)
    p <- summary(glm(y ~ ipm10, family = poisson,
                     offset = log(M$nacidos_vivos)))$coefficients[2, 4]
    p < 0.05
  }))
  cat(sprintf("  potencia IRR %.2f -> %.1f%% (muertes esperadas conservadas: %.1f)\n",
              irr, h / 3000 * 100, sum(lam)))
}
