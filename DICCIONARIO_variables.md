# Diccionario de variables

Descripción de las columnas de las bases agregadas de `Datos derivados/`. Todas las tasas usan nacidos vivos del DANE como denominador. Las muertes maternas tempranas son las defunciones con causa básica en los códigos CIE-10 O00 a O95, O98 a O99 o A34 (muerte durante el embarazo, el parto o dentro de los 42 días).

## `Serie anual numeradores y denominadores.csv`
- `anio` — año.
- `d_nac33` — muertes maternas tempranas en los 33 territorios del país.
- `d_cor` — muertes maternas tempranas en Córdoba.
- `b_nac33` — nacidos vivos en los 33 territorios.
- `b_cor` — nacidos vivos en Córdoba.
- `b_nacall` — nacidos vivos del país, todas las residencias.

## `Serie anual con intervalos.csv`
- `anio` — año.
- `d_cor`, `b_cor` — muertes tempranas y nacidos vivos en Córdoba.
- `rmm_cor`, `rmm_cor_lo`, `rmm_cor_hi` — razón de mortalidad materna de Córdoba con su intervalo exacto de Poisson (x100.000 NV).
- `d_nac`, `b_nac`, `rmm_nac`, `rmm_nac_lo`, `rmm_nac_hi` — lo mismo a nivel nacional.
- `mme_cor`, `mmer_cor`, `mmer_cor_lo`, `mmer_cor_hi` — morbilidad materna extrema y su razón (x1.000 NV), Córdoba, con intervalo.
- `mme_nac`, `mmer_nac`, `mmer_nac_lo`, `mmer_nac_hi` — lo mismo a nivel nacional.
- `d550_cor`, `d550_nac` — muertes notificadas (evento 550) en Córdoba y en el país.
- `idx_cor`, `idx_cor_lo`, `idx_cor_hi`, `idx_nac`, `idx_nac_lo`, `idx_nac_hi` — índice de mortalidad basado en notificaciones (%), con intervalo.

## `Muertes tempranas y tardias por anio.csv`
- `anio` — año.
- `early` — muertes maternas tempranas (numerador de la razón).
- `late` — muertes maternas tardías (código O96; se excluyen de la razón y se analizan aparte).
- `b_cor` — nacidos vivos en Córdoba.
- `total` — suma de tempranas y tardías.
- `yc` — año centrado en 2015.
- `covid` — indicador del periodo pandémico (1 en 2020 y 2021).

## `Causas de las 162 muertes.csv`
- `anio` — año.
- `C_BAS1`, `cod` — causa básica (código CIE-10).
- `grupo` — grupo de causa (hipertensivos, hemorragia, aborto, infección, embolia, otras directas, indirectas, no determinada).
- `periodo` — 2015-2021 o 2022-2024.

## `Causas por grupo y periodo.csv`
- `grupo` — grupo de causa.
- `2015-2021`, `2022-2024` — número de muertes en cada periodo.
- `total` — total en la década.

## `Base municipal con suavizado y LISA.csv`
- `name`, `municipio` — nombre del municipio.
- `cod`, `cod3`, `cod_municipio` — códigos del municipio.
- `subregion` — subregión del departamento.
- `nacidos_vivos`, `muertes_tempranas` — denominador y numerador municipales de la década.
- `IPM` — índice de pobreza multidimensional municipal (2018).
- `RMM_cruda` — razón de mortalidad materna cruda (x100.000 NV).
- `RMM_EBglobal` — razón suavizada por bayes empírico global (Marshall).
- `lisa_p`, `lisa_q` — valor p e índice del conglomerado local de asociación espacial.

## `Base municipal con MME para la decision de pobreza.csv`
Igual que la anterior, con columnas adicionales: `mme` (morbilidad materna extrema municipal), `IPM10` (IPM dividido entre 10) y `razon_mme`. Base de apoyo para el análisis ecológico de pobreza.

## `Estandarizacion por edad materna.csv`
- `grupo_edad` — quinquenio de edad de la madre (10-14 a 45-49).
- `muertes_cordoba`, `nacidos_cordoba` — muertes tempranas y nacidos vivos en Córdoba por grupo de edad.
- `muertes_resto`, `nacidos_resto` — lo mismo en el resto del país.

## `Ranking departamental RMM.csv`
- `dpto` — código del departamento.
- `Departamento` — nombre.
- `muertes`, `nacidos` — muertes tempranas y nacidos vivos de la década.
- `RMM` — razón de mortalidad materna acumulada (x100.000 NV).
- `rank` — puesto en el orden descendente.

## `Notificaciones de mortalidad materna por anio.csv`
- `anio` — año.
- `notificadas` — muertes maternas notificadas a la vigilancia (evento 550), para la comparación de totales.

## `Flujo de referencia.csv`
- `flujo` — dirección del flujo (entradas o salidas del departamento).
- `departamento` — departamento de residencia.
- `muertes` — número de muertes.

## `cordoba_municipios.geojson`
Geometría de los 30 municipios de Córdoba (geoportal DANE), con el código municipal para el enlace espacial.
