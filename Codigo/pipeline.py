# -*- coding: utf-8 -*-
"""
Pipeline reproducible del articulo de mortalidad materna, Cordoba 2015 a 2024.

QUE HACE ESTE ARCHIVO
    Recalcula, a partir de las bases derivadas incluidas en la carpeta "Datos derivados",
    las cifras estadisticas que aparecen en el manuscrito y en el material suplementario,
    y las escribe en "Salidas/resultados.txt" junto con los CSV de cada bloque.
    Concretamente reproduce:

      1. Razones de mortalidad del departamento y del resto del pais, con intervalos
         exactos de Poisson, y razones de tasas anuales y por periodo (Tabla 4, Tabla S1).
      2. Razon de razones de tasas entre 2015-2021 y 2022-2024 (Tabla 4).
      3. Posicion del departamento en el orden de los 33 territorios (Figura S1).
      4. Modelos de tendencia de Poisson para las cuatro series, con indicador de pandemia,
         reajuste omitiendo un anio cada vez y escala quasi-Poisson (Tablas S5 y S7).
      5. Modelo segmentado con nodos en 2019 y 2021 (Tabla S2).
      6. Modelos conjuntos con terminos de interaccion: causa por anio, causa por pandemia,
         ventana de defuncion por anio y ventana por pandemia (Resultados y Tabla S6).
      7. Descomposicion por causa con tasas e intervalos por periodo (Tabla 3).
      8. Moran global y LISA sobre las tasas municipales, modelo ecologico sobre el indice
         de pobreza y potencia por simulacion (Tablas S3 y S4).
      9. Flujo de ocurrencia de las muertes dentro y fuera del departamento, y dependencia
         de la capital por subregion de residencia (Figura 6 y Tabla S10).
     10. Perfil de las 162 mujeres por variable sociodemografica, con no informados (Tabla 1).
     11. Conteo anual de las notificaciones de mortalidad materna al sistema de vigilancia
         y su comparacion agregada con el registro civil (Tabla 2).
     12. Deteccion de registros identicos en todas las variables disponibles.

QUE NO HACE
    No reconstruye las bases derivadas desde los microdatos anuales del DANE y del SIVIGILA.
    Ese paso esta documentado, con las funciones de lectura y depuracion, en
    "Codigo/Codigo reproducible mortalidad.py"; requiere descargar los archivos originales,
    que no se redistribuyen en este paquete.
    No genera las figuras. Las figuras se producen con "Codigo/Figuras del articulo.py",
    que lee las mismas bases derivadas. El suplemento se arma con
    "Codigo/Generador del suplemento.py".

USO
    python pipeline.py

REQUISITOS
    Ver requirements.txt. Los bloques espacial y ecologico requieren geopandas, libpysal y
    esda; si no estan instalados el bloque se omite y queda constancia en el registro.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm, beta
import statsmodels.api as sm

BASE = os.path.dirname(os.path.abspath(__file__))
if not os.path.isdir(os.path.join(BASE, 'Datos derivados')):
    # permite ejecutar el archivo desde una subcarpeta como "Codigo"
    BASE = os.path.dirname(BASE)
DATA = os.path.join(BASE, 'Datos derivados')
OUT = os.path.join(BASE, 'Salidas')
os.makedirs(OUT, exist_ok=True)
LOG = []

SEED_SPATIAL, SEED_POWER = 42, 7


def log(x):
    print(x)
    LOG.append(str(x))


# ------------------------------------------------------------------ utilidades
def ic(k, n, m=1e5):
    """Intervalo exacto de Poisson para una tasa."""
    lo = chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0
    hi = chi2.ppf(0.975, 2 * (k + 1)) / 2
    return k / n * m, lo / n * m, hi / n * m


def rr(k1, n1, k2, n2):
    """Razon de tasas con intervalo log-normal."""
    r = (k1 / n1) / (k2 / n2)
    se = np.sqrt(1 / k1 + 1 / k2)
    return r, r * np.exp(-1.96 * se), r * np.exp(1.96 * se)


def rr_exacto(k1, n1, k2, n2, conf=0.95):
    """Intervalo exacto condicional para la razon de tasas, sin aproximacion normal."""
    a = (1 - conf) / 2
    t = n1 / n2
    lo_p = beta.ppf(a, k1, k2 + 1) if k1 > 0 else 0.0
    hi_p = beta.ppf(1 - a, k1 + 1, k2) if k2 > 0 else 1.0
    f = lambda p: (p / (1 - p)) / t if p < 1 else float('inf')
    return (k1 / n1) / (k2 / n2), f(lo_p), f(hi_p)


def ror(k1a, n1a, k1b, n1b, k2a, n2a, k2b, n2b):
    """Razon entre dos razones de tasas, con intervalo log-normal."""
    r1 = (k1a / n1a) / (k1b / n1b)
    r2 = (k2a / n2a) / (k2b / n2b)
    q = r2 / r1
    se = np.sqrt(1 / k1a + 1 / k1b + 1 / k2a + 1 / k2b)
    z = np.log(q) / se
    return q, q * np.exp(-1.96 * se), q * np.exp(1.96 * se), 2 * (1 - norm.cdf(abs(z)))


def poisson_trend(y, off, anios, scale=None):
    df = pd.DataFrame({'y': list(y), 'anio': list(anios)})
    df['yc'] = df.anio - df.anio.min()
    df['covid'] = df.anio.isin([2020, 2021]).astype(int)
    glm = sm.GLM(df.y, sm.add_constant(df[['yc', 'covid']]),
               family=sm.families.Poisson(), offset=np.log(off))
    mp = glm.fit()
    disp = mp.pearson_chi2 / mp.df_resid
    if scale == 'quasi' and disp > 1:
        m = glm.fit(scale='X2', use_t=True)
    else:
        m = mp
    apc = (np.exp(m.params['yc']) - 1) * 100
    ci = (np.exp(m.conf_int().loc['yc']) - 1) * 100
    covid = np.exp(m.params['covid'])
    cci = np.exp(m.conf_int().loc['covid'])
    return dict(apc=apc, apc_lo=ci[0], apc_hi=ci[1], p=m.pvalues['yc'],
                covid=covid, covid_lo=cci[0], covid_hi=cci[1], p_covid=m.pvalues['covid'],
                dispersion=m.pearson_chi2 / m.df_resid)


def omit_one(y, off, anios):
    out = []
    for k in anios:
        sel = [i for i, a in enumerate(anios) if a != k]
        r = poisson_trend([y[i] for i in sel], [off[i] for i in sel], [anios[i] for i in sel])
        out.append(dict(omitido=k, apc=r['apc'], p=r['p']))
    return pd.DataFrame(out)


# ------------------------------------------------------------------ carga
def cargar():
    R = pd.read_csv(os.path.join(DATA, 'Serie anual numeradores y denominadores.csv'))
    R['d_resto'] = R.d_nac33 - R.d_cor
    R['b_resto'] = R.b_nac33 - R.b_cor
    C = pd.read_csv(os.path.join(DATA, 'Causas de las 162 muertes.csv'))
    C['tipo'] = np.where(C.grupo == 'Indirect causes', 'Indirect', np.where(C.grupo == 'Undetermined causes', 'Undetermined', 'Direct'))
    Td = pd.read_csv(os.path.join(DATA, 'Muertes tempranas y tardias por anio.csv'))
    M = pd.read_csv(os.path.join(DATA, 'Base municipal con suavizado y LISA.csv'))
    return R, C, Td, M


# ------------------------------------------------------------------ 1. magnitud
def bloque_magnitud(R):
    log('\n== 1. MAGNITUD Y COMPARACION CON EL RESTO DEL PAIS ==')
    a = ic(R.d_cor.sum(), R.b_cor.sum())
    b = ic(R.d_resto.sum(), R.b_resto.sum())
    q = rr(R.d_cor.sum(), R.b_cor.sum(), R.d_resto.sum(), R.b_resto.sum())
    log('  departamento %.1f (%.1f a %.1f) | resto del pais %.1f (%.1f a %.1f) | RR %.2f (%.2f a %.2f)'
        % (a + b + q))
    filas = []
    for _, r0 in R.iterrows():
        q = rr(r0.d_cor, r0.b_cor, r0.d_resto, r0.b_resto)
        e = rr_exacto(int(r0.d_cor), int(r0.b_cor), int(r0.d_resto), int(r0.b_resto))
        filas.append(dict(anio=int(r0.anio), rr=q[0], lo=q[1], hi=q[2],
                          lo_exacto=e[1], hi_exacto=e[2],
                          distinguible=bool(q[1] > 1 or q[2] < 1),
                          distinguible_exacto=bool(e[1] > 1 or e[2] < 1)))
    F = pd.DataFrame(filas)
    F.to_csv(os.path.join(OUT, 'RR anual contra el resto del pais.csv'), index=False)
    log('  anios distinguibles del resto del pais: %s (log-normal) | %s (exacto condicional)'
        % (F[F.distinguible].anio.tolist(), F[F.distinguible_exacto].anio.tolist()))
    P = {}
    for lab, yy in [('2015-2021', range(2015, 2022)), ('2022-2024', range(2022, 2025))]:
        s = R[R.anio.isin(list(yy))]
        P[lab] = (s.d_cor.sum(), s.b_cor.sum(), s.d_resto.sum(), s.b_resto.sum())
        q = rr(*P[lab])
        t = ic(s.d_cor.sum(), s.b_cor.sum())
        log('  %s  RMM %.1f (%.1f a %.1f)  RR %.2f (%.2f a %.2f)' % ((lab,) + t + q))
    g = ror(*P['2015-2021'], *P['2022-2024'])
    log('  razon de razones 2022-2024 frente a 2015-2021: %.2f (%.2f a %.2f) p=%.3f' % g)
    log('  lectura: el cambio del exceso relativo entre periodos no es distinguible del azar')
    return F


# ------------------------------------------------------------------ 2. ranking
def bloque_ranking():
    f = os.path.join(DATA, 'Ranking departamental RMM.csv')
    if not os.path.exists(f):
        log('\n== 2. ORDEN DE LOS TERRITORIOS: base no disponible ==')
        return
    log('\n== 2. ORDEN DE LOS 33 TERRITORIOS POR RMM ACUMULADA ==')
    K = pd.read_csv(f).sort_values('rank')
    fila = K[K.Departamento.str.contains('rdoba', case=False, na=False)]
    if len(fila):
        r = fila.iloc[0]
        log('  Cordoba: RMM %.1f, posicion %d de %d' % (r.RMM, int(r['rank']), len(K)))
    log('  cinco primeros: %s' % ', '.join('%s %.1f' % (x.Departamento, x.RMM)
                                           for x in K.head(5).itertuples()))
    log('  mediana de los 33 territorios: %.1f' % K.RMM.median())


# ------------------------------------------------------------------ 3. tendencias
def series_anuales(R, C, Td):
    yrs = R.anio.tolist()
    return yrs, R.b_cor.tolist(), {
        'Todas las tempranas': R.d_cor.tolist(),
        'Causas directas': [int(((C.tipo == 'Direct') & (C.anio == y)).sum()) for y in yrs],
        'Causas indirectas': [int(((C.tipo == 'Indirect') & (C.anio == y)).sum()) for y in yrs],
        'Muertes tardias': Td.late.tolist()}


def bloque_tendencias(R, C, Td):
    log('\n== 3. TENDENCIAS ANUALES ==')
    yrs, off, series = series_anuales(R, C, Td)
    filas = []
    for lab, y in series.items():
        r = poisson_trend(y, off, yrs)
        q = poisson_trend(y, off, yrs, scale='quasi')
        log('  %-20s APC %+.1f%% (%+.1f a %+.1f) p=%.3f | pandemia RR %.2f (%.2f-%.2f) p=%.3f | dispersion %.2f'
            % (lab, r['apc'], r['apc_lo'], r['apc_hi'], r['p'],
               r['covid'], r['covid_lo'], r['covid_hi'], r['p_covid'], r['dispersion']))
        log('  %-20s quasi-Poisson APC %+.1f%% (%+.1f a %+.1f)'
            % ('', q['apc'], q['apc_lo'], q['apc_hi']))
        filas.append(dict(serie=lab, **r, apc_quasi_lo=q['apc_lo'], apc_quasi_hi=q['apc_hi'],
                          p_apc_quasi=q['p'], pandemia_quasi_lo=q['covid_lo'],
                          pandemia_quasi_hi=q['covid_hi'], p_pandemia_quasi=q['p_covid']))
        oo = omit_one(y, off, yrs)
        oo['serie'] = lab
        oo.to_csv(os.path.join(OUT, 'Omitir un anio - %s.csv' % lab.replace(' ', '_')), index=False)
        log('  %-20s pierde significacion en %d de %d omisiones'
            % ('', int((oo.p >= 0.05).sum()), len(oo)))
    pd.DataFrame(filas).to_csv(os.path.join(OUT, 'Modelos de tendencia.csv'), index=False)


def bloque_segmentado(R):
    """Modelo segmentado con cambios de pendiente en 2019 y 2021 (Tabla S2)."""
    log('\n== 4. MODELO SEGMENTADO CON NODOS EN 2019 Y 2021 ==')
    yrs = np.array(R.anio.tolist())
    y = np.array(R.d_cor.tolist())
    off = np.log(R.b_cor.values)
    x0 = yrs - 2015
    x1 = np.clip(yrs - 2019, 0, None)
    x2 = np.clip(yrs - 2021, 0, None)
    X = sm.add_constant(np.column_stack([x0, x1, x2]))
    m = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit()
    cov = np.asarray(m.cov_params())
    tramos = [('Pre-pandemia 2015-2019', np.array([0, 1, 0, 0]), int(y[(yrs >= 2015) & (yrs <= 2019)].sum())),
              ('Pandemia 2020-2021', np.array([0, 1, 1, 0]), int(y[(yrs >= 2020) & (yrs <= 2021)].sum())),
              ('Pos-pandemia 2022-2024', np.array([0, 1, 1, 1]), int(y[yrs >= 2022].sum()))]
    filas = []
    for lab, c, n in tramos:
        coef = float(c @ m.params)
        se = float(np.sqrt(c @ cov @ c))
        apc = (np.exp(coef) - 1) * 100
        lo = (np.exp(coef - 1.96 * se) - 1) * 100
        hi = (np.exp(coef + 1.96 * se) - 1) * 100
        p = 2 * (1 - norm.cdf(abs(coef / se)))
        log('  %-24s n=%3d  APC %+.1f%% (%+.1f a %+.1f) p=%.3f' % (lab, n, apc, lo, hi, p))
        filas.append(dict(tramo=lab, muertes=n, apc=apc, lo=lo, hi=hi, p=p))
    pd.DataFrame(filas).to_csv(os.path.join(OUT, 'Modelo segmentado.csv'), index=False)


def bloque_interacciones(R, C, Td):
    """Modelos conjuntos: prueban si dos series se distinguen entre si."""
    log('\n== 5. MODELOS CONJUNTOS CON INTERACCION ==')
    yrs, off, series = series_anuales(R, C, Td)
    yc = np.array(yrs) - 2015
    covid = np.isin(yrs, [2020, 2021]).astype(int)
    YC = np.r_[yc, yc]
    CV = np.r_[covid, covid]
    OF = np.log(np.r_[off, off])

    # a) causas directas frente a indirectas
    y = np.array(series['Causas directas'] + series['Causas indirectas'])
    grupo = np.r_[np.zeros(len(yrs)), np.ones(len(yrs))]
    X = sm.add_constant(np.column_stack([grupo, YC, CV, grupo * YC, grupo * CV]))
    m = sm.GLM(y, X, family=sm.families.Poisson(), offset=OF).fit()
    _d = m.pearson_chi2 / m.df_resid
    mq = sm.GLM(y, X, family=sm.families.Poisson(), offset=OF).fit(scale='X2', use_t=True) if _d > 1 else m
    ci = m.conf_int(); ciq = mq.conf_int()
    for k, lab in [(4, 'causa x anio'), (5, 'causa x pandemia')]:
        log('  %-18s RR %.2f (%.2f a %.2f) p=%.4f | quasi-Poisson (%.2f a %.2f)'
            % (lab, np.exp(m.params[k]), np.exp(ci[k][0]), np.exp(ci[k][1]), m.pvalues[k],
               np.exp(ciq[k][0]), np.exp(ciq[k][1])))

    # b) muertes tempranas frente a tardias
    #    el modelo lleva las dos interacciones: si solo se incluyera ventana x anio, las dos
    #    ventanas quedarian obligadas a compartir el mismo efecto del periodo pandemico
    y2 = np.array(series['Todas las tempranas'] + series['Muertes tardias'])
    vent = np.r_[np.zeros(len(yrs)), np.ones(len(yrs))]
    X2 = sm.add_constant(np.column_stack([vent, YC, CV, vent * YC, vent * CV]))
    m2 = sm.GLM(y2, X2, family=sm.families.Poisson(), offset=OF).fit()
    _d2 = m2.pearson_chi2 / m2.df_resid
    m2q = sm.GLM(y2, X2, family=sm.families.Poisson(), offset=OF).fit(scale='X2', use_t=True) if _d2 > 1 else m2
    ci2 = m2.conf_int(); ci2q = m2q.conf_int()
    log('  ventana x anio     %+.1f%% por anio (%+.1f a %+.1f) p=%.4f'
        % ((np.exp(m2.params[4]) - 1) * 100, (np.exp(ci2[4][0]) - 1) * 100,
           (np.exp(ci2[4][1]) - 1) * 100, m2.pvalues[4]))
    log('  ventana x pandemia RR %.2f (%.2f a %.2f) p=%.4f | quasi-Poisson (%.2f a %.2f)'
        % (np.exp(m2.params[5]), np.exp(ci2[5][0]), np.exp(ci2[5][1]), m2.pvalues[5],
           np.exp(ci2q[5][0]), np.exp(ci2q[5][1])))
    log('     pendiente tempranas %+.1f%% | tardias %+.1f%%'
        % ((np.exp(m2.params[2]) - 1) * 100, (np.exp(m2.params[2] + m2.params[4]) - 1) * 100))
    log('     pandemia tempranas RR %.2f | tardias RR %.2f'
        % (np.exp(m2.params[3]), np.exp(m2.params[3] + m2.params[5])))


# ------------------------------------------------------------------ 6. causas
def bloque_causas(R, C):
    log('\n== 6. CAUSAS ==')
    n1 = int(R[R.anio <= 2021].b_cor.sum())
    n2 = int(R[R.anio >= 2022].b_cor.sum())
    filas = []
    for g in sorted(C.grupo.unique()):
        k = int((C.grupo == g).sum())
        k1 = int(((C.grupo == g) & (C.anio <= 2021)).sum())
        k2 = int(((C.grupo == g) & (C.anio >= 2022)).sum())
        a = ic(k1, n1)
        b = ic(k2, n2)
        filas.append(dict(grupo=g, n=k, pct=k / 162 * 100,
                          tasa1=a[0], lo1=a[1], hi1=a[2], tasa2=b[0], lo2=b[1], hi2=b[2]))
        log('  %-30s n=%3d  %.1f (%.1f-%.1f) -> %.1f (%.1f-%.1f)'
            % (g, k, a[0], a[1], a[2], b[0], b[1], b[2]))
    pd.DataFrame(filas).to_csv(os.path.join(OUT, 'Causas con IC por periodo.csv'), index=False)
    sube = [f['grupo'] for f in filas if f['grupo'] != 'Indirect causes' and f['tasa2'] > f['tasa1']]
    log('  grupos directos con tasa mayor en el segundo periodo: %s' % sube)


# ------------------------------------------------------------------ 7. espacial
def bloque_espacial(M):
    log('\n== 7. ESPACIAL Y ECOLOGICO ==')
    try:
        import geopandas as gpd
        import libpysal
        import esda
        gj = os.path.join(DATA, 'cordoba_municipios.geojson')
        if os.path.exists(gj):
            g = gpd.read_file(gj)
            g['cod'] = g['cod'].astype(str)
            M2 = M.copy()
            M2['cod_municipio'] = M2.cod_municipio.astype(str)
            gdf = (g.merge(M2, left_on='cod', right_on='cod_municipio', suffixes=('_geo', ''))
                   .sort_values('cod_municipio').reset_index(drop=True))
            w = libpysal.weights.Queen.from_dataframe(gdf, use_index=False)
            w.transform = 'r'
            for col in ['RMM_cruda', 'RMM_EBglobal']:
                np.random.seed(SEED_SPATIAL)
                mi = esda.Moran(gdf[col].values, w, permutations=999)
                log('  %-14s I=%+.3f E=%+.3f p=%.3f p95=%+.3f'
                    % (col, mi.I, mi.EI, mi.p_sim, np.percentile(mi.sim, 95)))
            np.random.seed(SEED_SPATIAL)
            mr = esda.Moran_Rate(gdf['muertes_tempranas'].values, gdf['nacidos_vivos'].values, w,
                                 adjusted=True, permutations=999)
            log('  %-14s I=%+.3f E=%+.3f p=%.3f' % ('EB ajustado', mr.I, mr.EI, mr.p_sim))
            np.random.seed(SEED_SPATIAL)
            li = esda.Moran_Local(gdf['RMM_cruda'].values, w, permutations=999)
            log('  LISA sobre la tasa cruda: municipios con p<0.05: %d de %d'
                % (int((li.p_sim < 0.05).sum()), len(gdf)))
        else:
            log('  geometria no disponible; se omite la prueba espacial')
    except ImportError:
        log('  geopandas, libpysal o esda no instalados; se omite la prueba espacial')
    d = M.muertes_tempranas.values
    nv = M.nacidos_vivos.values
    ipm10 = (M.IPM / 10).values
    X = sm.add_constant(pd.Series(ipm10, name='ipm10'))
    off = np.log(nv)
    mod = sm.GLM(d, X, offset=off, family=sm.families.Poisson()).fit()
    log('  ecologico IRR %.2f (%.2f a %.2f) p=%.3f dispersion %.2f'
        % (np.exp(mod.params['ipm10']), np.exp(mod.conf_int().loc['ipm10'])[0],
           np.exp(mod.conf_int().loc['ipm10'])[1], mod.pvalues['ipm10'],
           mod.pearson_chi2 / mod.df_resid))
    sel = (M.cod_municipio.astype(str) != '23001').values
    mod2 = sm.GLM(d[sel], sm.add_constant(pd.Series(ipm10[sel], name='ipm10')),
                  offset=np.log(nv[sel]), family=sm.families.Poisson()).fit()
    log('  sin la capital IRR %.2f (%.2f a %.2f) p=%.3f'
        % (np.exp(mod2.params['ipm10']), np.exp(mod2.conf_int().loc['ipm10'])[0],
           np.exp(mod2.conf_int().loc['ipm10'])[1], mod2.pvalues['ipm10']))
    obs = d.sum()
    filas = []
    for irr in [1.05, 1.10, 1.15, 1.20]:
        rng = np.random.default_rng(SEED_POWER)   # semilla por escenario
        lam = nv * np.exp(np.log(irr) * (ipm10 - ipm10.mean()))
        lam = lam * obs / lam.sum()
        h = 0
        for _ in range(3000):
            y = rng.poisson(lam)
            try:
                if sm.GLM(y, X, offset=off, family=sm.families.Poisson()).fit().pvalues['ipm10'] < 0.05:
                    h += 1
            except Exception:
                pass
        filas.append(dict(irr=irr, potencia=h / 3000, esperadas=lam.sum()))
        log('  potencia IRR %.2f -> %.1f%% (muertes esperadas conservadas: %.1f)'
            % (irr, h / 3000 * 100, lam.sum()))
    pd.DataFrame(filas).to_csv(os.path.join(OUT, 'Potencia del modelo ecologico.csv'), index=False)


# ------------------------------------------------------------------ 8. ocurrencia
def bloque_ocurrencia(R):
    f = os.path.join(DATA, 'Flujo de referencia.csv')
    if not os.path.exists(f):
        log('\n== 8. FLUJO DE OCURRENCIA: base no disponible ==')
        return
    log('\n== 8. LUGAR DE OCURRENCIA Y FLUJO ENTRE DEPARTAMENTOS ==')
    FL = pd.read_csv(f)
    for k, s in FL.groupby('flujo'):
        log('  %-26s %d' % (k, int(s.muertes.sum())))
        if len(s) > 1:
            log('     %s' % ', '.join('%s %d' % (x.departamento, int(x.muertes))
                                      for x in s.itertuples()))
    dentro = int(FL[FL.flujo.str.startswith('Ocurrieron')].muertes.sum())
    salen = int(FL[FL.flujo.str.startswith('Salieron')].muertes.sum())
    entran = int(FL[FL.flujo.str.startswith('Entraron')].muertes.sum())
    n = int(R.b_cor.sum())
    a = ic(dentro + salen, n)
    b = ic(dentro + entran, n)
    log('  numerador por residencia %d (RMM %.1f) | por ocurrencia %d (RMM %.1f) | diferencia %+.0f%%'
        % (dentro + salen, a[0], dentro + entran, b[0], (b[0] / a[0] - 1) * 100))



# ------------------------------------------------------------------ 8b. dependencia del nodo
def bloque_subregion():
    """Cuantas muertes de cada subregion ocurrieron en la capital departamental."""
    f = os.path.join(DATA, 'Base 162 muertes tempranas.csv')
    if not os.path.exists(f):
        log('  [omitido: requiere la base a nivel de registro, disponible por solicitud]'); return
    if not os.path.exists(f):
        return
    B = pd.read_csv(f, dtype=str)
    if 'COD_MUNIC' not in B.columns:
        log('\n== 8b. DEPENDENCIA DEL NODO: la base no trae el municipio de ocurrencia ==')
        return
    log('\n== 8b. LUGAR DE OCURRENCIA POR SUBREGION DE RESIDENCIA ==')
    B['en_capital'] = ((pd.to_numeric(B.COD_DPTO, errors='coerce') == 23) &
                       (B.COD_MUNIC.astype(str).str.zfill(3) == '001'))
    B['fuera_dpto'] = pd.to_numeric(B.COD_DPTO, errors='coerce') != 23
    T = B.groupby('subregion').agg(muertes=('cod3', 'size'), en_capital=('en_capital', 'sum'),
                                   fuera_dpto=('fuera_dpto', 'sum')).reset_index()
    T['pct_capital'] = (T.en_capital / T.muertes * 100).round(1)
    for x in T.sort_values('pct_capital', ascending=False).itertuples():
        log('  %-12s n=%3d | en la capital %3d (%.1f%%) | fuera del departamento %d'
            % (x.subregion, x.muertes, x.en_capital, x.pct_capital, x.fuera_dpto))
    fuera = B[B.cod3.astype(str).str.zfill(3) != '001']
    log('  residentes fuera de la capital: %d, de ellas %d fallecieron en la capital (%.1f%%)'
        % (len(fuera), int(fuera.en_capital.sum()), fuera.en_capital.mean() * 100))
    T.to_csv(os.path.join(OUT, 'Ocurrencia por subregion.csv'), index=False)


# ------------------------------------------------------------------ 9. perfil
def bloque_perfil():
    f = os.path.join(DATA, 'Base 162 muertes tempranas.csv')
    if not os.path.exists(f):
        log('  [omitido: requiere la base a nivel de registro, disponible por solicitud]'); return
    if not os.path.exists(f):
        log('\n== 9. PERFIL DE LAS MUJERES: base no disponible ==')
        return
    log('\n== 9. PERFIL DE LAS 162 MUJERES ==')
    B = pd.read_csv(f, dtype=str)
    n = len(B)
    MAPAS = {
        'GRU_ED1': ('Grupo de edad', {'10': '10 a 14', '11': '15 a 19', '12': '20 a 24', '13': '25 a 29',
                                      '14': '30 a 34', '15': '35 a 39', '16': '40 a 44', '17': '45 a 49'}),
        'SEG_SOCIAL': ('Regimen de afiliacion', {'1': 'Contributivo', '2': 'Subsidiado', '5': 'No asegurada'}),
        'AREA_RES': ('Area de residencia', {'1': 'Cabecera', '2': 'Centro poblado', '3': 'Rural disperso'}),
        'IDPERTET': ('Pertenencia etnica', {'1': 'Indigena', '5': 'Afrodescendiente',
                                            '6': 'Ninguna de las anteriores'}),
        'ASIS_MED': ('Asistencia medica', {'1': 'Si', '2': 'No'}),
        'SIT_DEFUN': ('Sitio de defuncion', {'1': 'Institucion de salud', '2': 'Domicilio', '3': 'Otro'}),
    }
    filas = []
    for col, (titulo, mapa) in MAPAS.items():
        if col not in B.columns:
            continue
        v = B[col].map(mapa)
        miss = int(v.isna().sum())
        log('  %s (sin dato: %d)' % (titulo, miss))
        for k, c in v.value_counts().items():
            log('     %-28s %3d (%.1f%%)' % (k, c, c / n * 100))
            filas.append(dict(variable=titulo, categoria=k, n=int(c), pct=c / n * 100, sin_dato=miss))
    if 'subregion' in B.columns:
        log('  Subregion de residencia')
        for k, c in B.subregion.value_counts().items():
            log('     %-28s %3d (%.1f%%)' % (k, c, c / n * 100))
            filas.append(dict(variable='Subregion', categoria=k, n=int(c), pct=c / n * 100, sin_dato=0))
    pd.DataFrame(filas).to_csv(os.path.join(OUT, 'Perfil de las 162 muertes.csv'), index=False)


# ------------------------------------------------------------------ 10. vigilancia
def bloque_vigilancia(R):
    f = os.path.join(DATA, 'Notificaciones de mortalidad materna por anio.csv')
    if not os.path.exists(f):
        log('\n== 10. NOTIFICACIONES AL SISTEMA DE VIGILANCIA: base no disponible ==')
        return
    log('\n== 10. NOTIFICACIONES DE MORTALIDAD MATERNA AL SISTEMA DE VIGILANCIA ==')
    V = pd.read_csv(f)
    tot = int(V.notificadas.sum())
    log('  total notificado 2015-2024: %d' % tot)
    log('  registro civil, muertes maternas tempranas y tardias: 242')
    log('  diferencia agregada entre los dos sistemas: %+.1f%%' % ((tot / 242 - 1) * 100))
    for x in V.itertuples():
        log('     %d  notificadas %3d  |  registro civil, tempranas %3d'
            % (x.anio, int(x.notificadas), int(R[R.anio == x.anio].d_cor.iloc[0])))
    log('  las definiciones no son equivalentes: la comparacion es de magnitud agregada,')
    log('  no caso por caso, porque los registros no comparten identificador')


# ------------------------------------------------------------------ 11. duplicados
def bloque_duplicados():
    f = os.path.join(DATA, 'Base 162 muertes tempranas.csv')
    if not os.path.exists(f):
        log('  [omitido: requiere la base a nivel de registro, disponible por solicitud]'); return
    if not os.path.exists(f):
        return
    B = pd.read_csv(f, dtype=str)
    cols = [c for c in B.columns if c not in ('cod3',)]
    dup = B.duplicated(subset=cols, keep=False)
    log('\n== 11. REGISTROS IDENTICOS EN TODAS LAS VARIABLES DISPONIBLES: %d ==' % int(dup.sum()))
    if dup.sum():
        log(B[dup].to_string(index=False))
        log('  sensibilidad tratandolos como uno solo: %.1f (%.1f a %.1f) por 100 000' % ic(161, 243779))


def bloque_edad_homogeneidad(M):
    from scipy.stats import chi2
    log('\n== HOMOGENEIDAD MUNICIPAL Y ESTANDARIZACION POR EDAD (sugerencias 6.1 y 6.4) ==')
    O = M.muertes_tempranas.values.astype(float); N = M.nacidos_vivos.values.astype(float)
    theta = O.sum() / N.sum(); E = theta * N
    x2 = ((O - E) ** 2 / E).sum(); dfr = len(O) - 1
    log('  Homogeneidad Poisson chi-cuadrado = %.1f, gl %d, p = %.3f' % (x2, dfr, chi2.sf(x2, dfr)))
    fe = os.path.join(DATA, 'Estandarizacion por edad materna.csv')
    if os.path.exists(fe):
        AE = pd.read_csv(fe)
        exp = float((AE.muertes_resto / AE.nacidos_resto * AE.nacidos_cordoba).sum())
        obs = int(AE.muertes_cordoba.sum()); smr = obs / exp
        lo = chi2.ppf(.025, 2 * obs) / 2 / exp; hi = chi2.ppf(.975, 2 * (obs + 1)) / 2 / exp
        log('  Estandarizacion indirecta por edad: observadas %d, esperadas %.1f, SMR = %.2f (IC95%% %.2f a %.2f)' % (obs, exp, smr, lo, hi))


def main():
    R, C, Td, M = cargar()
    bloque_magnitud(R)
    bloque_ranking()
    bloque_tendencias(R, C, Td)
    bloque_segmentado(R)
    bloque_interacciones(R, C, Td)
    bloque_causas(R, C)
    bloque_espacial(M)
    bloque_edad_homogeneidad(M)
    bloque_ocurrencia(R)
    bloque_subregion()
    bloque_perfil()
    bloque_vigilancia(R)
    bloque_duplicados()
    with open(os.path.join(OUT, 'resultados.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(LOG))
    log('\nSalidas escritas en %s' % OUT)


if __name__ == '__main__':
    main()
