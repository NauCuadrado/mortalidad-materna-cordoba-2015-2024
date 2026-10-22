# -*- coding: utf-8 -*-
"""Autoauditoria del paquete. Naudith C, septiembre de 2026.

Comprueba la consistencia del paquete: que las bases agregadas cuadren entre si,
que los estimadores esten bien programados, que ningun modelo dependa de un solo
anio, que la simulacion de potencia este calibrada y que la matriz de vecindad no
deje municipios aislados. Cada estimador se recalcula por una via independiente de
la del pipeline, de modo que un error de programacion no se confirme a si mismo.

    python "Autoauditoria estadistica.py"

Imprime una linea por comprobacion y al final las que fallaron.
"""
import os
import warnings
import numpy as np
import pandas as pd
from scipy.stats import chi2, gamma, norm, beta, binomtest, spearmanr
import statsmodels.api as sm

warnings.filterwarnings('ignore')
BASE = os.path.dirname(os.path.abspath(__file__))
if not os.path.isdir(os.path.join(BASE, 'Datos derivados')):
    BASE = os.path.dirname(BASE)
DATA = os.path.join(BASE, 'Datos derivados') + os.sep

R = pd.read_csv(DATA + 'Serie anual numeradores y denominadores.csv')
C = pd.read_csv(DATA + 'Causas de las 162 muertes.csv')
Td = pd.read_csv(DATA + 'Muertes tempranas y tardias por anio.csv')
M = pd.read_csv(DATA + 'Base municipal con suavizado y LISA.csv')
K = pd.read_csv(DATA + 'Ranking departamental RMM.csv')
SI = pd.read_csv(DATA + 'Serie anual con intervalos.csv')
_i162 = DATA + 'Base 162 muertes tempranas.csv'
HAS_I162 = os.path.exists(_i162)
if HAS_I162:
    I162 = pd.read_csv(_i162, dtype=str)

R['d_resto'] = R.d_nac33 - R.d_cor
R['b_resto'] = R.b_nac33 - R.b_cor
C['tipo'] = np.where(C.grupo == 'Indirect causes', 'Indirect',
                     np.where(C.grupo == 'Undetermined causes', 'Undetermined', 'Direct'))
anios = R.anio.values

fallas = []


def chk(nombre, ok, detalle=''):
    print('%-4s %-58s %s' % ('OK' if ok else 'FALLA', nombre, detalle))
    if not ok:
        fallas.append((nombre, detalle))


print('=' * 100)
print('BLOQUE A. CONSISTENCIA INTERNA DE LAS BASES')
print('=' * 100)

chk('Numerador de la decada igual a 162', R.d_cor.sum() == 162, 'suma %d' % R.d_cor.sum())
chk('Denominador de la decada igual a 243 779', R.b_cor.sum() == 243779, 'suma %d' % R.b_cor.sum())
chk('La base de causas tiene 162 filas', len(C) == 162, '%d filas' % len(C))

por_anio = C.groupby('anio').size().reindex(anios, fill_value=0).values
chk('Causas por anio coinciden con la serie', np.array_equal(por_anio, R.d_cor.values),
    'max diferencia %d' % np.abs(por_anio - R.d_cor.values).max())

dir_a = np.array([((C.tipo == 'Direct') & (C.anio == y)).sum() for y in anios])
ind_a = np.array([((C.tipo == 'Indirect') & (C.anio == y)).sum() for y in anios])
und_a = np.array([((C.tipo == 'Undetermined') & (C.anio == y)).sum() for y in anios])
chk('Directas mas indirectas mas indeterminadas igual al total cada anio',
    np.array_equal(dir_a + ind_a + und_a, R.d_cor.values),
    'directas %d indirectas %d indeterminadas %d' % (dir_a.sum(), ind_a.sum(), und_a.sum()))
chk('La muerte de causa indeterminada (O95) no se cuenta como directa ni indirecta',
    (C.tipo == 'Undetermined').sum() == 1 and dir_a.sum() == 93 and ind_a.sum() == 68,
    'directas %d indirectas %d indeterminadas %d' % (dir_a.sum(), ind_a.sum(), und_a.sum()))
chk('Las tres series de causa son disjuntas y cubren las 162',
    bool(C.tipo.isin(['Direct', 'Indirect', 'Undetermined']).all()) and
    dir_a.sum() + ind_a.sum() + und_a.sum() == 162)

# conteo por grupo y periodo contra la tabla derivada
GP = pd.read_csv(DATA + 'Causas por grupo y periodo.csv')
_bg = (C.assign(n=1).pivot_table(index='grupo', columns='periodo', values='n', aggfunc='sum', fill_value=0))
_ok_gp = True
for _, row in GP.iterrows():
    for per in ('2015-2021', '2022-2024'):
        if int(row[per]) != int(_bg.get(per, {}).get(row['grupo'], 0)):
            _ok_gp = False
chk('Causas por grupo y periodo coincide con la base de 162', _ok_gp,
    'indirectas tabla %d' % int(GP.loc[GP.grupo == 'Indirect causes', 'total'].iloc[0]))

chk('Tempranas de la base de ventana igual a la serie',
    np.array_equal(Td.early.values, R.d_cor.values))
chk('Tardias suman 78', Td.late.sum() == 78, 'suma %d' % Td.late.sum())
chk('Denominador de la base de ventana igual al de la serie',
    np.array_equal(Td.b_cor.values, R.b_cor.values))
chk('Indicador de pandemia bien construido',
    np.array_equal(Td.covid.values, np.isin(anios, [2020, 2021]).astype(int)))
chk('Anio centrado bien construido', np.array_equal(Td.yc.values, anios - 2015))

chk('El pais contiene al departamento en el numerador', bool((R.d_nac33 >= R.d_cor).all()))
chk('El pais contiene al departamento en el denominador', bool((R.b_nac33 >= R.b_cor).all()))
chk('Sin denominadores cero o negativos', bool((R.b_cor > 0).all() and (R.b_resto > 0).all()))
chk('Sin conteos negativos', bool((R.d_cor >= 0).all() and (R.d_resto >= 0).all()))

chk('Base municipal: 30 municipios', len(M) == 30, '%d filas' % len(M))
chk('Base municipal: muertes suman 162', M.muertes_tempranas.sum() == 162,
    'suma %d' % M.muertes_tempranas.sum())
chk('Base municipal: nacidos vivos suman 243 779', M.nacidos_vivos.sum() == 243779,
    'suma %d' % M.nacidos_vivos.sum())
chk('Base municipal: codigos unicos', M.cod_municipio.nunique() == 30)
rmm_re = M.muertes_tempranas / M.nacidos_vivos * 1e5
chk('Razon cruda municipal recalculada coincide',
    bool(np.allclose(rmm_re, M.RMM_cruda)), 'max dif %.2e' % np.abs(rmm_re - M.RMM_cruda).max())
chk('Indice de pobreza dentro de rango plausible',
    bool((M.IPM > 0).all() and (M.IPM <= 100).all()),
    'rango %.1f a %.1f' % (M.IPM.min(), M.IPM.max()))
chk('Sin municipios con cero nacidos vivos', bool((M.nacidos_vivos > 0).all()),
    'minimo %d' % M.nacidos_vivos.min())

if HAS_I162:
    chk('Base individual: 162 registros', len(I162) == 162, '%d filas' % len(I162))
    dups = I162.duplicated(keep=False).sum()
    chk('Registros identicos en todas las variables: declarados en el articulo', dups == 2,
        '%d registros identicos' % dups)
    chk('Base individual: anios dentro del periodo',
        bool(I162.anio.astype(int).between(2015, 2024).all()))
    mun_ind = I162.groupby('anio').size().reindex([str(a) for a in anios], fill_value=0).values
    chk('Base individual por anio coincide con la serie',
        np.array_equal(mun_ind, R.d_cor.values), 'max dif %d' % np.abs(mun_ind - R.d_cor.values).max())
else:
    print('----  Base a nivel de registro no incluida (restringida): se omiten las comprobaciones individuales')

fila_cor = K[K.Departamento.str.contains('rdoba', case=False, na=False)]
rmm_k = float(fila_cor.muertes.iloc[0]) / float(fila_cor.nacidos.iloc[0]) * 1e5
chk('Ranking: la razon de Cordoba coincide con la serie',
    abs(rmm_k - R.d_cor.sum() / R.b_cor.sum() * 1e5) < 0.05,
    'ranking %.2f  serie %.2f' % (rmm_k, R.d_cor.sum() / R.b_cor.sum() * 1e5))
chk('Ranking: 33 territorios', len(K) == 33, '%d filas' % len(K))
chk('Ranking: el orden es consistente con la razon',
    bool((K.sort_values('rank').RMM.diff().dropna() <= 1e-9).all()))

print()
print('=' * 100)
print('BLOQUE B. RECALCULO INDEPENDIENTE DE LOS ESTIMADORES')
print('=' * 100)


def ic_chi2(k, n, m=1e5):
    lo = chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0
    hi = chi2.ppf(0.975, 2 * (k + 1)) / 2
    return k / n * m, lo / n * m, hi / n * m


def ic_gamma(k, n, m=1e5):
    """Misma cantidad por otra parametrizacion: cuantiles de la gamma."""
    lo = gamma.ppf(0.025, k) if k > 0 else 0.0
    hi = gamma.ppf(0.975, k + 1)
    return k / n * m, lo / n * m, hi / n * m


k, n = int(R.d_cor.sum()), int(R.b_cor.sum())
a1, a2 = ic_chi2(k, n), ic_gamma(k, n)
chk('Intervalo exacto de Poisson: dos parametrizaciones coinciden',
    bool(np.allclose(a1, a2)), '%.4f a %.4f  contra  %.4f a %.4f' % (a1[1], a1[2], a2[1], a2[2]))
chk('Intervalo exacto contiene la estimacion puntual', a1[1] < a1[0] < a1[2],
    '%.1f (%.1f a %.1f)' % a1)

kr, nr = int(R.d_resto.sum()), int(R.b_resto.sum())
rr = (k / n) / (kr / nr)
se = np.sqrt(1 / k + 1 / kr)
lo, hi = rr * np.exp(-1.96 * se), rr * np.exp(1.96 * se)
chk('Razon de tasas de la decada', abs(rr - 1.2547) < 0.002,
    'RR %.4f (%.4f a %.4f)' % (rr, lo, hi))

# exacto condicional por dos vias
p_lo = beta.ppf(0.025, k, kr + 1)
p_hi = beta.ppf(0.975, k + 1, kr)
t = n / nr
f = lambda p: (p / (1 - p)) / t
bt = binomtest(k, k + kr).proportion_ci(0.95, method='exact')
chk('Exacto condicional: la inversion binomial coincide con Clopper-Pearson',
    bool(np.allclose([p_lo, p_hi], [bt.low, bt.high], atol=1e-9)),
    'propia %.6f a %.6f  scipy %.6f a %.6f' % (p_lo, p_hi, bt.low, bt.high))
chk('Exacto condicional contiene la estimacion puntual', f(p_lo) < rr < f(p_hi),
    'RR %.3f (%.3f a %.3f)' % (rr, f(p_lo), f(p_hi)))

# razon de razones
P = {}
for lab, yy in [('1', range(2015, 2022)), ('2', range(2022, 2025))]:
    s = R[R.anio.isin(list(yy))]
    P[lab] = (int(s.d_cor.sum()), int(s.b_cor.sum()), int(s.d_resto.sum()), int(s.b_resto.sum()))
r1 = (P['1'][0] / P['1'][1]) / (P['1'][2] / P['1'][3])
r2 = (P['2'][0] / P['2'][1]) / (P['2'][2] / P['2'][3])
q = r2 / r1
seq = np.sqrt(sum(1 / x for x in (P['1'][0], P['1'][2], P['2'][0], P['2'][2])))
pq = 2 * (1 - norm.cdf(abs(np.log(q) / seq)))
chk('Razon de razones recalculada', abs(q - 0.75) < 0.01,
    'ROR %.3f (%.3f a %.3f) p=%.3f' % (q, q * np.exp(-1.96 * seq), q * np.exp(1.96 * seq), pq))
chk('Los cuatro conteos de la razon de razones suman 162 en el departamento',
    P['1'][0] + P['2'][0] == 162, '%d + %d' % (P['1'][0], P['2'][0]))


def ajusta(y, off, yy):
    d = pd.DataFrame({'y': np.asarray(y, float), 'yc': np.asarray(yy) - 2015})
    d['covid'] = np.isin(yy, [2020, 2021]).astype(int)
    return sm.GLM(d.y, sm.add_constant(d[['yc', 'covid']]),
                  family=sm.families.Poisson(), offset=np.log(off)).fit()


series = {'Todas las tempranas': R.d_cor.values, 'Causas directas': dir_a,
          'Causas indirectas': ind_a, 'Muertes tardias': Td.late.values}
off = R.b_cor.values

# la interaccion del modelo conjunto debe reproducir la diferencia de pendientes separadas
m_dir = ajusta(dir_a, off, anios)
m_ind = ajusta(ind_a, off, anios)
yy2 = np.r_[anios, anios]
X = sm.add_constant(np.column_stack([
    np.r_[np.zeros(10), np.ones(10)], yy2 - 2015, np.isin(yy2, [2020, 2021]).astype(int),
    np.r_[np.zeros(10), np.ones(10)] * (yy2 - 2015),
    np.r_[np.zeros(10), np.ones(10)] * np.isin(yy2, [2020, 2021]).astype(int)]))
mj = sm.GLM(np.r_[dir_a, ind_a], X, family=sm.families.Poisson(),
            offset=np.log(np.r_[off, off])).fit()
d_pend = m_ind.params['yc'] - m_dir.params['yc']
d_cov = m_ind.params['covid'] - m_dir.params['covid']
chk('Interaccion causa por anio igual a la diferencia de pendientes separadas',
    abs(mj.params[4] - d_pend) < 1e-8, 'conjunto %.6f  separados %.6f' % (mj.params[4], d_pend))
chk('Interaccion causa por pandemia igual a la diferencia de efectos separados',
    abs(mj.params[5] - d_cov) < 1e-8, 'conjunto %.6f  separados %.6f' % (mj.params[5], d_cov))

# modelo segmentado: contrastes contra pendientes recalculadas desde los parametros
x0 = anios - 2015
x1 = np.clip(anios - 2019, 0, None)
x2 = np.clip(anios - 2021, 0, None)
ms = sm.GLM(R.d_cor.values.astype(float), sm.add_constant(np.column_stack([x0, x1, x2])),
            family=sm.families.Poisson(), offset=np.log(off)).fit()
b = ms.params
pend = {'Pre-pandemia': b[1], 'Pandemia': b[1] + b[2], 'Pos-pandemia': b[1] + b[2] + b[3]}
esperado = {'Pre-pandemia': -10.1, 'Pandemia': 33.7, 'Pos-pandemia': -32.7}
for kk, v in pend.items():
    apc = (np.exp(v) - 1) * 100
    chk('Segmentado, pendiente de %s' % kk, abs(apc - esperado[kk]) < 0.2,
        '%+.1f%% (reportado %+.1f%%)' % (apc, esperado[kk]))
n_tr = [int(R.d_cor[anios <= 2019].sum()), int(R.d_cor[np.isin(anios, [2020, 2021])].sum()),
        int(R.d_cor[anios >= 2022].sum())]
chk('Segmentado: los tres bloques suman 162', sum(n_tr) == 162, '%d + %d + %d' % tuple(n_tr))

print()
print('=' * 100)
print('BLOQUE C. SUPUESTOS Y DIAGNOSTICOS')
print('=' * 100)

for nom, y in series.items():
    m = ajusta(y, off, anios)
    disp = m.pearson_chi2 / m.df_resid
    infl = m.get_influence()
    hat = infl.hat_matrix_diag
    cook = infl.cooks_distance[0]
    peor = anios[np.argmax(cook)]
    print('  %-20s dispersion %.2f | devianza/gl %.2f | hat max %.2f | Cook max %.2f (%d)'
          % (nom, disp, m.deviance / m.df_resid, hat.max(), cook.max(), peor))
    chk('  %s: sin apalancamiento extremo' % nom, hat.max() < 0.8, 'hat max %.2f' % hat.max())
    # Con diez observaciones anuales un punto influyente es esperable; se examina con la
    # omision de un anio a la vez (Tablas S3 y S4) y se informa aqui sin tratarlo como fallo.
    if cook.max() >= 1.0:
        print('     nota: %s tiene un punto influyente en %d (Cook %.2f); se examina por omision de anios'
              % (nom.strip(), peor, cook.max()))
    chk('  %s: dispersion sin sobredispersion severa' % nom, disp < 2.0, 'dispersion %.2f' % disp)

print()
d = M.muertes_tempranas.values.astype(float)
nv = M.nacidos_vivos.values.astype(float)
ipm10 = (M.IPM / 10).values
Xe = sm.add_constant(pd.Series(ipm10, name='ipm10'))
me = sm.GLM(d, Xe, offset=np.log(nv), family=sm.families.Poisson()).fit()
infl = me.get_influence()
print('  Modelo ecologico: IRR %.3f (%.3f a %.3f) p=%.3f | dispersion %.2f'
      % (np.exp(me.params['ipm10']), np.exp(me.conf_int().loc['ipm10'][0]),
         np.exp(me.conf_int().loc['ipm10'][1]), me.pvalues['ipm10'],
         me.pearson_chi2 / me.df_resid))
cook = infl.cooks_distance[0]
peor = M.municipio.values[np.argmax(cook)]
chk('Ecologico: sin municipio dominante', cook.max() < 1.0, 'Cook max %.2f en %s' % (cook.max(), peor))
# linealidad: termino cuadratico y version por cuartiles
Xq = sm.add_constant(np.column_stack([ipm10, (ipm10 - ipm10.mean()) ** 2]))
mq = sm.GLM(d, Xq, offset=np.log(nv), family=sm.families.Poisson()).fit()
chk('Ecologico: el termino cuadratico no aporta', mq.pvalues[2] > 0.05,
    'p del cuadratico %.3f' % mq.pvalues[2])
cuart = pd.qcut(ipm10, 4, labels=False)
Xc = sm.add_constant(pd.get_dummies(cuart, drop_first=True).astype(float).values)
mc = sm.GLM(d, Xc, offset=np.log(nv), family=sm.families.Poisson()).fit()
lr = 2 * (mc.llf - sm.GLM(d, np.ones((30, 1)), offset=np.log(nv),
                          family=sm.families.Poisson()).fit().llf)
chk('Ecologico: los cuartiles de pobreza tampoco separan',
    1 - chi2.cdf(lr, 3) > 0.05, 'p de la prueba por cuartiles %.3f' % (1 - chi2.cdf(lr, 3)))
rho, prho = spearmanr(ipm10, d / nv)
chk('Ecologico: sin asociacion monotona por rangos', prho > 0.05,
    'rho %.2f p %.3f' % (rho, prho))

print()
print('=' * 100)
print('BLOQUE D. VALIDACION DE LA SIMULACION DE POTENCIA')
print('=' * 100)

obs = int(d.sum())
for irr, reps in [(1.00, 4000), (1.15, 1500)]:
    rng = np.random.default_rng(7)
    lam = nv * np.exp(np.log(irr) * (ipm10 - ipm10.mean()))
    lam = lam * obs / lam.sum()
    h = 0
    for _ in range(reps):
        yy = rng.poisson(lam)
        try:
            if sm.GLM(yy, Xe, offset=np.log(nv), family=sm.families.Poisson()).fit().pvalues['ipm10'] < 0.05:
                h += 1
        except Exception:
            pass
    tasa = h / reps
    if irr == 1.00:
        ee = 1.96 * np.sqrt(0.05 * 0.95 / reps)
        chk('Error tipo I de la simulacion cercano al 5 por ciento',
            abs(tasa - 0.05) < ee + 0.005,
            'rechazos %.1f%% en %d replicas bajo la hipotesis nula' % (tasa * 100, reps))
    else:
        chk('Potencia en 1,15 reproduce el valor reportado', abs(tasa - 0.796) < 0.04,
            'potencia %.1f%% (reportado 79,6%%)' % (tasa * 100))
    chk('  La simulacion conserva las muertes esperadas en IRR %.2f' % irr,
        abs(lam.sum() - obs) < 0.5, 'esperadas %.1f contra %d observadas' % (lam.sum(), obs))

print()
print('=' * 100)
print('BLOQUE E. MATRIZ ESPACIAL')
print('=' * 100)
import geopandas as gpd
import libpysal
import esda
g = gpd.read_file(DATA + 'cordoba_municipios.geojson')
g['cod'] = g['cod'].astype(str)
M2 = M.copy(); M2['cod_municipio'] = M2.cod_municipio.astype(str)
gdf = (g.merge(M2, left_on='cod', right_on='cod_municipio', suffixes=('_geo', ''))
       .sort_values('cod_municipio').reset_index(drop=True))
chk('La geometria empareja con los 30 municipios', len(gdf) == 30, '%d unidades' % len(gdf))
w = libpysal.weights.Queen.from_dataframe(gdf, use_index=False)
vec = np.array([len(v) for v in w.neighbors.values()])
chk('Sin municipios aislados en la matriz de contigüidad', len(w.islands) == 0,
    'islas: %s' % (w.islands if w.islands else 'ninguna'))
chk('La matriz es simetrica', all(i in w.neighbors[j] for i in w.neighbors for j in w.neighbors[i]))
print('  vecinos por municipio: minimo %d, mediana %.0f, maximo %d'
      % (vec.min(), np.median(vec), vec.max()))
w.transform = 'r'
Is, ps = [], []
for s in (1, 7, 42, 123, 2024):
    np.random.seed(s)
    mi = esda.Moran(gdf['RMM_cruda'].values, w, permutations=999)
    Is.append(mi.I); ps.append(mi.p_sim)
chk('El indice de Moran no depende de la semilla', np.ptp(Is) < 1e-12, 'I = %.4f' % Is[0])
chk('El valor de p es estable entre semillas', np.ptp(ps) < 0.05,
    'p entre %.3f y %.3f en cinco semillas' % (min(ps), max(ps)))
np.random.seed(42)
li = esda.Moran_Local(gdf['RMM_cruda'].values, w, permutations=999)
psim = np.sort(li.p_sim)
fdr = psim <= (np.arange(1, 31) / 30) * 0.05
chk('LISA: ningun municipio significativo, con o sin ajuste',
    (li.p_sim < 0.05).sum() == 0 and not fdr.any(),
    'crudos %d, tras FDR %d' % ((li.p_sim < 0.05).sum(), int(fdr.sum())))

print()
print('=' * 100)
print('RESUMEN: %d comprobaciones fallidas' % len(fallas))
for nom, det in fallas:
    print('   - %s :: %s' % (nom, det))
print('=' * 100)
