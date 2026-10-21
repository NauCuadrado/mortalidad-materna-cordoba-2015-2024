# -*- coding: utf-8 -*-
"""
Código reproducible del artículo de mortalidad materna en Córdoba, 2015 a 2024.
Autor: Naudith C. Fecha: 31 de julio de 2026.

Reconstruye las bases de "Datos derivados" a partir de los microdatos públicos del DANE
(estadísticas vitales) y del INS (SIVIGILA, evento de mortalidad materna). No requiere
datos restringidos, pero sí que esos archivos originales estén descargados: el paquete
no los redistribuye.

Uso:
    python "Codigo reproducible mortalidad.py" RUTA_DE_LOS_MICRODATOS

    Sin argumento toma la variable de entorno BASE_MORTALIDAD. Si la carpeta no existe o
    no contiene los archivos anuales, el script lo dice y termina sin escribir nada.

    Las bases se escriben siempre en la carpeta "Datos derivados" del paquete, calculada a
    partir de la ubicacion de este archivo, sea cual sea el directorio desde el que se ejecute.
    Para escribir en otro sitio, defina BASE_MORTALIDAD apuntando a la carpeta que debe
    contener "Datos derivados".

Una vez reconstruidas las bases derivadas, las cifras del manuscrito se obtienen con
"pipeline.py", que es el script que produce Salidas/resultados.txt.

Entradas esperadas (rutas relativas a BASE):
  01_Datos_fuente/DANE/          archivos anuales de defunciones no fetales y de nacimientos
  01_Datos_fuente/SIVIGILA/      Datos_<año>_550.xls o .xlsx
  01_Datos_fuente/Covariables_geometria/cordoba_municipios.geojson
  02_Bases_analiticas/Ecologico_RMM_IPM_Cordoba_Naudith C.csv   (índice de pobreza municipal 2018)

Salidas: los CSV de la carpeta "Datos derivados" del paquete de trazabilidad.
"""
import os, glob, json, sys
import numpy as np
import pandas as pd
from scipy.stats import chi2
import statsmodels.api as sm

# La salida se define respecto de la carpeta del paquete, no del directorio desde el que se
# ejecute el programa. La variable de entorno BASE_MORTALIDAD permite redirigirla.
_AQUI = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get('BASE_MORTALIDAD') or (
    _AQUI if os.path.isdir(os.path.join(_AQUI, 'Datos derivados')) else os.path.dirname(_AQUI))
OUT = os.path.join(BASE, 'Datos derivados')
os.makedirs(OUT, exist_ok=True)
SEED = 42

DEPN = [5, 8, 11, 13, 15, 17, 18, 19, 20, 23, 25, 27, 41, 44, 47, 50, 52, 54, 63, 66, 68, 70, 73, 76,
        81, 85, 86, 88, 91, 94, 95, 97, 99]


NOMBRE_DPTO = {5: 'Antioquia', 8: 'Atlántico', 11: 'Bogotá', 13: 'Bolívar', 15: 'Boyacá', 17: 'Caldas',
               18: 'Caquetá', 19: 'Cauca', 20: 'Cesar', 23: 'Córdoba', 25: 'Cundinamarca', 27: 'Chocó',
               41: 'Huila', 44: 'La Guajira', 47: 'Magdalena', 50: 'Meta', 52: 'Nariño',
               54: 'N. Santander', 63: 'Quindío', 66: 'Risaralda', 68: 'Santander', 70: 'Sucre',
               73: 'Tolima', 76: 'Valle', 81: 'Arauca', 85: 'Casanare', 86: 'Putumayo',
               88: 'San Andrés', 91: 'Amazonas', 94: 'Guainía', 95: 'Guaviare', 97: 'Vaupés',
               99: 'Vichada'}


# subregion de cada municipio del departamento, segun la division subregional oficial
SUBREGION = {'001': 'Medio Sinu', '068': 'San Jorge', '079': 'San Jorge', '090': 'Costanera', 
             '162': 'Medio Sinu', '168': 'Bajo Sinu', '182': 'Sabanas', '189': 'Medio Sinu', 
             '300': 'Bajo Sinu', '350': 'San Jorge', '417': 'Bajo Sinu', '419': 'Costanera', 
             '464': 'Bajo Sinu', '466': 'San Jorge', '500': 'Costanera', '555': 'San Jorge', 
             '570': 'San Jorge', '574': 'Costanera', '580': 'San Jorge', '586': 'Bajo Sinu', 
             '660': 'Sabanas', '670': 'Sabanas', '672': 'Costanera', '675': 'Costanera', 
             '678': 'Medio Sinu', '682': 'San Jorge', '686': 'Medio Sinu', '807': 'Alto Sinu', 
             '815': 'Sabanas', '855': 'Alto Sinu'}


# etiqueta en ingles de cada grupo de causa, que es la que consumen el pipeline y los generadores
GRUPO_EN = {'Trastornos hipertensivos': 'Hypertensive disorders',
            'Hemorragia obstetrica': 'Obstetric hemorrhage',
            'Aborto': 'Abortive outcomes',
            'Infeccion relacionada con el embarazo': 'Pregnancy-related infection',
            'Embolia obstetrica': 'Obstetric embolism',
            'Otras causas directas': 'Other direct obstetric causes',
            'Causas indirectas': 'Indirect causes',
            'Muerte materna no especificada': 'Unspecified maternal death'}


def _sin_tildes(t):
    import unicodedata
    return ''.join(c for c in unicodedata.normalize('NFD', t) if unicodedata.category(c) != 'Mn')

# ---------------------------------------------------------------- 1. definiciones
def clase_muerte(cb):
    """Devuelve 'temprana', 'tardia' o None segun la causa basica CIE-10."""
    cb = str(cb).strip().upper().replace('"', '')
    if cb == 'A34':
        return 'temprana'
    if cb[:1] != 'O':
        return None
    if cb[:3] == 'O96':
        return 'tardia'
    if cb[:3] == 'O97':
        return None  # secuelas a un anio o mas: no es muerte materna tardia (CIE-10), se excluye
    return 'temprana'

def grupo_causa(cb):
    cb = str(cb).strip().upper().replace('"', '')
    if cb == 'A34':
        return 'Infeccion relacionada con el embarazo'
    n = cb[1:3]
    try:
        n = int(n)
    except ValueError:
        return 'Otras causas directas'
    if 0 <= n <= 8:
        return 'Aborto'
    if 10 <= n <= 16:
        return 'Trastornos hipertensivos'
    if n in (20, 44, 45, 46, 67, 72) or cb[:4] in ('O710', 'O711', 'O712', 'O713'):
        return 'Hemorragia obstetrica'
    if n in (85, 86, 91):
        return 'Infeccion relacionada con el embarazo'
    if n == 88:
        return 'Embolia obstetrica'
    if cb[:4] in ('O240', 'O241', 'O242', 'O243'):
        return 'Causas indirectas'  # diabetes preexistente: indirecta segun CIE-MM (grupo 7)
    if n in (98, 99):
        return 'Causas indirectas'
    if n == 95:
        return 'Muerte materna no especificada'
    return 'Otras causas directas'

def ic_poisson(k, n, mult=1e5):
    lo = chi2.ppf(0.025, 2 * k) / 2 if k > 0 else 0.0
    hi = chi2.ppf(0.975, 2 * (k + 1)) / 2
    return k / n * mult, lo / n * mult, hi / n * mult

def razon_tasas(k1, n1, k2, n2):
    r = (k2 / n2) / (k1 / n1)
    se = np.sqrt(1 / k1 + 1 / k2)
    return r, r * np.exp(-1.96 * se), r * np.exp(1.96 * se)

# ---------------------------------------------------------------- 2. lectura de microdatos
def leer_defunciones(anio, ruta):
    cols = ('C_BAS1', 'SEXO', 'CODPTORE', 'CODMUNRE', 'COD_DPTO', 'COD_MUNIC', 'GRU_ED1', 'EST_CIVIL',
            'NIVEL_EDU', 'IDPERTET', 'AREA_RES', 'SEG_SOCIAL', 'SIT_DEFUN', 'ASIS_MED', 'EMB_FAL')
    ruta = _csv_de(ruta)
    df = pd.read_csv(ruta, usecols=lambda c: str(c).strip().replace('"', '') in cols,
                     dtype=str, encoding='latin-1', low_memory=False)
    df.columns = [c.strip().replace('"', '') for c in df.columns]
    for c in df.columns:
        df[c] = df[c].astype(str).str.replace('"', '').str.strip()
    df = df[df['SEXO'] == '2']
    df['anio'] = anio
    df['clase'] = df['C_BAS1'].map(clase_muerte)
    return df

def leer_nacimientos(anio, ruta):
    ruta = _csv_de(ruta)
    df = pd.read_csv(ruta, usecols=['CODPTORE', 'CODMUNRE'], dtype=str,
                     encoding='latin-1', low_memory=False)
    for c in df.columns:
        df[c] = df[c].astype(str).str.replace('"', '').str.strip()
    df['anio'] = anio
    return df


# ---------------------------------------------------------------- localizacion de archivos
def _csv_de(ruta):
    """Devuelve una ruta de CSV legible. Si el archivo es un ZIP, extrae el CSV que contiene."""
    if ruta.lower().endswith('.csv'):
        return ruta
    if ruta.lower().endswith('.zip'):
        import zipfile, tempfile
        with zipfile.ZipFile(ruta) as z:
            miembros = [m for m in z.namelist() if m.lower().endswith('.csv')]
            if not miembros:
                raise SystemExit('El archivo %s no contiene ningun CSV. Descomprimalo un nivel.' % ruta)
            destino = os.path.join(tempfile.mkdtemp(prefix='dane_'), os.path.basename(miembros[0]))
            with z.open(miembros[0]) as src, open(destino, 'wb') as dst:
                dst.write(src.read())
            return destino
    raise SystemExit('Formato no reconocido: %s' % ruta)


def _buscar(carpeta, patrones):
    """Busca de forma recursiva, sin distinguir mayusculas, el primer archivo que cumpla
       alguno de los patrones. Los nombres del portal del DANE varian entre anios."""
    import fnmatch
    todos = []
    for raiz, _, archivos in os.walk(carpeta):
        todos += [os.path.join(raiz, a) for a in archivos]
    todos.sort()
    for pat in patrones:
        hit = [f for f in todos if fnmatch.fnmatch(os.path.basename(f).lower(), pat.lower())]
        if hit:
            return hit[0]
    return None

# ---------------------------------------------------------------- 3. serie anual
def serie_anual(DEF, NAC):
    filas = []
    for y in range(2015, 2025):
        d = DEF[DEF.anio == y]; b = NAC[NAC.anio == y]
        dp = pd.to_numeric(d.CODPTORE, errors='coerce'); bp = pd.to_numeric(b.CODPTORE, errors='coerce')
        filas.append(dict(
            anio=y,
            d_cor=int(((dp == 23) & (d.clase == 'temprana')).sum()),
            tardias_cor=int(((dp == 23) & (d.clase == 'tardia')).sum()),
            d_nac33=int((dp.isin(DEPN) & (d.clase == 'temprana')).sum()),
            b_cor=int((bp == 23).sum()),
            b_nac33=int(bp.isin(DEPN).sum()),
            b_nacall=int(len(b))))
    S = pd.DataFrame(filas)
    for lado, k, n, m in [('cor', 'd_cor', 'b_cor', 1e5), ('nac', 'd_nac33', 'b_nac33', 1e5)]:
        r = [ic_poisson(a, c, m) for a, c in zip(S[k], S[n])]
        S['rmm_' + lado] = [x[0] for x in r]
        S['rmm_%s_lo' % lado] = [x[1] for x in r]
        S['rmm_%s_hi' % lado] = [x[2] for x in r]
    S.to_csv(os.path.join(OUT, 'Serie anual mortalidad.csv'), index=False)
    return S

# ---------------------------------------------------------------- 4. modelos de tendencia
def tendencia(y, offset, anios, etiqueta):
    df = pd.DataFrame({'y': y, 'anio': anios})
    df['yc'] = df.anio - df.anio.min()
    df['covid'] = df.anio.isin([2020, 2021]).astype(int)
    X = sm.add_constant(df[['yc', 'covid']])
    m = sm.GLM(df.y, X, family=sm.families.Poisson(), offset=np.log(offset)).fit()
    apc = (np.exp(m.params['yc']) - 1) * 100
    ci = (np.exp(m.conf_int().loc['yc']) - 1) * 100
    rr = np.exp(m.params['covid']); rci = np.exp(m.conf_int().loc['covid'])
    return dict(modelo=etiqueta, apc=apc, apc_lo=ci[0], apc_hi=ci[1], p_apc=m.pvalues['yc'],
                covid_rr=rr, covid_lo=rci[0], covid_hi=rci[1], p_covid=m.pvalues['covid'],
                dispersion=m.pearson_chi2 / m.df_resid)

def tendencia_segmentada(y, offset, anios):
    df = pd.DataFrame({'y': y, 'anio': anios}); df['yc'] = df.anio - 2015
    df['s1'] = np.minimum(df.yc, 4)
    df['s2'] = np.clip(df.yc - 4, 0, 2)
    df['s3'] = np.clip(df.yc - 6, 0, None)
    X = sm.add_constant(df[['s1', 's2', 's3']])
    m = sm.GLM(df.y, X, family=sm.families.Poisson(), offset=np.log(offset)).fit()
    out = []
    for seg, lab in [('s1', 'Prepandemia 2015-2019'), ('s2', 'Pandemia 2020-2021'), ('s3', 'Pospandemia 2022-2024')]:
        a = (np.exp(m.params[seg]) - 1) * 100; ci = (np.exp(m.conf_int().loc[seg]) - 1) * 100
        out.append(dict(segmento=lab, apc=a, lo=ci[0], hi=ci[1], p=m.pvalues[seg]))
    return pd.DataFrame(out)

# ---------------------------------------------------------------- 5. espacial y ecologico
def espacial_y_ecologico(MUN, geojson, permutaciones=999):
    import geopandas as gpd, libpysal, esda
    g = gpd.read_file(geojson); g['cod'] = g['cod'].astype(str)
    gdf = g.merge(MUN, left_on='cod', right_on='cod_municipio', how='left').sort_values('cod').reset_index(drop=True)
    dd = gdf.muertes.values.astype(float); nn = gdf.nacidos_vivos.values.astype(float)
    m = dd.sum() / nn.sum()
    s2 = max((nn * ((dd / nn) - m) ** 2).sum() / nn.sum() - m / nn.mean(), 0)
    C = s2 / (s2 + m / nn)
    gdf['rmm_cruda'] = dd / nn * 1e5
    gdf['rmm_eb_global'] = (m + C * ((dd / nn) - m)) * 1e5
    w = libpysal.weights.Queen.from_dataframe(gdf, use_index=False); w.transform = 'r'
    res = {}
    for col in ['rmm_cruda', 'rmm_eb_global']:
        np.random.seed(SEED)
        mi = esda.Moran(gdf[col].values, w, permutations=permutaciones)
        res[col] = dict(I=mi.I, EI=mi.EI, p=mi.p_sim, p95=float(np.percentile(mi.sim, 95)))
    np.random.seed(SEED)
    lisa = esda.Moran_Local(gdf['rmm_cruda'].values, w, permutations=permutaciones)
    res['lisa_significativos'] = int((lisa.p_sim < 0.05).sum())
    gdf['ipm10'] = gdf.IPM / 10
    X = sm.add_constant(gdf['ipm10']); off = np.log(gdf.nacidos_vivos)
    mod = sm.GLM(gdf.muertes, X, offset=off, family=sm.families.Poisson()).fit()
    res['ecologico'] = dict(irr=np.exp(mod.params['ipm10']),
                            lo=np.exp(mod.conf_int().loc['ipm10'])[0],
                            hi=np.exp(mod.conf_int().loc['ipm10'])[1],
                            p=mod.pvalues['ipm10'],
                            dispersion=mod.pearson_chi2 / mod.df_resid)
    return gdf, res

def potencia_ecologica(MUN, irr, replicas=3000, semilla=7):
    rng = np.random.default_rng(semilla)
    base = MUN.muertes.sum() / MUN.nacidos_vivos.sum()
    ipm10 = (MUN.IPM / 10).values; nv = MUN.nacidos_vivos.values
    X = sm.add_constant(pd.Series(ipm10, name='ipm10')); off = np.log(nv)
    lam = nv * base * np.exp(np.log(irr) * (ipm10 - ipm10.mean()))
    hits = 0
    for _ in range(replicas):
        y = rng.poisson(lam)
        try:
            if sm.GLM(y, X, offset=off, family=sm.families.Poisson()).fit().pvalues['ipm10'] < 0.05:
                hits += 1
        except Exception:
            pass
    return hits / replicas


# ---------------------------------------------------------------- 7. ranking territorial
def ranking_departamental(DEF, NAC):
    """RMM acumulada de la decada por territorio, con su posicion (Figura S1)."""
    dp = pd.to_numeric(DEF.CODPTORE, errors='coerce')
    bp = pd.to_numeric(NAC.CODPTORE, errors='coerce')
    filas = []
    for cod in DEPN:
        k = int(((dp == cod) & (DEF.clase == 'temprana')).sum())
        n = int((bp == cod).sum())
        if n:
            filas.append(dict(dpto=cod, Departamento=NOMBRE_DPTO.get(cod, str(cod)),
                              muertes=k, nacidos=n, RMM=round(k / n * 1e5, 1)))
    K = pd.DataFrame(filas).sort_values('RMM', ascending=False).reset_index(drop=True)
    K['rank'] = np.arange(1, len(K) + 1)
    return K


# ---------------------------------------------------------------- 8. base municipal
def base_municipal(DEF, NAC, ipm_csv, geojson=None):
    """Muertes y nacimientos por municipio de residencia del departamento, con el indice de
       pobreza, el suavizado empirico de Bayes global y, si hay geometria, el LISA."""
    dp = pd.to_numeric(DEF.CODPTORE, errors='coerce')
    bp = pd.to_numeric(NAC.CODPTORE, errors='coerce')
    E = DEF[(dp == 23) & (DEF.clase == 'temprana')]
    N = NAC[bp == 23]
    mu = E.CODMUNRE.astype(str).str.zfill(3).value_counts()
    nv = N.CODMUNRE.astype(str).str.zfill(3).value_counts()
    M = pd.DataFrame({'cod3': nv.index, 'nacidos_vivos': nv.values})
    M['muertes_tempranas'] = M.cod3.map(mu).fillna(0).astype(int)
    M['cod_municipio'] = '23' + M.cod3
    IPM = pd.read_csv(ipm_csv)
    IPM['cod_municipio'] = IPM.cod_municipio.astype(str)
    M = M.merge(IPM[['cod_municipio', 'municipio', 'subregion', 'IPM']], on='cod_municipio', how='left')
    d = M.muertes_tempranas.values.astype(float); n = M.nacidos_vivos.values.astype(float)
    m = d.sum() / n.sum()
    s2 = max((n * ((d / n) - m) ** 2).sum() / n.sum() - m / n.mean(), 0)   # Marshall global
    C = s2 / (s2 + m / n)
    M['RMM_cruda'] = d / n * 1e5
    M['RMM_EBglobal'] = (m + C * ((d / n) - m)) * 1e5
    M['lisa_p'] = np.nan
    M['lisa_q'] = np.nan
    if geojson and os.path.exists(geojson):
        try:
            import geopandas as gpd, libpysal, esda
            g = gpd.read_file(geojson); g['cod'] = g['cod'].astype(str)
            gdf = g.merge(M, left_on='cod', right_on='cod_municipio').sort_values('cod_municipio')
            gdf = gdf.reset_index(drop=True)
            w = libpysal.weights.Queen.from_dataframe(gdf, use_index=False); w.transform = 'r'
            np.random.seed(SEED)
            li = esda.Moran_Local(gdf.RMM_cruda.values, w, permutations=999)
            lp = pd.DataFrame({'cod_municipio': gdf.cod_municipio.values,
                               'p': li.p_sim, 'q': li.q})
            M = M.drop(columns=['lisa_p', 'lisa_q']).merge(lp, on='cod_municipio', how='left')
            M = M.rename(columns={'p': 'lisa_p', 'q': 'lisa_q'})
        except ImportError:
            print('   geopandas, libpysal o esda no instalados; la base municipal sale sin LISA')
    return M.sort_values('cod_municipio').reset_index(drop=True)


# ---------------------------------------------------------------- 6. flujo de referencia
def flujo_referencia(DEF):
    dp = pd.to_numeric(DEF.CODPTORE, errors='coerce'); oc = pd.to_numeric(DEF.COD_DPTO, errors='coerce')
    prop = DEF[(dp == 23) & (DEF.clase == 'temprana')]
    salen = prop[pd.to_numeric(prop.COD_DPTO, errors='coerce') != 23]
    entran = DEF[(oc == 23) & (dp != 23) & (DEF.clase == 'temprana')]
    return dict(propias=len(prop), dentro=len(prop) - len(salen), salen=len(salen), entran=len(entran),
                destino_salen=salen.COD_DPTO.value_counts().to_dict(),
                origen_entran=entran.CODPTORE.value_counts().to_dict(),
                numerador_por_ocurrencia=len(prop) - len(salen) + len(entran))


def flujo_referencia_tabla(DEF):
    """El mismo flujo en la forma larga de "Flujo de referencia.csv"."""
    f = flujo_referencia(DEF)
    filas = [dict(flujo='Ocurrieron en Cordoba', departamento='Cordoba', muertes=f['dentro'])]
    for etiqueta, conteos in [('Salieron de Cordoba', f['destino_salen']),
                              ('Entraron a Cordoba', f['origen_entran'])]:
        for cod, k in sorted(conteos.items(), key=lambda x: -x[1]):
            c = int(pd.to_numeric(cod, errors='coerce'))
            filas.append(dict(flujo=etiqueta,
                              departamento=_sin_tildes(NOMBRE_DPTO.get(c, str(cod))),
                              muertes=int(k)))
    return pd.DataFrame(filas)

# ---------------------------------------------------------------- 7. verificacion con SIVIGILA
def conteo_sivigila_550(carpeta):
    out = {}
    for y in range(2015, 2025):
        f = glob.glob(os.path.join(carpeta, 'Datos_%d_550.*' % y))[0]
        df = pd.read_excel(f, engine='calamine', dtype=str)
        df.columns = [c.strip() for c in df.columns]
        dp = pd.to_numeric(df['COD_DPTO_R'], errors='coerce')
        out[y] = dict(nacional=int(len(df)), cordoba=int((dp == 23).sum()))
    return out

# ---------------------------------------------------------------- ejecucion
def reconstruir(base):
    """Ejecuta el flujo completo desde los microdatos y escribe las bases derivadas."""
    dane = os.path.join(base, '01_Datos_fuente', 'DANE')
    sivi = os.path.join(base, '01_Datos_fuente', 'SIVIGILA')
    geo = os.path.join(base, '01_Datos_fuente', 'Covariables_geometria', 'cordoba_municipios.geojson')
    if not os.path.isdir(dane):
        print('No se encuentra "%s".' % dane)
        print('Descargue los microdatos anuales del DANE y del SIVIGILA y vuelva a ejecutar:')
        print('    python "Codigo reproducible mortalidad.py" RUTA_DE_LOS_MICRODATOS')
        return 1

    print('1. Leyendo defunciones no fetales y nacimientos, 2015 a 2024')
    DEF, NAC = [], []
    for y in range(2015, 2025):
        fd = _buscar(dane, ['nofetal%d.csv' % y, 'nofetal%d.zip' % y,
                            '*nofetal*%d*.csv' % y, '*nofetal*%d*.zip' % y])
        fn = _buscar(dane, ['nac%d.csv' % y, 'nac%d.zip' % y, 'nacimientos%d.csv' % y,
                            'nacimientos%d.zip' % y, '*nacimiento*%d*.csv' % y,
                            '*nacimiento*%d*.zip' % y])
        if not fd or not fn:
            print('   falta el archivo anual de %d. Los archivos deben llamarse nofetal<anio> y '
                  'nac<anio>, en .csv o en .zip, dentro de la carpeta DANE (los ZIP anidados del '
                  'portal deben descomprimirse un nivel).' % y)
            return 1
        print('   %d: %s | %s' % (y, os.path.basename(fd), os.path.basename(fn)))
        DEF.append(leer_defunciones(y, fd))
        NAC.append(leer_nacimientos(y, fn))
    DEF = pd.concat(DEF, ignore_index=True)
    NAC = pd.concat(NAC, ignore_index=True)

    print('2. Serie anual de numeradores y denominadores')
    S = serie_anual(DEF, NAC)
    S[['anio', 'd_nac33', 'd_cor', 'b_nac33', 'b_cor', 'b_nacall']].to_csv(
        os.path.join(OUT, 'Serie anual numeradores y denominadores.csv'), index=False)
    S.to_csv(os.path.join(OUT, 'Serie anual con intervalos.csv'), index=False)

    print('3. Base de las muertes tempranas y agrupacion por causa')
    dp = pd.to_numeric(DEF.CODPTORE, errors='coerce')
    E = DEF[(dp == 23) & (DEF.clase == 'temprana')].copy()
    E['grupo'] = E.C_BAS1.map(grupo_causa).map(lambda g: GRUPO_EN.get(g, g))
    E['cod3'] = E.CODMUNRE.astype(str).str.zfill(3)
    E['subregion'] = E.cod3.map(SUBREGION)
    E['periodo'] = np.where(E.anio <= 2021, '2015-2021', '2022-2024')
    COLS = ['anio', 'CODMUNRE', 'cod3', 'subregion', 'GRU_ED1', 'EST_CIVIL', 'NIVEL_EDU', 'IDPERTET',
            'AREA_RES', 'SEG_SOCIAL', 'SIT_DEFUN', 'ASIS_MED', 'EMB_FAL', 'COD_DPTO', 'COD_MUNIC', 'C_BAS1']
    E[COLS].to_csv(os.path.join(OUT, 'Base 162 muertes tempranas.csv'), index=False)
    E.assign(cod=E.C_BAS1)[['anio', 'C_BAS1', 'cod', 'grupo', 'periodo']].to_csv(
        os.path.join(OUT, 'Causas de las 162 muertes.csv'), index=False)
    G = (E.pivot_table(index='grupo', columns='periodo', aggfunc='size', fill_value=0)
           .reindex(columns=['2015-2021', '2022-2024'], fill_value=0))
    G['total'] = G.sum(axis=1)
    G.to_csv(os.path.join(OUT, 'Causas por grupo y periodo.csv'))

    print('4. Muertes tempranas y tardias por anio')
    T = S[['anio']].copy()
    T['early'] = S.d_cor.values
    T['late'] = S.tardias_cor.values
    T['b_cor'] = S.b_cor.values
    T['total'] = T.early + T.late
    T['yc'] = T.anio - 2015
    T['covid'] = T.anio.isin([2020, 2021]).astype(int)
    T.to_csv(os.path.join(OUT, 'Muertes tempranas y tardias por anio.csv'), index=False)

    print('5. Flujo de ocurrencia')
    flujo_referencia_tabla(DEF).to_csv(os.path.join(OUT, 'Flujo de referencia.csv'), index=False)

    print('6. Orden de los 33 territorios')
    ranking_departamental(DEF, NAC).to_csv(
        os.path.join(OUT, 'Ranking departamental RMM.csv'), index=False)

    print('7. Base municipal con suavizado y LISA')
    ipm = os.path.join(base, '02_Bases_analiticas', 'Ecologico_RMM_IPM_Cordoba_Naudith C.csv')
    if os.path.exists(ipm):
        base_municipal(DEF, NAC, ipm, geo).to_csv(
            os.path.join(OUT, 'Base municipal con suavizado y LISA.csv'), index=False)
    else:
        print('   falta "%s": el indice de pobreza municipal de 2018 es una entrada externa '
              'que no se deriva de los microdatos vitales; sin el no se arma la base municipal.' % ipm)

    print('8. Notificaciones al sistema de vigilancia')
    if os.path.isdir(sivi):
        v = conteo_sivigila_550(sivi)
        pd.DataFrame([dict(anio=y, notificadas=v[y]['cordoba']) for y in sorted(v)]).to_csv(
            os.path.join(OUT, 'Notificaciones de mortalidad materna por anio.csv'), index=False)
    else:
        print('   carpeta SIVIGILA no disponible; se omite este bloque')

    print('9. Geometria municipal')
    if os.path.exists(geo):
        import shutil
        shutil.copy(geo, os.path.join(OUT, 'cordoba_municipios.geojson'))
    else:
        print('   geojson no disponible; se omite')

    print('\nBases derivadas escritas en "%s".' % OUT)
    print('Ahora ejecute "pipeline.py" para reproducir las cifras del manuscrito.')
    return 0


# ---------------------------------------------------------------- ejecucion
if __name__ == '__main__':
    ruta = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('BASE_MORTALIDAD', '')
    if not ruta:
        print('Indique la carpeta de los microdatos:')
        print('    python "Codigo reproducible mortalidad.py" RUTA_DE_LOS_MICRODATOS')
        print('\nCada funcion de este archivo reproduce un bloque del manuscrito:')
        print('  serie_anual            -> Tabla 4, Figura 2 y Tabla S1')
        print('  tendencia              -> tendencia global y por tipo de causa')
        print('  tendencia_segmentada   -> Tabla S2')
        print('  espacial_y_ecologico   -> Tablas S3 y S4, Figuras 4 y 5')
        print('  potencia_ecologica     -> filas de potencia de la Tabla S4')
        print('  flujo_referencia       -> Figura 6 y sensibilidad por ocurrencia')
        print('  conteo_sivigila_550    -> Tabla 2, verificacion de 242 frente a 247')
        sys.exit(2)
    sys.exit(reconstruir(ruta))
