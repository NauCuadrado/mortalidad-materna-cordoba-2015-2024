# -*- coding: utf-8 -*-
"""Por que el suavizado espacial no entra en la inferencia.
Naudith C, agosto de 2026.

Genero conteos municipales sin ninguna estructura espacial y le paso el indice de
Moran a la tasa cruda y a la tasa suavizada con el estimador espacial empirico de
Bayes. Si el procedimiento fuera valido, con datos sin estructura deberia rechazar
alrededor del 5 por ciento de las veces.
"""
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import libpysal
import esda

BASE = os.path.dirname(os.path.abspath(__file__))
if not os.path.isdir(os.path.join(BASE, 'Datos derivados')):
    BASE = os.path.dirname(BASE)
DATA = os.path.join(BASE, 'Datos derivados')
OUT = os.path.join(BASE, 'Salidas')

M = pd.read_csv(os.path.join(DATA, 'Base municipal con suavizado y LISA.csv'))
g = gpd.read_file(os.path.join(DATA, 'municipios de Cordoba.geojson'))
g['cod'] = g['cod'].astype(str)
M['cod_municipio'] = M.cod_municipio.astype(str)
gdf = (g.merge(M, left_on='cod', right_on='cod_municipio', suffixes=('_geo', ''))
       .sort_values('cod_municipio').reset_index(drop=True))
w = libpysal.weights.Queen.from_dataframe(gdf, use_index=False)
w.transform = 'r'
w.id_order = list(range(len(gdf)))

nv = gdf.nacidos_vivos.values.astype(float)
tasa = M.muertes_tempranas.sum() / nv.sum()
REP = 1000

rng = np.random.default_rng(42)
rechazos = {'Tasa cruda': 0, 'Suavizado espacial': 0}
for _ in range(REP):
    y = rng.poisson(nv * tasa)
    np.random.seed(42)
    if esda.Moran(y / nv * 1e5, w, permutations=999).p_sim < 0.05:
        rechazos['Tasa cruda'] += 1
    s = esda.smoothing.Spatial_Empirical_Bayes(y, nv, w)
    np.random.seed(42)
    if esda.Moran(np.nan_to_num(s.r), w, permutations=999).p_sim < 0.05:
        rechazos['Suavizado espacial'] += 1

filas = [dict(entrada=k, replicas=REP, rechazos=v, falsos_positivos=v / REP)
         for k, v in rechazos.items()]
pd.DataFrame(filas).to_csv(os.path.join(OUT, 'Simulacion del suavizado espacial.csv'), index=False)
for f in filas:
    print('%-20s rechazos %4d de %d  (%.1f%%)'
          % (f['entrada'], f['rechazos'], REP, f['falsos_positivos'] * 100))
