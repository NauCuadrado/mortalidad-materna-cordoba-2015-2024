# -*- coding: utf-8 -*-
"""Genera las figuras 1 a 6 y la figura S1 del articulo de mortalidad materna.
   Lee unicamente las bases de la carpeta "Datos derivados" y escribe PNG y PDF
   en la carpeta "Figuras". Uso: python "Figuras del articulo.py"."""
import os
import numpy as np, pandas as pd, geopandas as gpd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import Normalize
from scipy.stats import chi2
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.edgecolor':'#555555'})
BASE = os.path.dirname(os.path.abspath(__file__))
if not os.path.isdir(os.path.join(BASE, 'Datos derivados')):
    # permite ejecutar el archivo desde una subcarpeta como "Codigo"
    BASE = os.path.dirname(BASE)
DATA = os.path.join(BASE, 'Datos derivados')
OUT  = os.path.join(BASE, 'Figuras') + os.sep
os.makedirs(OUT, exist_ok=True)
def _d(nombre): return os.path.join(DATA, nombre)
def _gdf():
    g = gpd.read_file(_d('cordoba_municipios.geojson')); g['cod'] = g['cod'].astype(str)
    M = pd.read_csv(_d('Base municipal con suavizado y LISA.csv')); M['cod_municipio'] = M.cod_municipio.astype(str)
    return (g.merge(M, left_on='cod', right_on='cod_municipio', suffixes=('_geo',''))
             .sort_values('cod_municipio').reset_index(drop=True))
BLUE='#1F4E79'; RED='#C00000'; GREY='#666666'
SRC='Source: DANE vital statistics.'
def save(fig,name):
    fig.savefig(OUT+name+'.png',dpi=300,bbox_inches='tight')
    fig.savefig(OUT+name+'.pdf',bbox_inches='tight'); plt.close(fig)

# --------------------------------------------------- Fig 1 flujo
fig,ax=plt.subplots(figsize=(9.8,7.4)); ax.set_xlim(0,100); ax.set_ylim(2,100); ax.axis('off')
def box(c,y,txt,bold=False,fc='white',w=34,h=8):
    ax.add_patch(FancyBboxPatch((c-w/2,y),w,h,boxstyle='round,pad=0.35,rounding_size=0.8',fc=fc,ec=BLUE,lw=1.0))
    ax.text(c,y+h/2,txt,ha='center',va='center',fontsize=8.4,fontweight='bold' if bold else 'normal',color='#111111')
def step(c,y1,y2,txt=None):
    xa=c-11
    ax.add_patch(FancyArrowPatch((xa,y1),(xa,y2),arrowstyle='-|>',mutation_scale=11,color=GREY,lw=1.0))
    if txt:
        x=c-9; w=26; h=5; y=(y1+y2)/2-h/2
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.25,rounding_size=0.5',fc='#F4F4F4',ec='#B0B0B0',lw=0.7))
        ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=7.3,color='#333333')
CA,CB=27,75
ax.text(CA,95,'Vital statistics (DANE)',ha='center',fontsize=9.8,fontweight='bold',color=BLUE)
ax.text(CB,95,'Surveillance notifications of maternal death',ha='center',fontsize=9.8,fontweight='bold',color=BLUE)
LV=[85,68,51,34]
box(CA,LV[0],'Non-fetal deaths recorded in\nColombia, 2015 to 2024\nn = 2 646 843')
step(CA,LV[0],LV[1]+8,'Excluded: men and women not\nresident in Córdoba\nn = 2 610 957')
box(CA,LV[1],'Deaths of women resident\nin Córdoba\nn = 35 886')
step(CA,LV[1],LV[2]+8,'Excluded: underlying cause outside\nO00 to O99 and A34\nn = 35 644')
box(CA,LV[2],'Maternal deaths by ICD-10 code\nn = 242')
step(CA,LV[2],LV[3]+8,'Late maternal deaths, O96 (n = 78),\nanalyzed separately; O97 (n = 2) excluded')
box(CA,LV[3],'Early maternal deaths, numerator\nof the mortality ratio\nn = 162',bold=True,fc='#EAF1F8')
box(CB,LV[0],'National notifications of\nmaternal death, 2015 to 2024\nn = 5 177')
step(CB,LV[0],LV[1]+8,'Excluded: residents outside\nCórdoba, n = 4 930')
box(CB,LV[1],'Notifications of residents\nof Córdoba\nn = 247')
step(CB,LV[1],LV[2]+8,'Discarded or reclassified: 0')
box(CB,LV[2],'Notified maternal deaths\nn = 247')
step(CB,LV[2],LV[3]+8)
box(CB,LV[3],'Aggregate counts under non-equivalent\ndefinitions: 242 and 247',bold=True,fc='#EAF1F8')
step(CA,LV[3],22)
box(CA,14,'Maternal mortality ratio\n162 / 243 779 live births',bold=True,fc='#EAF1F8')
ax.text(0.5,7,'Deaths are attributed by department and municipality of residence, consistent with the live-birth denominator. Surveillance notifications are counted\nunder a different case definition and are reported for comparison of totals only. Sources: national vital statistics and national surveillance microdata.',
        fontsize=7.2,color='#666666',va='top')
save(fig,'Figure S2')

# --------------------------------------------------- Fig 2 razon anual
T=pd.read_csv(_d('Serie anual con intervalos.csv')); R=pd.read_csv(_d('Serie anual numeradores y denominadores.csv'))
R['d_resto']=R.d_nac33-R.d_cor; R['b_resto']=R.b_nac33-R.b_cor
def ic(k,n,m=1e5):
    lo=chi2.ppf(0.025,2*k)/2 if k>0 else 0; hi=chi2.ppf(0.975,2*(k+1))/2
    return k/n*m, lo/n*m, hi/n*m
rest=[ic(k,n) for k,n in zip(R.d_resto,R.b_resto)]
yrs=T.anio.values
fig,ax=plt.subplots(figsize=(7.6,4.6))
ax.axvspan(2019.5,2021.5,color='#EDEDED',zorder=0)
ax.fill_between(yrs,[x[1] for x in rest],[x[2] for x in rest],color=BLUE,alpha=0.15,lw=0)
ax.plot(yrs,[x[0] for x in rest],'-o',color=BLUE,lw=1.9,ms=4.5,label='Rest of the country')
ax.fill_between(yrs,T.rmm_cor_lo,T.rmm_cor_hi,color=RED,alpha=0.15,lw=0)
ax.plot(yrs,T.rmm_cor,'-s',color=RED,lw=1.9,ms=4.5,label='Córdoba')
for x,y,n in zip(yrs,T.rmm_cor,T.d_cor): ax.annotate(str(int(n)),(x,y),fontsize=7,color=RED,xytext=(0,7),textcoords='offset points',ha='center')
ax.set_xticks(yrs); ax.set_xlabel('Year'); ax.set_ylabel('Deaths per 100 000 live births')
ax.set_ylim(0,160); ax.grid(alpha=0.2,lw=0.6); ax.set_axisbelow(True); ax.legend(fontsize=8.6,frameon=False)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.text(2020.5,152,'Pandemic',ha='center',fontsize=8,color='#777777')
fig.text(0.012,0.022,'Shaded bands are exact 95% Poisson confidence intervals. Annual counts of deaths in the department are shown above each point.',fontsize=7.2,color='#777777')
fig.text(0.012,0.004,SRC,fontsize=7.2,color='#777777')
fig.tight_layout(rect=[0,0.05,1,1]); save(fig,'Figure 1')

# --------------------------------------------------- Fig 3 causas
E=pd.read_csv(_d('Causas de las 162 muertes.csv')); E['tipo']=np.where(E.grupo=='Indirect causes','Indirect',np.where(E.grupo=='Undetermined causes','Undetermined','Direct'))
B=R.set_index('anio').b_cor; yl=list(range(2015,2025))
dirs=[((E.tipo=='Direct')&(E.anio==y)).sum() for y in yl]; inds=[((E.tipo=='Indirect')&(E.anio==y)).sum() for y in yl]
rd=[d/B[y]*1e5 for d,y in zip(dirs,yl)]; ri=[i/B[y]*1e5 for i,y in zip(inds,yl)]
fig,axes=plt.subplots(1,2,figsize=(11.4,4.9),gridspec_kw={'width_ratios':[1.45,1]})
ax=axes[0]; ax.axvspan(2019.5,2021.5,color='#EDEDED',zorder=0)
ax.plot(yl,rd,'-o',color=BLUE,lw=2,ms=5,label='Direct obstetric causes')
ax.plot(yl,ri,'-s',color=RED,lw=2,ms=5,label='Indirect causes')
ax.set_xticks(yl); ax.set_xlabel('Year'); ax.set_ylabel('Deaths per 100 000 live births'); ax.set_ylim(-2,84)
ax.set_title('A. Annual rates by cause group',fontsize=10.5,fontweight='bold',loc='left',pad=6)
ax.text(2020.5,80,'Pandemic',ha='center',fontsize=8,color='#777777')
ax.grid(alpha=0.2,lw=0.6); ax.set_axisbelow(True); ax.legend(fontsize=8.5,frameon=False,loc='upper left')
for s in ('top','right'): ax.spines[s].set_visible(False)
for x,y,n in zip(yl,rd,dirs): ax.annotate(str(n),(x,y),fontsize=7,color=BLUE,xytext=(0,-11),textcoords='offset points',ha='center')
for x,y,n in zip(yl,ri,inds): ax.annotate(str(n),(x,y),fontsize=7,color=RED,xytext=(0,7),textcoords='offset points',ha='center')
ax=axes[1]; n1=int(B.loc[2015:2021].sum()); n2=int(B.loc[2022:2024].sum())
G=['Indirect causes','Hypertensive disorders','Obstetric hemorrhage','Abortive outcomes']
LB=['Indirect causes','Hypertensive\ndisorders','Obstetric\nhemorrhage','Abortive\noutcomes']
ypos=np.arange(len(G))[::-1]
for i,(g,yp) in enumerate(zip(G,ypos)):
    k1=int(((E.grupo==g)&(E.anio<=2021)).sum()); k2=int(((E.grupo==g)&(E.anio>=2022)).sum())
    r1=ic(k1,n1); r2=ic(k2,n2)
    ax.plot([r1[1],r1[2]],[yp+0.17]*2,color=BLUE,lw=1.4)
    ax.plot(r1[0],yp+0.17,'o',color=BLUE,ms=6,label='2015 to 2021' if i==0 else None)
    ax.plot([r2[1],r2[2]],[yp-0.17]*2,color=RED,lw=1.4)
    ax.plot(r2[0],yp-0.17,'s',color=RED,ms=6,label='2022 to 2024' if i==0 else None)
ax.set_yticks(ypos); ax.set_yticklabels(LB,fontsize=8.6); ax.set_ylim(-0.65,len(G)-0.35); ax.set_xlim(0,46)
ax.set_xlabel('Deaths per 100 000 live births (95% CI)')
ax.set_title('B. Cause-specific rates by period',fontsize=10.5,fontweight='bold',loc='left',pad=6)
ax.grid(axis='x',alpha=0.2,lw=0.6); ax.set_axisbelow(True); ax.legend(fontsize=8.5,frameon=False,loc='lower right')
for s in ('top','right','left'): ax.spines[s].set_visible(False)
fig.text(0.012,0.012,'Counts of deaths are shown next to each point in panel A. '+SRC,fontsize=7.2,color='#777777')
fig.tight_layout(rect=[0,0.035,1,1]); save(fig,'Figure 2')

# --------------------------------------------------- Fig 4 mapas
gdf=_gdf()
vmax=max(gdf.RMM_cruda.max(),gdf.RMM_EBglobal.max()); norm=Normalize(0,vmax)
fig,axes=plt.subplots(1,2,figsize=(11.6,6.4))
for ax,col,tit in [(axes[0],'RMM_cruda','A. Crude ratio'),(axes[1],'RMM_EBglobal','B. Empirical Bayes smoothed ratio')]:
    gdf.plot(column=col,cmap='OrRd',norm=norm,edgecolor='#666666',lw=0.4,ax=ax)
    ax.set_title(tit,fontsize=10.5,fontweight='bold',loc='left',color='#000000',pad=6); ax.axis('off')
    for _,r in gdf.iterrows():
        c=r.geometry.representative_point(); ax.annotate(r.municipio,(c.x,c.y),fontsize=5.4,ha='center',color='#111111')
sm_=plt.cm.ScalarMappable(cmap='OrRd',norm=norm); sm_.set_array([])
cb=fig.colorbar(sm_,ax=axes,orientation='horizontal',fraction=0.045,pad=0.04,aspect=45)
cb.set_label('Maternal mortality ratio per 100 000 live births (same scale in both panels)',fontsize=8.5)
fig.suptitle('Municipal maternal mortality, Córdoba, 2015 to 2024',fontsize=12,fontweight='bold',color='#000000',x=0.012,ha='left',y=0.97)
fig.text(0.012,0.015,'Global empirical Bayes smoothing. Crude rates are the input for the spatial autocorrelation tests. '+SRC,fontsize=7.2,color='#777777')
fig.savefig(OUT+'Figure S3.png',dpi=300,bbox_inches='tight')
fig.savefig(OUT+'Figure S3.pdf',bbox_inches='tight'); plt.close(fig)

# --------------------------------------------------- Fig 5 funnel
theta=gdf.muertes_tempranas.sum()/gdf.nacidos_vivos.sum()*1e5
n=np.linspace(1000,gdf.nacidos_vivos.max()*1.1,800)
up95=theta+1.959964*np.sqrt(theta*1e5/n); lo95=np.maximum(theta-1.959964*np.sqrt(theta*1e5/n),0)
up998=theta+3.090232*np.sqrt(theta*1e5/n); lo998=np.maximum(theta-3.090232*np.sqrt(theta*1e5/n),0)
g=gdf.copy(); g['u']=np.interp(g.nacidos_vivos,n,up95); g['l']=np.interp(g.nacidos_vivos,n,lo95)
g['out']=(g.RMM_cruda>g.u)|(g.RMM_cruda<g.l)
fig,ax=plt.subplots(figsize=(7.6,5.2))
ax.plot(n,lo95,color='#777777',lw=1,ls='--'); ax.plot(n,up95,color='#777777',lw=1,ls='--',label='95% limits')
ax.plot(n,lo998,color='#BBBBBB',lw=1,ls=':'); ax.plot(n,up998,color='#BBBBBB',lw=1,ls=':',label='99.8% limits')
ax.axhline(theta,color=BLUE,lw=1.6,label='Departmental ratio (%.1f)'%theta)
ax.scatter(g.nacidos_vivos[~g.out],g.RMM_cruda[~g.out],s=34,color=BLUE,zorder=5,edgecolor='white',lw=0.6)
ax.scatter(g.nacidos_vivos[g.out],g.RMM_cruda[g.out],s=52,color=RED,zorder=6,edgecolor='white',lw=0.6)
offs={'Monteria':(-38,8),'Lorica':(8,-2),'Tierralta':(8,-10),'La Apartada':(9,2),'Moñitos':(9,2)}
for _,r in g[g.out | (g.nacidos_vivos>=16000)].iterrows():
    dx,dy=offs.get(r.municipio,(8,3))
    ax.annotate(r.municipio,(r.nacidos_vivos,r.RMM_cruda),fontsize=7.6,color='#333333',xytext=(dx,dy),textcoords='offset points')
ax.set_xscale('log'); ax.set_xlabel('Live births, 2015 to 2024 (log scale)')
ax.set_ylabel('Maternal mortality ratio per 100 000 live births')
ax.set_ylim(-12,300); ax.grid(alpha=0.2,lw=0.6); ax.set_axisbelow(True)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.legend(fontsize=8,frameon=False,loc='upper right')
ax.text(0,-0.165,SRC,transform=ax.transAxes,fontsize=7.2,color='#777777')
fig.tight_layout(rect=[0,0.035,1,1]); save(fig,'Figure S4')

# --------------------------------------------------- Fig 6 flujo
fig,axes=plt.subplots(1,2,figsize=(11.6,4.8),gridspec_kw={'width_ratios':[1,1.05]})
ax=axes[0]
lab=['Montería','Lorica','Other municipalities\nof Córdoba','Outside Córdoba']
val=[121,17,12,12]; y=np.arange(len(lab))[::-1]
ax.barh(y,val,color=[BLUE,BLUE,'#9EC3E0',RED],height=0.6)
for yy,v in zip(y,val): ax.text(v+2,yy,str(v),va='center',fontsize=9,color='#333333')
ax.set_yticks(y); ax.set_yticklabels(lab,fontsize=9)
ax.set_xlabel('Deaths of women resident in Córdoba'); ax.set_xlim(0,140)
ax.set_title('A. Municipality of occurrence of the 162 deaths',fontsize=10.5,fontweight='bold',loc='left',pad=6)
ax.grid(axis='x',alpha=0.2,lw=0.6); ax.set_axisbelow(True)
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.text(0.0,-0.32,'Montería accounts for 81% of the deaths that occurred inside the department,\nand 48 of those 121 women were resident there.',transform=ax.transAxes,fontsize=7.8,color='#555555')
ax=axes[1]
ins=[('Antioquia',32),('Sucre',6),('Bolívar',4),('Chocó',2),('Other',3)]
outs=[('Sucre',6),('Atlántico',3),('Antioquia',2),('Bogotá',1)]
yin=np.arange(len(ins))[::-1]+1.6
ax.barh(yin,[v for _,v in ins],color=BLUE,height=0.55)
for yy,(nm,v) in zip(yin,ins): ax.text(v+0.8,yy,str(v),va='center',fontsize=8.6,color='#333333')
yout=np.arange(len(outs))[::-1]-len(outs)+0.4
ax.barh(yout,[-v for _,v in outs],color=RED,height=0.55)
for yy,(nm,v) in zip(yout,outs): ax.text(-v-0.8,yy,str(v),va='center',ha='right',fontsize=8.6,color='#333333')
ax.set_yticks(list(yin)+list(yout)); ax.set_yticklabels([n for n,_ in ins]+[n for n,_ in outs],fontsize=8.2)
ax.axvline(0,color='#888888',lw=1); ax.set_xlim(-16,40); ax.set_ylim(min(yout)-1.0,max(yin)+1.5)
ax.set_xticks([-10,-5,0,10,20,30]); ax.set_xticklabels(['10','5','0','10','20','30'])
ax.set_xlabel('Number of maternal deaths')
ax.set_title('B. Deaths crossing departmental boundaries: 47 in, 12 out',fontsize=10.5,fontweight='bold',loc='left',pad=6)
ax.grid(axis='x',alpha=0.2,lw=0.6); ax.set_axisbelow(True)
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.text(-15.4,max(yin)+0.9,'Occurred in Córdoba, resident elsewhere',fontsize=8.2,color=BLUE,fontweight='bold')
ax.text(-15.4,min(yout)-0.75,'Resident in Córdoba, occurred elsewhere',fontsize=8.2,color=RED,fontweight='bold')
fig.text(0.012,0.012,'Deaths are attributed by department of residence throughout the analysis. '+SRC,fontsize=7.2,color='#777777')
fig.tight_layout(rect=[0,0.07,1,1]); save(fig,'Figure S5')
# --------------------------------------------------- Fig S1 orden de los territorios
K = pd.read_csv(_d('Ranking departamental RMM.csv')).sort_values('RMM', ascending=True).reset_index(drop=True)
lo, hi = [], []
for k, n in zip(K.muertes, K.nacidos):
    a = ic(k, n); lo.append(a[1]); hi.append(a[2])
col = [RED if 'rdoba' in str(x) else BLUE for x in K.Departamento]
fig, ax = plt.subplots(figsize=(7.4, 8.6))
y = np.arange(len(K))
for i in range(len(K)):
    ax.plot([lo[i], hi[i]], [y[i]] * 2, color=col[i], lw=1.1, alpha=0.55)
ax.scatter(K.RMM, y, s=26, color=col, zorder=3)
ax.set_yticks(y); ax.set_yticklabels(K.Departamento, fontsize=7.6)
for t, c in zip(ax.get_yticklabels(), col):
    if c == RED: t.set_color(RED); t.set_fontweight('bold')
ax.set_xlabel('Deaths per 100 000 live births, 2015 to 2024 (95% CI)')
ax.set_xlim(0, 340); ax.grid(axis='x', alpha=0.2, lw=0.6); ax.set_axisbelow(True)
for s_ in ('top', 'right', 'left'): ax.spines[s_].set_visible(False)
fig.text(0.012, 0.008, 'Exact Poisson confidence intervals. The department under study is highlighted. ' + SRC,
         fontsize=7.2, color='#777777')
fig.tight_layout(rect=[0, 0.025, 1, 1]); save(fig, 'Figure S1')

print('figuras listas')
