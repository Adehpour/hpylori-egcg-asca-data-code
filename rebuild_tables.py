"""Rebuild all computable contrasts and all tabular figure inputs from the bundled CSVs.
Run from any working directory. Outputs are written under rebuilt/.
Per-strain host comparisons use the five-group within-strain error. Only
computable tests based on experimental replicate data are included; no
approximate P values are reconstructed from rounded manuscript summaries.
"""
from pathlib import Path
import csv, json, hashlib, itertools as it
import numpy as np
import pandas as pd
from scipy import stats as st
from collections import Counter
B=Path(__file__).resolve().parent
raw=B/'archive_csv'; sub=raw/'assay_exports'; repro=B/'rebuilt';repro.mkdir(exist_ok=True)
strains=['26695','CLR-resistant']
groups=['Untreated','Blank NP','Free EGCG','EGCG-NP','Clarithromycin']
genes=['cagA','vacA','ureA','flaA','babA','sabA']
def csvr(path):
 with open(path,encoding='utf-8-sig',newline='') as f: return list(csv.reader(f))
def loadwide(path):
 rs=csvr(path); start=7 if rs[0][1]=='No bacteria' else 1
 return np.array([[float(v) for v in row[start:start+30]] for row in rs[2:] if row[0] in strains]).reshape(2,5,6)
def find(path,pattern):return next(path.glob(pattern))
sourcefiles={}
endpoints={}
def store(key,base,glob):
 f=find(base,glob); sourcefiles[key]=str(f.relative_to(raw));endpoints[key]=loadwide(f)
store('Urease_pct',sub,'*Residual urease*');store('Biofilm_pct',sub,'*Biofilm biomass*')
for g in genes:store(g+'_dCt',sub,'*model*'+g+'.csv')
store('Associated_log10_CFU_well',raw,'Fig3B*')
# The laboratory-confirmed CagA table is the canonical source. Obsolete exports
# that contained an incorrect third blank-NP ratio are not distributed.
cagfile=sub/'phospho-CagA _ total CagA.csv'
endpoints['Phospho_total_CagA']=np.array([r[1:] for r in csvr(cagfile)[1:]],float).reshape(5,2,6).transpose(1,0,2)
sourcefiles['Phospho_total_CagA']=str(cagfile.relative_to(raw))
store('Nuc_cyto_p65',raw,'Fig4A*');store('IL8_pg_mL',raw,'Fig4B*')
rows=[]
for s,t,j in it.product(range(2),range(5),range(6)):
 rows.append({'strain':strains[s],'treatment':groups[t],'replicate':j+1,**{k:v[s,t,j] for k,v in endpoints.items()}})
pd.DataFrame(rows).to_csv(repro/'asca_input.csv',index=False)
print('Same-layout crosschecks')
for k,nm in [('Associated_log10_CFU_well','AGS_adhesion.csv'),('Phospho_total_CagA','phospho-CagA _ total CagA.csv'),('Nuc_cyto_p65','Nuclear _ cytosolic p65 ratio.csv')]:
 rs=csvr(sub/nm); x=np.array([[float(v) for v in row[1:13]] for row in rs[1:6]]).reshape(5,2,6).transpose(1,0,2);print(k,abs(x-endpoints[k]).max())
for s,sn in enumerate(strains):
 rs=csvr(find(sub,'*biological mean DeltaCt _ '+sn+'.csv'))
 for rr in rs[2:]:
  v=np.array(rr[1:],float).reshape(5,6);g=rr[0];assert np.allclose(v,endpoints[g+'_dCt'][s],atol=2e-5)
print('batch')
features=['Size_nm','PDI','Zeta_mV','EE_pct','DL_pct']
bdata=[np.array(csvr(find(raw,pat))[1:],float).T for pat in ['*Particle size*','*Polydispersity*','*Zeta potential*','*Encapsulation*','*Drug loading*']]
brows=[]
for t,j in it.product(range(5),range(3)):
 brows.append({'formulation':'F'+str(t+1),'batch_ordinal':j+1,**{k:float(x[t,j]) for k,x in zip(features,bdata)}})
pd.DataFrame(brows).to_csv(repro/'formulation_batches.csv',index=False)
print(pd.DataFrame(brows).groupby('formulation').agg(['mean','std']).round(3))
times=[0,2,4,8,12,24]
tk=[]
for sn in strains:
 rs=csvr(find(sub,'*Time-kill*'+sn+'.csv'))
 assert [int(a[0]) for a in rs[2:]]==times
 tk.append(np.array([a[1:] for a in rs[2:]],float).reshape(6,5,6))
tk=np.array(tk)
time_rows=[];kr=[]
for s,h,t,j in it.product(range(2),range(6),range(5),range(6)):
 l0=tk[s,0,t,j];lt=tk[s,h,t,j]; diff=l0-lt;red=100*(1-10**(-diff))
 time_rows.append([strains[s],groups[t],j+1,times[h],l0,lt,float(diff),float(red)])
 if h==5: kr.append(time_rows[-1])
tdf=pd.DataFrame(time_rows,columns=['Strain','Treatment','Run','Time_h','Baseline_log10','Observed_log10','Log_reduction','Net_reduction_pct'])
tdf.to_csv(repro/'time_kill.csv',index=False)
kst=tdf[tdf.Time_h==24].groupby(['Strain','Treatment'],sort=False).Net_reduction_pct.agg(['mean','std']).reset_index();print('24h',kst.round(5).to_string(index=False))
mrs=[];bound=[]
for sn in strains:
 rs=csvr(find(raw,'Fig2A*'+sn+'.csv')); mat=np.array([a[1:] for a in rs[2:]],float).reshape(2,4,4)
 for ti,grp in enumerate(groups[1:]):
  for j in range(4):mrs.append([sn,grp,j+1,mat[0,ti,j],mat[1,ti,j],bool(mat[0,ti,j]==mat[1,ti,j])])
  bound.append([sn,grp,int(sum(mat[0,ti]==mat[1,ti])),4])
mdf=pd.DataFrame(mrs,columns=['Strain','Treatment','Culture','MIC','MBC','MIC_equals_MBC'])
mdf.to_csv(repro/'mic_mbc.csv',index=False);print('MBC bound',bound)
powerdata=np.array(list(endpoints.values()))
def holm(p):
 p=np.asarray(p); ii=np.argsort(p);o=p.copy();o[ii]=np.minimum(1.,np.maximum.accumulate(p[ii]*np.arange(p.size,0,-1)));return o
def tuk(q,k,df):return float(st.studentized_range.sf(abs(q),k,df))
pairs=[]; rowdict={}
def pair(panel,endpoint,s,context,a,b,effect,scale,calc,adjust,family,origin):
 row=dict(Panel=panel,Endpoint=endpoint,Stratum=s,Context=context,Group_A=a,Group_B=b,Effect_A_minus_B=float(effect),Scale=scale,P_adjusted=float(calc),Method=adjust,Family=family,Source=origin)
 pairs.append(row);rowdict[panel,endpoint,s,context,a,b]=row
means={};sexp={}
for key,arr in endpoints.items():
 if key=='Phospho_total_CagA': tr=np.log(arr);scale='ln ratio'
 elif key=='IL8_pg_mL':tr=np.log10(arr);scale='log10 pg/mL'
 else:tr=arr.copy();scale={'Urease_pct':'percentage points','Biofilm_pct':'percentage points, normalized','Associated_log10_CFU_well':'log10 CFU/well','Nuc_cyto_p65':'ratio difference'}.get(key,'Ct cycles')
 sm=tr.mean(2);ms=((tr-sm[:,:,None])**2).sum()/50
 panel='Fig 4' if key.endswith('_dCt') else {'Urease_pct':'Fig 3C','Biofilm_pct':'Fig 3D','Associated_log10_CFU_well':'Fig 5A','Phospho_total_CagA':'Fig 5B','Nuc_cyto_p65':'Fig 6A','IL8_pg_mL':'Fig 6B'}[key]
 means[key]=[arr.mean(2).tolist(),arr.std(2,ddof=1).tolist()]
 for s in range(2):
  for i,j in it.combinations(range(5),2):
   eff=sm[s,j]-sm[s,i];host=key in ['Associated_log10_CFU_well','Phospho_total_CagA','Nuc_cyto_p65','IL8_pg_mL']
   local_ms=np.sum((tr[s]-sm[s,:,None])**2)/25 if host else ms
   local_df=25 if host else 50
   p=tuk(eff/np.sqrt(local_ms/6),5,local_df)
   pair(panel,key,strains[s],'Within strain',groups[j],groups[i],eff,scale,p,('Tukey; five groups within strain, df=25' if host else 'Tukey; pooled full factorial residual, df=50'),'10 pairs among 5 treatments, within strain and endpoint','CSV reconstruction')
 pp=holm(st.t.sf(np.abs(sm[1]-sm[0])/np.sqrt(2*ms/6),50)*2)
 for j,p in enumerate(pp):
  pair(panel,key,groups[j],'Between strains',strains[1],strains[0],sm[1,j]-sm[0,j],scale,p,'Holm; pooled full factorial residual, df=50','5 strain comparisons within endpoint','Exploratory CSV reanalysis')
 if key in ['Urease_pct','Biofilm_pct','cagA_dCt','ureA_dCt','babA_dCt','Phospho_total_CagA','Associated_log10_CFU_well','Nuc_cyto_p65','IL8_pg_mL']:
  print(key,'LF effects',sm[:,3]-sm[:,2], 'Tukey LF', [round(rowdict[panel,key,s,'Within strain',groups[3],groups[2]]['P_adjusted'],6) for s in strains])
  if key=='Phospho_total_CagA':print('CagA geom',np.exp(sm[:,[2,3]]))
 # Primary Holm families, made explicit in manuscript, distinct from all-pairs family
 if key in ('Urease_pct','Biofilm_pct'):
  con=[(s,j,i,sm[s,j]-sm[s,i]) for s in range(2) for j,i in [(3,2),(1,0)]]
  pp=holm([2*st.t.sf(abs(v)/np.sqrt(2*ms/6),50) for s,j,i,v in con])
  for (s,j,i,v),p in zip(con,pp):
   pair(panel,key,strains[s],'Planned / normalized source',groups[j],groups[i],v,scale,p,'Holm; pooled factorial residual, df=50','4 planned pairs across 2 strains','Urease: CSV; biofilm: normalized CSV sensitivity')
for s in range(2):
 vals=tk[s]; cmean=vals.mean(2);runbias=(vals-cmean[:,:,None]).mean((0,1))
 # Balance gives a common pooled error term for within-time, between-treatment contrasts.
 res=vals-cmean[:,:,None]-runbias[None,None,:];df=res.size-6*5-5;ms=np.sum(res**2)/df
 for h,i,j in it.product(range(6),range(5),range(5)):
  if i>=j:continue
  eff=cmean[h,j]-cmean[h,i];p=tuk(eff/np.sqrt(ms/6),5,df)
  pair('Fig 3A/B','Viable_log10_CFU_mL',strains[s],f'{times[h]} h',groups[j],groups[i],eff,'log10 CFU/mL',p,'Tukey; observational run-block model, df=145','10 pairs within strain and time','Exploratory CSV reanalysis')
 for h in [3,5]:print('TK',s,times[h],cmean[h,3]-cmean[h,2])
# Factor sweep on formulation batches; complete Tukey set, exact from supplied batch exports.
for key,mat in zip(features,bdata):
 av=mat.mean(1);ms=((mat-av[:,None])**2).sum()/10
 for i,j in it.combinations(range(5),2):
  eff=av[j]-av[i];p=tuk(eff/np.sqrt(ms/3),5,10)
  pair('Table 1A',key,'F1–F5','Between formulations','F'+str(j+1),'F'+str(i+1),eff,key,p,'Tukey; one-way residual, df=10','10 formulation pairs per outcome','CSV reconstruction')
# Gene family of six planned two-sided loaded-free effects.
gp=[]
for g in ('cagA','ureA','babA'):
 x=endpoints[g+'_dCt'];m=x.mean(2);err=((x-m[:,:,None])**2).sum()/50
 for s in range(2):
  eff=m[s,3]-m[s,2];gp.append((g,s,eff,st.t.sf(abs(eff)/np.sqrt(err/3),50)*2))
for (g,s,eff,p0),p in zip(gp,holm([v[3] for v in gp])):
 pair('Fig 4',g+'_dCt',strains[s],'Planned gene contrast',groups[3],groups[2],eff,'Ct cycles',p,'Holm; pooled factorial residual, df=50','6 pairs across 3 targets and 2 strains','CSV reconstruction')
print('gene LF',[(r['Endpoint'],r['Stratum'],round(r['P_adjusted'],8)) for r in pairs[-6:]])
# Reconstruct Dunnett control contrasts in each target/strain with the pooled model error.
from scipy.special import roots_hermitenorm, roots_genlaguerre, gamma, ndtr
z,w=roots_hermitenorm(64); z=z[:,None]; w=w/np.sqrt(2*np.pi)
u,v=roots_genlaguerre(60,24); v=v/gamma(25)
def dunnett4(t):
 width=np.sqrt(2)*abs(t)*np.sqrt(2*u/50)
 within=np.clip(ndtr(z+width)-ndtr(z-width),0,1)
 return float(w @ (1-within**4) @ v)
for g in genes:
 key=g+'_dCt';x=endpoints[key];av=x.mean(2);mse=((x-av[:,:,None])**2).sum()/50
 for si,j in it.product(range(2),range(1,5)):
  eff=av[si,j]-av[si,0]
  pair('Fig 4',key,strains[si],'Control comparison',groups[j],groups[0],eff,'Ct cycles',dunnett4(eff/np.sqrt(mse/3)),'Dunnett, two-sided, 4 vs common control; df=50','4 treatments versus untreated within target and strain','CSV reconstruction')
# Known p-values reported by the primary record are also transcribed on a separate worksheet.
for r in pairs:
 r['P_bound']='<' if r['P_adjusted'] < 1e-12 else '='
 if r['P_bound']=='<':r['P_adjusted']=1e-12
pd.DataFrame(pairs).to_csv(repro/'adjusted_contrasts.csv',index=False)
base={k:v.tolist() for k,v in endpoints.items()}
uninf={k:[float(v) for v in csvr(raw/sourcefiles[k])[2][1:7]] for k in ['Nuc_cyto_p65','IL8_pg_mL']}
files=[]
for p in sorted(raw.rglob('*.csv')):
 files.append([str(p.relative_to(raw)),hashlib.sha256(p.read_bytes()).hexdigest(),len(csvr(p))])
json.dump(dict(endpoints=base,means=means,uninfected=uninf,rows=rows,sourcefiles=sourcefiles,formulation=brows,formfeatures=features,pairs=pairs,time=tk.tolist(),time_rows=time_rows,mic=mrs,mbc_bound=bound,kill24=kst.to_dict('records'),files=files),open(repro/'figure_data.json','w'),ensure_ascii=False,allow_nan=False)

# Simultaneous Tukey intervals for the loaded/free CagA geometric mean ratio.
ci=[]
for si,sn in enumerate(strains):
 z=np.log(endpoints['Phospho_total_CagA'][si]);av=z.mean(1)
 mse=((z-av[:,None])**2).sum()/25
 half=st.studentized_range.ppf(.95,5,25)*np.sqrt(mse/6)
 diff=av[3]-av[2]
 ci.append([sn,np.exp(av[2]),np.exp(av[3]),np.exp(diff),np.exp(diff-half),np.exp(diff+half)])
pd.DataFrame(ci,columns=['strain','free_geometric_mean','loaded_geometric_mean','ratio','lower95','upper95']).to_csv(repro/'caga_ratio_intervals.csv',index=False)
