# INSPIRE arthroplasty anesthesia analysis (v1.02) - 07_corrected_icu_and_composite.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os, math
src=open(os.path.join(HERE,"04_aki_composite_sensitivity.py")).read().split("# propensity / IPTW on full cohort")[0]
exec(src)  # gives a (arthroplasty) with aki, comp(old), bmi, wlogit, fmt, XB, XB_nok, pnorm, Z
def auc(y,p):
    r=pd.Series(p).rank().values; n1=y.sum(); n0=len(y)-n1; return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def corr_icu(d):
    return (d.icuin_time.notna()&(d.icuin_time>=d.admission_time)&(d.icuin_time<=d.discharge_time)).astype(int)
a["icu_idx"]=corr_icu(a)
a["d_inh_idx"]=(a.inhosp_death_time.notna()&(a.inhosp_death_time>=a.admission_time)&(a.inhosp_death_time<=a.discharge_time+1440)).astype(int)
print("inhosp deaths total",int(a.d_inh.sum()),"within index (+1d):",int(a.d_inh_idx.sum()))
print("ICU_idx GA",int(a[a.neur==0].icu_idx.sum()),"/269 NA",int(a[a.neur==1].icu_idx.sum()),"/2271")
for g in [0,1]:
    s=a[a.neur==g]; print(" pct",round(100*s.icu_idx.mean(),1))
y=a.icu_idx.values.astype(float); X=XB(a)
b,se,p=wlogit(X,y); print("ICU_idx adjusted:"); 
for i,nm in enumerate(["const","neuraxial","age","female","asa3","emergency","knee"]):
    if i: print(f"   {nm}: {fmt(b,se,i)}")
print("   c-stat",round(auc(a.icu_idx.values,p),3))
# SMD for ICU_idx in table1
p0=a[a.neur==0].icu_idx.mean(); p1=a[a.neur==1].icu_idx.mean(); print("   SMD",round((p0-p1)/math.sqrt((p0*(1-p0)+p1*(1-p1))/2),3))
# IPTW
Xp=np.column_stack([np.ones(len(a)),a.age,a.female,a.asa3,a.emop,a.knee]).astype(float)
_,_,ps=wlogit(Xp,a.neur.values.astype(float)); prev=a.neur.mean()
iw=np.where(a.neur==1,prev/ps,(1-prev)/(1-ps)); a["iptw"]=np.clip(iw,None,np.quantile(iw,.99))
b,se,_=wlogit(X,y,a.iptw.values); print("ICU_idx IPTW:",fmt(b,se))
# first op per patient
f1=a.sort_values("orin_time").groupby("subject_id").head(1); b,se,_=wlogit(XB(f1),f1.icu_idx.values.astype(float)); print("ICU_idx first op n=",len(f1),fmt(b,se))
# BMI+cluster arthroplasty
ab=a.dropna(subset=["bmi"]); b,se,_=wlogit(XB(ab,True),ab.icu_idx.values.astype(float),groups=ab.subject_id); print("ICU_idx BMI+clust n=",len(ab),fmt(b,se))
# subgroups
for kn,l in [(0,"hip"),(1,"knee")]:
    s=a[a.knee==kn].dropna(subset=["bmi"]); ev=int(s.icu_idx.sum())
    try: b,se,_=wlogit(XB_nok(s),s.icu_idx.values.astype(float),groups=s.subject_id); print(f"   {l} n={len(s)} ev={ev}:",fmt(b,se))
    except Exception as e: print(l,"fail",e)
    print(f"   {l} events by group GA/NA:",int(s[s.neur==0].icu_idx.sum()),int(s[s.neur==1].icu_idx.sum()))
Xi=np.column_stack([X,a.neur*a.knee]); b,se,_=wlogit(Xi,y); print("   interaction p",round(pnorm(b[-1]/se[-1]),3))
# early (<=24h after OR exit) within index
a["icu_early"]=(a.icu_idx==1)&(((a.icuin_time-a.orout_time)/60)<=24)
print("early ICU GA/NA",int(a[a.neur==0].icu_early.sum()),int(a[a.neur==1].icu_early.sum()))
b,se,_=wlogit(X,a.icu_early.values.astype(float)); print("ICU early<=24h:",fmt(b,se))
# composite corrected
k=a.dropna(subset=["aki"]).copy(); k["comp_c"]=((k.aki==1)|(k.icu_idx==1)|(k.d_inh_idx==1)).astype(float)
for g,l in [(0,"GA"),(1,"NA")]:
    s=k[k.neur==g]; print(f"COMP {l}: {int(s.comp_c.sum())}/{len(s)} ({100*s.comp_c.mean():.1f})")
b,se,_=wlogit(XB(k),k.comp_c.values); print("COMP adj:",fmt(b,se))
b,se,_=wlogit(XB(k),k.comp_c.values,k.iptw.values); print("COMP IPTW:",fmt(b,se))
kb=k.dropna(subset=["bmi"]); b,se,_=wlogit(XB(kb,True),kb.comp_c.values,groups=kb.subject_id); print("COMP BMI+clust n=",len(kb),fmt(b,se))
k1=k.sort_values("orin_time").groupby("subject_id").head(1); b,se,_=wlogit(XB(k1),k1.comp_c.values); print("COMP first op n=",len(k1),fmt(b,se))
for kn,l in [(0,"hip"),(1,"knee")]:
    s=k[k.knee==kn].dropna(subset=["bmi"]); b,se,_=wlogit(XB_nok(s),s.comp_c.values,groups=s.subject_id); print(f"   COMP {l} n={len(s)} ev={int(s.comp_c.sum())}:",fmt(b,se))
Xi=np.column_stack([XB(k),k.neur*k.knee]); b,se,_=wlogit(Xi,k.comp_c.values); print("   COMP interaction p",round(pnorm(b[-1]/se[-1]),3))
# AKI unaffected; AKI IPTW reuse check
b,se,_=wlogit(XB(k),k.aki.values,k.iptw.values); print("AKI IPTW check:",fmt(b,se))
# all-orthopedic cohort
ops=pd.read_csv(os.path.join(DATA,"operations.csv.gz"))
o=ops[ops.department.eq("OS")&ops.antype.isin(["General","Neuraxial"])].copy()
o["los_raw"]=(o.discharge_time-o.admission_time)/1440; o=o[o.los_raw.gt(0)&o.los_raw.lt(365)].copy(); o["los"]=o.los_raw.clip(lower=.5)
o["neur"]=o.antype.eq("Neuraxial").astype(int); o["female"]=o.sex.eq("F").astype(int)
o["asa3"]=np.where(o.asa.isna(),np.nan,(o.asa>2).astype(float)); o["knee"]=o.icd10_pcs.astype(str).str[:4].isin({"0SRC","0SRD"}).astype(int)
code=o.icd10_pcs.astype(str)
o["proc"]=np.select([code.str[:4].isin({"0SRC","0SRD"}),code.str[:4].isin({"0SR9","0SRB","0SRA","0SRE"}),code.str[:3].eq("0SG")],["knee","hip","spine"],"other")
o=o.dropna(subset=["los","age","asa3","emop"]).copy()
o["icu_old"]=o.icuin_time.notna().astype(int); o["icu_idx"]=corr_icu(o)
o["bmi"]=o.weight/(o.height/100)**2; o.loc[~o.bmi.between(10,70),"bmi"]=np.nan
def XO(d,bmi=False):
    cols=[np.ones(len(d)),d.neur,d.age,d.female,d.asa3,d.emop,(d.proc=="knee").astype(int),(d.proc=="other").astype(int),(d.proc=="spine").astype(int)]+([d.bmi] if bmi else [])
    return np.column_stack(cols).astype(float)
print("allortho n",len(o))
ob=o.dropna(subset=["bmi"]); print("allortho bmi n",len(ob))
b,se,_=wlogit(XO(ob,True),ob.icu_old.values.astype(float),groups=ob.subject_id); print("VALID allortho ICU old (target 0.664 0.563-0.782):",fmt(b,se))
b,se,_=wlogit(XO(ob,True),ob.icu_idx.values.astype(float),groups=ob.subject_id); print("allortho ICU_idx:",fmt(b,se), "events",int(ob.icu_idx.sum()))
a.to_pickle(os.path.join(WORK,"arth2.pkl"))
