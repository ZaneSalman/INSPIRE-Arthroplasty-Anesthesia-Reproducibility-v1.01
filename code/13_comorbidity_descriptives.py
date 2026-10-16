# INSPIRE arthroplasty anesthesia analysis (v1.02) - 13_comorbidity_descriptives.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os, math
exec(open(os.path.join(HERE,"12_comorbidity_adjusted_models.py")).read().split("# g-computation")[0].replace("print(","(lambda *x,**y:None)("))
for g,l in [(0,"GA"),(1,"NA")]:
    s=a[a.neur==g]; n=len(s)
    print(l,"CCI mean",round(s.cci.mean(),2),"sd",round(s.cci.std(),2),"cci>=2",int(s.cci2.sum()),round(100*s.cci2.mean(),1))
cc=pd.read_pickle(os.path.join(WORK,"cci.pkl")); x=a[["op_id","neur"]].merge(cc,on="op_id")
for k in ["dm","renal","cvd","cancer","chf","pulm"]:
    r=[]
    for g in [0,1]:
        s=x[x.neur==g]; r.append(f"{int(s[k].sum())}/{len(s)} ({100*s[k].mean():.1f})")
    p0=x[x.neur==0][k].mean(); p1=x[x.neur==1][k].mean(); smd=(p0-p1)/math.sqrt((p0*(1-p0)+p1*(1-p1))/2)
    print(k,r,round(smd,3))
# IPTW2 balance for each covariate
for v in ["age","female","asa3","emop","knee","cci1","cci2"]: print("bal",v,round(smd_w(v,"iptw2"),3))
# p values for forest: already in fmt. E-values for CCI-adjusted
def ev(rr): rr=1/rr if rr<1 else rr; return rr+math.sqrt(rr*(rr-1))
for nm,o,hi in [("ICU",0.518,1.013),("AKI",0.691,1.242),("COMP",0.569,0.917)]: print("E cci",nm,round(ev(o),2), (round(ev(hi),2) if hi<1 else 1.0))
# coded AKI N17 validation already: print lab AKI vs coded
