# INSPIRE arthroplasty anesthesia analysis (v1.02) - 02_extract_creatinine.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os
a=pd.read_pickle(os.path.join(WORK,"arth.pkl"))
print("admission_time==0 frac:",float((a.admission_time==0).mean()), "admission_time quantiles:",a.admission_time.quantile([.1,.5,.9]).tolist())
print("opdate quantiles:",a.opdate.quantile([.1,.5,.9]).tolist())
print("orin-admission days q:",((a.orin_time-a.admission_time)/1440).quantile([.1,.5,.9]).round(2).tolist())
subs=set(a.subject_id)
cr=[]
for ch in pd.read_csv(os.path.join(DATA,"labs.csv.gz"),chunksize=2_000_000):
    ch=ch[(ch.item_name=="creatinine")&ch.subject_id.isin(subs)]
    cr.append(ch)
cr=pd.concat(cr); cr.to_pickle(os.path.join(WORK,"cr.pkl"))
m=cr.merge(a[["op_id","subject_id","orin_time","admission_time","discharge_time"]],on="subject_id")
m["dt"]=(m.chart_time-m.orin_time)/1440
print("creatinine rows:",len(cr),"subjects w/ any:",cr.subject_id.nunique(),"of",len(subs))
print("dt (days from OR entry) quantiles:",m.dt.quantile([.05,.25,.5,.75,.95]).round(1).tolist())
print("frac ops w/ cr in [-30,0):",float(m[(m.dt>=-30)&(m.dt<0)].op_id.nunique()/len(a)),
      " in (0,7]:",float(m[(m.dt>0)&(m.dt<=7)].op_id.nunique()/len(a)),
      " in (0,2]:",float(m[(m.dt>0)&(m.dt<=2)].op_id.nunique()/len(a)))
print("deaths inhosp:",int(a.inhosp_death_time.notna().sum()),"allcause any:",int(a.allcause_death_time.notna().sum()),
      "allcause<=30d post:",int(((a.allcause_death_time-a.orin_time)/1440).between(0,30).sum()),
      "<=90d:",int(((a.allcause_death_time-a.orin_time)/1440).between(0,90).sum()))
