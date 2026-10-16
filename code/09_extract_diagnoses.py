# INSPIRE arthroplasty anesthesia analysis (v1.02) - 09_extract_diagnoses.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import pandas as pd, numpy as np, os
a=pd.read_pickle(os.path.join(WORK,"arth.pkl"))
subs=set(a.subject_id)
dx=[]
for ch in pd.read_csv(os.path.join(DATA,"diagnosis.csv.gz"),chunksize=2_000_000):
    dx.append(ch[ch.subject_id.isin(subs)])
dx=pd.concat(dx); dx.to_pickle(os.path.join(WORK,"dx.pkl"))
print("rows",len(dx),"subjects with any dx",dx.subject_id.nunique(),"of",len(subs))
print("code examples:",dx.icd10_cm.astype(str).str.len().value_counts().head().to_dict(), dx.icd10_cm.astype(str).sample(8,random_state=1).tolist())
m=dx.merge(a[["op_id","subject_id","orin_time","discharge_time"]],on="subject_id")
m["dd"]=(m.chart_time-m.orin_time)/1440
print("dx timing rel to op (days) quantiles",m.dd.quantile([.05,.25,.5,.75,.95]).round(0).tolist())
print("ops with any pre-op dx:",m[m.dd<0].op_id.nunique(),"any dx during index stay after op:",m[(m.dd>0)&(m.chart_time<=m.discharge_time)].op_id.nunique(),"any dx within 30d after:",m[(m.dd>0)&(m.dd<=30)].op_id.nunique())
# distinct chart_times per subject: are diagnoses timestamped at admission?
print("frac dx at exactly admission-ish? n unique times per subject median:",dx.groupby("subject_id").chart_time.nunique().median())
