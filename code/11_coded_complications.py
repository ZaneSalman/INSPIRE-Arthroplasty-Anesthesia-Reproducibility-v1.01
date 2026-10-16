# INSPIRE arthroplasty anesthesia analysis (v1.02) - 11_coded_complications.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import pandas as pd, numpy as np, os
a=pd.read_pickle(os.path.join(WORK,"arth2.pkl"))
dx=pd.read_pickle(os.path.join(WORK,"dx.pkl")); dx["c"]=dx.icd10_cm.astype(str).str.upper()
m=dx.merge(a[["op_id","subject_id","orin_time","admission_time","discharge_time","neur"]],on="subject_id")
cats={"pneumonia":[f"J{i}" for i in range(12,19)]+["J69"],"PE":["I26"],"DVT/VTE":["I80","I82"],"MI":["I21","I22"],"stroke":["I61","I63","I64"],
"delirium":["F05"],"sepsis":["A40","A41"],"resp_failure":["J96"],"proc_compl_T81":["T81"],"prosthesis_compl_T84":["T84"],"AKI_N17":["N17"],"ileus":["K56"],"UTI":["N39"],"arrhythmia_I48":["I48"],"hf_I50":["I50"]}
pre=m[m.chart_time<m.orin_time]; post=m[(m.chart_time>=m.orin_time)&(m.chart_time<=m.discharge_time)]
post30=m[(m.chart_time>=m.orin_time)&(m.chart_time<=m.orin_time+30*1440)]
for k,codes in cats.items():
    had=set(pre[pre.c.isin(codes)].op_id); 
    ev=set(post[post.c.isin(codes)].op_id)-had; ev30=set(post30[post30.c.isin(codes)].op_id)-had
    g=a.set_index("op_id").neur
    def cnt(s): s=list(s); return int((g.loc[s]==0).sum()), int((g.loc[s]==1).sum())
    print(f"{k:22s} index-stay new: GA/NA {cnt(ev)}   30d new: {cnt(ev30)}   (pre-existing ops {len(had)})")
# timing of post-op dx within index
print("postop index dx timing (days after OR) quantiles:",((post.chart_time-post.orin_time)/1440).quantile([.1,.5,.9]).round(1).tolist(), "n",len(post))
print("dx at exactly discharge time frac:",round((post.chart_time==post.discharge_time).mean(),3))
