# INSPIRE arthroplasty anesthesia analysis (v1.02) - 06_icu_timing_audit.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os
a=pd.read_pickle(os.path.join(WORK,"arth.pkl"))
ic=a[a.icu==1].copy()
ic["in_index"]=(ic.icuin_time>=ic.admission_time)&(ic.icuin_time<=ic.discharge_time)
ic["after_or"]=ic.icuin_time>=ic.orout_time-60
ic["h"]=(ic.icuin_time-ic.orout_time)/60
print("ICU ops:",len(ic)," within index admission:",int(ic.in_index.sum())," before index admission:",int((ic.icuin_time<ic.admission_time).sum())," after discharge:",int((ic.icuin_time>ic.discharge_time).sum()))
print("within index & <=24h of OR exit:",int((ic.in_index&(ic.h<=24)).sum()), " within index >24h:",int((ic.in_index&(ic.h>24)).sum()))
for g,l in [(0,"GA"),(1,"NA")]:
    s=ic[ic.neur==g]; print(l,"total",len(s),"in-index",int(s.in_index.sum()),"post-discharge",int((s.icuin_time>s.discharge_time).sum()))
# other ops for same subject: does icuin_time equal another op's window?
print("days from discharge to ICU (post-discharge cases) quantiles:",((ic[ic.icuin_time>ic.discharge_time].icuin_time-ic[ic.icuin_time>ic.discharge_time].discharge_time)/1440).quantile([.1,.5,.9]).round(0).tolist())
ops=pd.read_csv(os.path.join(DATA,"operations.csv.gz"),usecols=["op_id","subject_id","icuin_time","orout_time","admission_time","discharge_time","department"])
ops["h"]=(ops.icuin_time-ops.orout_time)/60
ops["in_index"]=(ops.icuin_time>=ops.admission_time)&(ops.icuin_time<=ops.discharge_time)
x=ops[ops.icuin_time.notna()]
print("ALL INSPIRE ops with icuin:",len(x),"in index admission frac:",round(x.in_index.mean(),3),"median h after OR:",round(x.h.median(),1))
