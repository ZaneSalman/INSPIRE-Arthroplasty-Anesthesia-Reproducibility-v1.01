# INSPIRE arthroplasty anesthesia analysis (v1.02) - 01_rebuild_cohort.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, math
D=DATA.rstrip("/")+"/"
ops=pd.read_csv(D+"operations.csv.gz")
AR={"0SRC","0SRD","0SR9","0SRB","0SRA","0SRE"}
o=ops[ops.department.eq("OS")]; o=o[o.antype.isin(["General","Neuraxial"])]
a=o[o.icd10_pcs.astype(str).str[:4].isin(AR)].copy()
a["los_raw"]=(a.discharge_time-a.admission_time)/1440
a=a[a.los_raw.gt(0)&a.los_raw.lt(365)].copy(); a["los"]=a.los_raw.clip(lower=.5)
a["icu"]=a.icuin_time.notna().astype(int); a["neur"]=a.antype.eq("Neuraxial").astype(int)
a["female"]=a.sex.eq("F").astype(int); a["asa3"]=np.where(a.asa.isna(),np.nan,(a.asa>2).astype(float))
a["knee"]=a.icd10_pcs.astype(str).str[:4].isin({"0SRC","0SRD"}).astype(int)
a=a.dropna(subset=["los","age","asa3","emop"]).copy()
print("N",len(a),"GA",(a.neur==0).sum(),"NA",(a.neur==1).sum(),"ICU GA/NA",a[a.neur==0].icu.sum(),a[a.neur==1].icu.sum())
# time-reference check
print("orin within opdate day:",float(((a.orin_time-a.opdate*1440).between(0,1440)).mean()),
      "admission<=orin:",float((a.admission_time<=a.orin_time).mean()))
a.to_pickle(os.path.join(WORK,"arth.pkl"))
