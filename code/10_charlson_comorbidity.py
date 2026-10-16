# INSPIRE arthroplasty anesthesia analysis (v1.02) - 10_charlson_comorbidity.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import pandas as pd, numpy as np, os
a=pd.read_pickle(os.path.join(WORK,"arth2.pkl"))
dx=pd.read_pickle(os.path.join(WORK,"dx.pkl")); dx["c"]=dx.icd10_cm.astype(str).str.upper()
def rng(L,a_,b_): return [f"{L}{i:02d}" for i in range(a_,b_+1)]
C={"mi":(["I21","I22"],1),"chf":(["I50","I43","I42"],1),"pvd":(["I70","I71","I73"],1),"cvd":(rng("I",60,69)+["G45","G46"],1),
"dementia":(rng("F",0,3)+["G30"],1),"pulm":(rng("J",40,47)+rng("J",60,67),1),"rheum":(["M05","M06","M32","M33","M34"],1),"pud":(rng("K",25,28),1),
"mildliver":(["B18","K73","K74"],1),"dm":(rng("E",10,14),1),"hemi":(["G81","G82"],2),"renal":(["N18","N19","Z49"],2),
"cancer":(rng("C",0,76)+rng("C",81,97),2),"sevliver":(["K72","I85"],3),"mets":(rng("C",77,80),6),"hiv":(rng("B",20,24),6)}
m=dx.merge(a[["op_id","subject_id","orin_time"]],on="subject_id"); pre=m[m.chart_time<m.orin_time]
flags=pd.DataFrame(index=a.op_id)
for k,(codes,w) in C.items(): flags[k]=flags.index.isin(set(pre[pre.c.isin(codes)].op_id)).astype(int)
flags.loc[flags.sevliver==1,"mildliver"]=0; flags.loc[flags.mets==1,"cancer"]=0
w={k:v[1] for k,v in C.items()}
flags["cci"]=sum(flags[k]*w[k] for k in C)
flags=flags.reset_index().rename(columns={"index":"op_id"})
flags.to_pickle(os.path.join(WORK,"cci.pkl"))
x=a.merge(flags,on="op_id")
for g,l in [(0,"GA"),(1,"NA")]:
    s=x[x.neur==g]; print(l,"CCI median",s.cci.median(),"IQR",s.cci.quantile(.25),s.cci.quantile(.75),"mean",round(s.cci.mean(),2),"CCI>=2 %",round(100*(s.cci>=2).mean(),1), {k:int(s[k].sum()) for k in ["dm","renal","chf","cvd","pulm","cancer"]})
g0=x[x.neur==0].cci; g1=x[x.neur==1].cci; import math
print("SMD cci",round((g0.mean()-g1.mean())/math.sqrt((g0.var()+g1.var())/2),3))
print("cci>=2 counts GA",int((g0>=2).sum()),"NA",int((g1>=2).sum()))
