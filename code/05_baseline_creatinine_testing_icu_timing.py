# INSPIRE arthroplasty anesthesia analysis (v1.02) - 05_baseline_creatinine_testing_icu_timing.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os, math
src=open(os.path.join(HERE,"04_aki_composite_sensitivity.py")).read().split("# propensity / IPTW on full cohort")[0]
exec(src)
# baseline creatinine descriptive
for g,l in [(0,"GA"),(1,"NA")]:
    s=a[a.neur==g].base.dropna(); print(f"{l} baseline cr n={len(s)} median {s.median():.2f} [{s.quantile(.25):.2f}-{s.quantile(.75):.2f}] mean {s.mean():.2f} sd {s.std():.2f} ; >=1.5: {(s>=1.5).sum()}")
g0=a[a.neur==0].base.dropna(); g1=a[a.neur==1].base.dropna()
print("SMD base", (g0.mean()-g1.mean())/math.sqrt((g0.var()+g1.var())/2))
print("baseline missing:",int(a.base.isna().sum()))
# postop creatinine testing intensity within 7 days
m2=cr.merge(a[["op_id","subject_id","opend_time","neur"]],on="subject_id")
cnt=m2[(m2.chart_time>m2.opend_time)&(m2.chart_time<=m2.opend_time+7*1440)].groupby("op_id").size()
a["ncr7"]=a.op_id.map(cnt).fillna(0)
for g,l in [(0,"GA"),(1,"NA")]:
    s=a[a.neur==g].ncr7; print(f"{l} postop cr tests/7d median {s.median():.0f} [{s.quantile(.25):.0f}-{s.quantile(.75):.0f}] mean {s.mean():.2f}; >=1 test {100*(s>=1).mean():.1f}%")
# AKI adjusted for baseline creatinine
k=a.dropna(subset=["aki"]).copy()
def X2(d,extra): return np.column_stack([np.ones(len(d)),d.neur,d.age,d.female,d.asa3,d.emop,d.knee]+extra).astype(float)
b,se,_=wlogit(X2(k,[k.base]),k.aki.values); print("AKI +baseline cr:",fmt(b,se))
b,se,_=wlogit(X2(k,[k.base,k.ncr7]),k.aki.values); print("AKI +baseline cr +n tests:",fmt(b,se))
# ICU timing
ic=a[a.icu==1].copy(); ic["h"]=(ic.icuin_time-ic.orout_time)/60
print("ICU hours after OR exit quantiles:",ic.h.quantile([0,.1,.25,.5,.75,.9,1]).round(1).tolist())
for thr in [6,24]:
    ic["late"]=ic.h>thr
    print(f"thr {thr}h: early GA/NA {((~ic.late)&(ic.neur==0)).sum()}/{((~ic.late)&(ic.neur==1)).sum()}  late GA/NA {(ic.late&(ic.neur==0)).sum()}/{(ic.late&(ic.neur==1)).sum()}")
a["h_icu"]=(a.icuin_time-a.orout_time)/60
for thr in [6,24]:
    for nm,cond in [("early",a.icu.eq(1)&(a.h_icu<=thr)),("late",a.icu.eq(1)&(a.h_icu>thr))]:
        y=cond.astype(float).values
        try: b,se,_=wlogit(X_(a),y); print(f"  {nm} ICU (<= / > {thr}h) events={int(y.sum())}:",fmt(b,se))
        except Exception as e: print(nm,e)
# E-values
def ev(rr): rr=1/rr if rr<1 else rr; return rr+math.sqrt(rr*(rr-1))
def evalue(or_,hi,common):
    f=(lambda x: math.sqrt(x)) if common else (lambda x:x)
    pt=ev(f(or_)); ci=1.0 if hi>=1 else ev(f(hi)); return round(pt,2),round(ci,2)
print("E ICU (rare approx):",evalue(0.494,0.775,False))
print("E AKI:",evalue(0.516,0.926,False))
print("E composite (common, sqrt):",evalue(0.482,0.709,True))
print("E ICU sqrt (conservative):",evalue(0.494,0.775,True))
# sex-disaggregated
a["comp2"]=a.comp
for sx,l in [(1,"Female"),(0,"Male")]:
    s=a[a.female==sx]
    for g,gl in [(0,"GA"),(1,"NA")]:
        t=s[s.neur==g]; kk=t.dropna(subset=["aki"])
        print(f"{l} {gl} n={len(t)} ICU {int(t.icu.sum())} ({100*t.icu.mean():.1f}) AKI {int(kk.aki.sum())}/{len(kk)} ({100*kk.aki.mean():.1f}) LOS med {t.los.median():.1f} inhosp death {int(t.d_inh.sum())}")
# sex interaction for ICU & AKI
Xi=np.column_stack([X_(a),a.neur*a.female]); b,se,_=wlogit(Xi,a.icu.values.astype(float)); print("ICU sex interaction p",round(pnorm(b[-1]/se[-1]),3))
Xi=np.column_stack([X_(k),k.neur*k.female]); b,se,_=wlogit(Xi,k.aki.values); print("AKI sex interaction p",round(pnorm(b[-1]/se[-1]),3))
