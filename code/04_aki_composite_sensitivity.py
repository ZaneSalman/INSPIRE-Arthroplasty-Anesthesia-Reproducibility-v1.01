# INSPIRE arthroplasty anesthesia analysis (v1.02) - 04_aki_composite_sensitivity.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os, math
exec(open(os.path.join(HERE,"03_aki_and_mortality_models.py")).read().split("# AKI (KDIGO creatinine)")[0].replace("print(","(lambda *x,**y:None)("))
a=pd.read_pickle(os.path.join(WORK,"arth.pkl")).reset_index(drop=True)
cr=pd.read_pickle(os.path.join(WORK,"cr.pkl")); cr["value"]=pd.to_numeric(cr.value,errors="coerce"); cr=cr.dropna(subset=["value"]); cr=cr[cr.value.between(0.1,20)]
m=cr.merge(a[["op_id","subject_id","orin_time","opend_time"]],on="subject_id")
pre=m[(m.chart_time<m.orin_time)&(m.chart_time>=m.orin_time-30*1440)].sort_values("chart_time").groupby("op_id").value.last().rename("base")
post7=m[(m.chart_time>m.opend_time)&(m.chart_time<=m.opend_time+7*1440)].groupby("op_id").value.max().rename("max7")
post48=m[(m.chart_time>m.opend_time)&(m.chart_time<=m.opend_time+2*1440)].groupby("op_id").value.max().rename("max48")
a=a.merge(pd.concat([pre,post7,post48],axis=1),left_on="op_id",right_index=True,how="left")
a["has_aki_data"]=a.base.notna()&a.max7.notna()
a["max48"]=a.max48.fillna(a.base); r=a.max7/a.base
a["aki"]=np.where(a.has_aki_data,((r>=1.5)|((a.max48-a.base)>=0.3)).astype(float),np.nan)
a["d_inh"]=a.inhosp_death_time.notna().astype(int)
a["comp"]=np.where(a.has_aki_data,((a.aki==1)|(a.d_inh==1)|(a.icu==1)).astype(float),np.nan)
a["bmi"]=a.weight/(a.height/100)**2; a.loc[~a.bmi.between(10,70),"bmi"]=np.nan

def wlogit(X,y,w=None,groups=None,mode="w"):
    w=np.ones(len(y)) if w is None else w
    b=np.zeros(X.shape[1])
    for _ in range(200):
        p=1/(1+np.exp(-X@b)); H=X.T@(X*(w*p*(1-p))[:,None]); nb=b+np.linalg.solve(H,X.T@(w*(y-p)))
        if np.max(abs(nb-b))<1e-11: b=nb;break
        b=nb
    p=1/(1+np.exp(-X@b)); H=X.T@(X*(w*p*(1-p))[:,None]); Hi=np.linalg.inv(H)
    if groups is None:
        f=(w if mode=="w" else w**2)*(y-p)**2; meat=X.T@(X*f[:,None])
    else:
        S=pd.DataFrame(X*(w*(y-p))[:,None]).groupby(groups.values).sum().values
        G=S.shape[0]; n,K=X.shape; meat=S.T@S*(G/(G-1))*((n-1)/(n-K))
    V=Hi@meat@Hi; return b,np.sqrt(np.diag(V)),p
def fmt(b,se,i=1): return f"{math.exp(b[i]):.3f} ({math.exp(b[i]-Z*se[i]):.3f}-{math.exp(b[i]+Z*se[i]):.3f}) p={pnorm(b[i]/se[i]):.4f}"
def XB(d,bmi=False):
    cols=[np.ones(len(d)),d.neur,d.age,d.female,d.asa3,d.emop,d.knee]+([d.bmi] if bmi else [])
    return np.column_stack(cols).astype(float)
def XB_nok(d,bmi=True): return np.column_stack([np.ones(len(d)),d.neur,d.age,d.female,d.asa3,d.emop]+([d.bmi] if bmi else [])).astype(float)
# propensity / IPTW on full cohort
Xp=np.column_stack([np.ones(len(a)),a.age,a.female,a.asa3,a.emop,a.knee]).astype(float)
bp,_,ps=wlogit(Xp,a.neur.values.astype(float)); prev=a.neur.mean()
iw=np.where(a.neur==1,prev/ps,(1-prev)/(1-ps)); a["iptw"]=np.clip(iw,None,np.quantile(iw,.99))
print("ESS",round(a.iptw.sum()**2/(a.iptw**2).sum()))
for mode in ["w","w2"]:
    b,se,_=wlogit(XB(a),a.icu.values.astype(float),a.iptw.values,mode=mode); print("VALID ICU IPTW",mode,fmt(b,se),"(target 0.440 0.272-0.710)")
h=a[(a.knee==0)].dropna(subset=["bmi"]); b,se,_=wlogit(XB_nok(h),h.icu.values.astype(float),groups=h.subject_id); print("VALID hip ICU bmi+clust",len(h),fmt(b,se),"(target 0.406 0.241-0.684, n=1288)")
for out in ["aki","comp"]:
    k=a.dropna(subset=[out]); y=k[out].values
    b,se,_=wlogit(XB(k),y); print(f"\n{out.upper()} adjusted n={len(k)}:",fmt(b,se))
    b,se,_=wlogit(XB(k),y,k.iptw.values,mode="w"); print(f"{out.upper()} IPTW:",fmt(b,se))
    kb=k.dropna(subset=["bmi"]); b,se,_=wlogit(XB(kb,True),kb[out].values,groups=kb.subject_id); print(f"{out.upper()} BMI+cluster n={len(kb)}:",fmt(b,se))
    k1=k.sort_values("orin_time").groupby("subject_id").head(1); b,se,_=wlogit(XB(k1),k1[out].values); print(f"{out.upper()} one op/patient n={len(k1)}:",fmt(b,se))
    for kn,lab in [(0,"hip"),(1,"knee")]:
        s=k[k.knee==kn].dropna(subset=["bmi"]); ev=int(s[out].sum())
        try: b,se,_=wlogit(XB_nok(s),s[out].values,groups=s.subject_id); print(f"  {lab} n={len(s)} events={ev}:",fmt(b,se))
        except Exception as e: print(lab,"fail",e)
    # interaction
    Xi=np.column_stack([XB(k),k.neur*k.knee]).astype(float); b,se,_=wlogit(Xi,y); print("  interaction p=",round(pnorm(b[-1]/se[-1]),3))
print("subjects:",a.subject_id.nunique(),"ops:",len(a))
# AKI subgroup counts by group for table
k=a.dropna(subset=["aki"])
print("AKI hip/knee by group:",k.groupby(["knee","neur"]).aki.agg(["sum","count"]).to_dict())
