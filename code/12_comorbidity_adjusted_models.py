# INSPIRE arthroplasty anesthesia analysis (v1.02) - 12_comorbidity_adjusted_models.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os, math
src=open(os.path.join(HERE,"04_aki_composite_sensitivity.py")).read().split("# propensity / IPTW on full cohort")[0]
exec(src)
a2=pd.read_pickle(os.path.join(WORK,"arth2.pkl"))[["op_id","icu_idx","d_inh_idx"]]
a=a.merge(a2,on="op_id").merge(pd.read_pickle(os.path.join(WORK,"cci.pkl"))[["op_id","cci"]],on="op_id")
a["cci1"]=(a.cci==1).astype(int); a["cci2"]=(a.cci>=2).astype(int)
a["comp_c"]=np.where(a.aki.notna(),((a.aki==1)|(a.icu_idx==1)|(a.d_inh_idx==1)).astype(float),np.nan)
def gamma_fit(X,y,w=None):
    w=np.ones(len(y)) if w is None else w
    b=np.zeros(X.shape[1]); b[0]=np.log(y.mean())
    for _ in range(300):
        mu=np.exp(X@b); z=X@b+(y-mu)/mu; W=np.sqrt(w); nb=np.linalg.lstsq(X*W[:,None],z*W,rcond=None)[0]
        if np.max(abs(nb-b))<1e-11: b=nb;break
        b=nb
    mu=np.exp(X@b); Hi=np.linalg.inv(X.T@(X*w[:,None])); meat=X.T@(X*(w*((y-mu)/mu)**2)[:,None]); V=Hi@meat@Hi
    return b,np.sqrt(np.diag(V))
def XC(d): return np.column_stack([np.ones(len(d)),d.neur,d.age,d.female,d.asa3,d.emop,d.knee,d.cci1,d.cci2]).astype(float)
def Xps(d,cci=True): return np.column_stack([np.ones(len(d)),d.age,d.female,d.asa3,d.emop,d.knee]+([d.cci1,d.cci2] if cci else [])).astype(float)
outs={"ICU":"icu_idx","AKI":"aki","COMP":"comp_c"}
print("=== + CCI adjustment")
for nm,v in outs.items():
    d=a.dropna(subset=[v]); b,se,_=wlogit(XC(d),d[v].values.astype(float)); print(nm,fmt(b,se), "| CCI>=2 OR:",fmt(b,se,8))
b,se=gamma_fit(XC(a),a.los.values); print("LOS",fmt(b,se))
# PS with CCI: IPTW and overlap weights
_,_,ps=wlogit(Xps(a),a.neur.values.astype(float)); prev=a.neur.mean()
iw=np.where(a.neur==1,prev/ps,(1-prev)/(1-ps)); a["iptw2"]=np.clip(iw,None,np.quantile(iw,.99))
a["ow"]=np.where(a.neur==1,1-ps,ps)
def smd_w(var,w):
    r=[]
    for g in [1,0]:
        z=a[a.neur==g]; ww=z[w].values; m=np.average(z[var],weights=ww); vv=np.average((z[var]-m)**2,weights=ww); r.append((m,vv))
    return abs(r[0][0]-r[1][0])/math.sqrt((r[0][1]+r[1][1])/2)
for w in ["iptw2","ow"]:
    print(w,"maxSMD",round(max(smd_w(v,w) for v in ["age","female","asa3","emop","knee","cci1","cci2"]),3),"ESS",round(a[w].sum()**2/(a[w]**2).sum()))
print("unweighted SMD cci1/cci2:",round(smd_w("cci1","ones") if False else 0,3))
for w in ["iptw2","ow"]:
    print("===",w)
    for nm,v in outs.items():
        d=a.dropna(subset=[v]); b,se,_=wlogit(XC(d),d[v].values.astype(float),d[w].values); print(" ",nm,fmt(b,se))
    b,se=gamma_fit(XC(a),a.los.values,a[w].values); print("  LOS",fmt(b,se))
# low-risk restriction: elective, ASA 1-2
r=a[(a.emop==0)&(a.asa3==0)]; print("=== elective ASA1-2 n",len(r),"GA",int((r.neur==0).sum()),"NA",int((r.neur==1).sum()))
def XR(d): return np.column_stack([np.ones(len(d)),d.neur,d.age,d.female,d.knee,d.cci1,d.cci2]).astype(float)
for nm,v in outs.items():
    d=r.dropna(subset=[v]); ev=(int(d[d.neur==0][v].sum()),int(d[d.neur==1][v].sum()))
    b,se,_=wlogit(XR(d),d[v].values.astype(float)); print(" ",nm,"events GA/NA",ev,fmt(b,se))
b,se=gamma_fit(XR(r),r.los.values); print("  LOS",fmt(b,se))
# g-computation adjusted risks & RD with bootstrap (model with CCI)
rng=np.random.default_rng(2026)
def gcomp(d,v):
    X=XC(d); y=d[v].values.astype(float); b,_,_=wlogit(X,y)
    X1=X.copy(); X1[:,1]=1; X0=X.copy(); X0[:,1]=0
    p1=(1/(1+np.exp(-X1@b))).mean(); p0=(1/(1+np.exp(-X0@b))).mean(); return p0,p1,p1-p0
print("=== g-computation (CCI-adjusted model), risks under general vs neuraxial")
for nm,v in outs.items():
    d=a.dropna(subset=[v]).reset_index(drop=True); p0,p1,rd=gcomp(d,v)
    subj=d.subject_id.unique(); idx=d.groupby("subject_id").indices; bs=[]
    for _ in range(1000):
        s=rng.choice(subj,len(subj),replace=True); ii=np.concatenate([idx[k] for k in s])
        try: bs.append(gcomp(d.iloc[ii],v)[2])
        except Exception: pass
    lo,hi=np.percentile(bs,[2.5,97.5])
    print(f"  {nm}: risk GA {100*p0:.1f}% NA {100*p1:.1f}% RD {100*rd:.1f} pp (95% CI {100*lo:.1f} to {100*hi:.1f}); NNT {abs(1/rd):.0f} ; boots {len(bs)}")
