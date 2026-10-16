# INSPIRE arthroplasty anesthesia analysis (v1.02) - 03_aki_and_mortality_models.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os, math
a=pd.read_pickle(os.path.join(WORK,"arth.pkl")).reset_index(drop=True)
cr=pd.read_pickle(os.path.join(WORK,"cr.pkl")); cr["value"]=pd.to_numeric(cr.value,errors="coerce"); cr=cr.dropna(subset=["value"]); cr=cr[cr.value.between(0.1,20)]
Z=1.959963984540054
def X_(d): return np.column_stack([np.ones(len(d)),d.neur,d.age,d.female,d.asa3,d.emop,d.knee]).astype(float)
names=["const","neuraxial","age","female","asa3","emergency","knee"]
def pnorm(z): return math.erfc(abs(z)/math.sqrt(2))
def logit(X,y):
    b=np.zeros(X.shape[1])
    for _ in range(100):
        p=1/(1+np.exp(-X@b)); W=p*(1-p); H=X.T@(X*W[:,None]); nb=b+np.linalg.solve(H,X.T@(y-p))
        if np.max(abs(nb-b))<1e-10: b=nb;break
        b=nb
    p=1/(1+np.exp(-X@b)); H=X.T@(X*(p*(1-p))[:,None]); Hi=np.linalg.inv(H)
    meat=X.T@(X*((y-p)**2)[:,None]); V=Hi@meat@Hi; return b,np.sqrt(np.diag(V)),p
def gamma_log(X,y):
    b=np.zeros(X.shape[1]); b[0]=np.log(y.mean())
    for _ in range(200):
        mu=np.exp(X@b); z=X@b+(y-mu)/mu; nb=np.linalg.lstsq(X,z,rcond=None)[0]
        if np.max(abs(nb-b))<1e-10: b=nb;break
        b=nb
    mu=np.exp(X@b); Hi=np.linalg.inv(X.T@X); meat=X.T@(X*(((y-mu)/mu)**2)[:,None]); V=Hi@meat@Hi
    return b,np.sqrt(np.diag(V))
def auc(y,p):
    r=pd.Series(p).rank().values; n1=y.sum(); n0=len(y)-n1; return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def show(label,b,se,rows=("neuraxial",)):
    for nm in rows:
        i=names.index(nm); e=math.exp(b[i]); print(f"  {label} {nm}: {e:.3f} ({math.exp(b[i]-Z*se[i]):.3f}-{math.exp(b[i]+Z*se[i]):.3f}) p={pnorm(b[i]/se[i]):.4f}")
X=X_(a)
b,se,p=logit(X,a.icu.values.astype(float)); print("VALIDATION ICU:"); show("ICU",b,se); print("  c=",round(auc(a.icu.values,p),3))
b,se=gamma_log(X,a.los.values); print("VALIDATION LOS:"); show("LOS",b,se)

# AKI (KDIGO creatinine)
m=cr.merge(a[["op_id","subject_id","orin_time","opend_time"]],on="subject_id")
pre=m[(m.chart_time<m.orin_time)&(m.chart_time>=m.orin_time-30*1440)].sort_values("chart_time").groupby("op_id").value.last().rename("base")
post7=m[(m.chart_time>m.opend_time)&(m.chart_time<=m.opend_time+7*1440)].groupby("op_id").value.max().rename("max7")
post48=m[(m.chart_time>m.opend_time)&(m.chart_time<=m.opend_time+2*1440)].groupby("op_id").value.max().rename("max48")
a=a.merge(pd.concat([pre,post7,post48],axis=1),left_on="op_id",right_index=True,how="left")
k=a.dropna(subset=["base","max7"]).copy()
k["max48"]=k.max48.fillna(k.base)
r=k.max7/k.base
k["aki"]=((r>=1.5)|((k.max48-k.base)>=0.3)).astype(int)
k["stage"]=np.select([(r>=3)|((k.max7>=4)&(k.aki==1)),r>=2,k.aki==1],[3,2,1],0)
print("\nAKI cohort n=",len(k),"excluded (no baseline or no postop cr):",len(a)-len(k), " base>=4.0:",int((k.base>=4).sum()))
for g,lab in [(0,"GA"),(1,"NA")]:
    s=k[k.neur==g]; print(f"  {lab}: AKI {s.aki.sum()}/{len(s)} ({100*s.aki.mean():.1f}%)  stage1/2/3: {(s.stage==1).sum()}/{(s.stage==2).sum()}/{(s.stage==3).sum()}  baseline cr median {s.base.median():.2f}")
b,se,p=logit(X_(k),k.aki.values.astype(float)); show("AKI",b,se,names[1:]); print("  c=",round(auc(k.aki.values,p),3))
# excluding baseline>=4 (possible ESKD)
k2=k[k.base<4]; b,se,p=logit(X_(k2),k2.aki.values.astype(float)); print(" sens excl baseline>=4, n=",len(k2)); show("AKI",b,se)

# mortality
def fisher(a_,b_,c_,d_):
    n=a_+b_+c_+d_; r1=a_+b_; c1=a_+c_
    pr=lambda x: math.comb(r1,x)*math.comb(n-r1,c1-x)/math.comb(n,c1)
    p0=pr(a_); return sum(pr(x) for x in range(max(0,c1-(n-r1)),min(r1,c1)+1) if pr(x)<=p0*(1+1e-7))
a["d_inh"]=a.inhosp_death_time.notna().astype(int)
a["d30"]=((a.allcause_death_time-a.orin_time)/1440).between(0,30).astype(int)
a["d90"]=((a.allcause_death_time-a.orin_time)/1440).between(0,90).astype(int)
for v in ["d_inh","d30","d90"]:
    g=a[a.neur==0][v]; n_=a[a.neur==1][v]
    print(f"\n{v}: GA {g.sum()}/{len(g)} ({100*g.mean():.2f}%)  NA {n_.sum()}/{len(n_)} ({100*n_.mean():.2f}%)  Fisher p={fisher(int(g.sum()),len(g)-int(g.sum()),int(n_.sum()),len(n_)-int(n_.sum())):.3f}")
# hip/knee for inhosp deaths
print("inhosp deaths by joint (knee=1):",a.groupby("knee").d_inh.sum().to_dict(), "emergency among deaths:",int(a[a.d_inh==1].emop.sum()))
# composite: inhosp death or AKI or ICU
k["d_inh"]=k.inhosp_death_time.notna().astype(int)
k["comp"]=((k.aki==1)|(k.d_inh==1)|(k.icu==1)).astype(int)
for g,lab in [(0,"GA"),(1,"NA")]:
    s=k[k.neur==g]; print(f"composite {lab}: {s.comp.sum()}/{len(s)} ({100*s.comp.mean():.1f}%)")
b,se,p=logit(X_(k),k.comp.values.astype(float)); show("Composite",b,se)
