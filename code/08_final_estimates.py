# INSPIRE arthroplasty anesthesia analysis (v1.02) - 08_final_estimates.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import numpy as np, pandas as pd, os, math
src=open(os.path.join(HERE,"04_aki_composite_sensitivity.py")).read().split("# propensity / IPTW on full cohort")[0]
exec(src)
a2=pd.read_pickle(os.path.join(WORK,"arth2.pkl"))[["op_id","icu_idx","d_inh_idx","iptw","icu_early"]]
a=a.merge(a2,on="op_id")
def pv(b,se,i=1): return pnorm(b[i]/se[i])
names=["const","neuraxial","age","female","asa3","emergency","knee"]
y=a.icu_idx.values.astype(float); X=XB(a)
b,se,_=wlogit(X,y)
for i in range(1,7): print("ICU",names[i],fmt(b,se,i))
# deaths index
for g,l in [(0,"GA"),(1,"NA")]:
    s=a[a.neur==g]; print(l,"inhosp death (index)",int(s.d_inh_idx.sum()),"/",len(s),round(100*s.d_inh_idx.mean(),2))
def fisher(a_,b_,c_,d_):
    n=a_+b_+c_+d_; r1=a_+b_; c1=a_+c_
    pr=lambda x: math.comb(r1,x)*math.comb(n-r1,c1-x)/math.comb(n,c1)
    p0=pr(a_); return sum(pr(x) for x in range(max(0,c1-(n-r1)),min(r1,c1)+1) if pr(x)<=p0*(1+1e-7))
g=a[a.neur==0].d_inh_idx; n_=a[a.neur==1].d_inh_idx
print("fisher inhosp idx", round(fisher(int(g.sum()),len(g)-int(g.sum()),int(n_.sum()),len(n_)-int(n_.sum())),3))
print("inhosp idx deaths hip/knee:",a.groupby("knee").d_inh_idx.sum().to_dict(),"emergency:",int(a[a.d_inh_idx==1].emop.sum()))
# LOS gamma extras
def gamma_fit(X,y,groups=None):
    b=np.zeros(X.shape[1]); b[0]=np.log(y.mean())
    for _ in range(300):
        mu=np.exp(X@b); z=X@b+(y-mu)/mu; nb=np.linalg.lstsq(X,z,rcond=None)[0]
        if np.max(abs(nb-b))<1e-11: b=nb;break
        b=nb
    mu=np.exp(X@b); Hi=np.linalg.inv(X.T@X); r=(y-mu)/mu
    if groups is None: meat=X.T@(X*(r**2)[:,None])
    else:
        S=pd.DataFrame(X*r[:,None]).groupby(groups.values).sum().values; G=S.shape[0]; n,K=X.shape; meat=S.T@S*(G/(G-1))*((n-1)/(n-K))
    V=Hi@meat@Hi; return b,np.sqrt(np.diag(V))
b,se=gamma_fit(X,a.los.values); print("LOS adj check",fmt(b,se))
ab=a.dropna(subset=["bmi"]); b,se=gamma_fit(XB(ab,True),ab.los.values,ab.subject_id); print("LOS BMI+clust",len(ab),fmt(b,se))
f1=a.sort_values("orin_time").groupby("subject_id").head(1); b,se=gamma_fit(XB(f1),f1.los.values); print("LOS first-op",len(f1),fmt(b,se))
# weighted gamma: weights via IRLS with weights
def wgamma(X,y,w):
    b=np.zeros(X.shape[1]); b[0]=np.log(y.mean())
    for _ in range(300):
        mu=np.exp(X@b); z=X@b+(y-mu)/mu; W=np.sqrt(w); nb=np.linalg.lstsq(X*W[:,None],z*W,rcond=None)[0]
        if np.max(abs(nb-b))<1e-11: b=nb;break
        b=nb
    mu=np.exp(X@b); Hi=np.linalg.inv(X.T@(X*w[:,None])); meat=X.T@(X*(w*((y-mu)/mu)**2)[:,None]); V=Hi@meat@Hi
    return b,np.sqrt(np.diag(V))
b,se=wgamma(X,a.los.values,a.iptw.values); print("LOS IPTW check (0.900 0.842-0.961)",fmt(b,se))
# Balance SMDs
def wsm(var,binary,w=None):
    vals=[]
    for gg in [1,0]:
        z=a[a.neur==gg]; ww=np.ones(len(z)) if w is None else z[w].values
        m=np.average(z[var],weights=ww); v=np.average((z[var]-m)**2,weights=ww) if not binary else m*(1-m)
        if not binary and w is None: v=z[var].var()
        vals.append((m,v))
    return abs((vals[0][0]-vals[1][0])/math.sqrt((vals[0][1]+vals[1][1])/2))
for lab,var,bi in [("Age","age",False),("Female","female",True),("Emergency","emop",True),("Knee","knee",True),("ASA 3+","asa3",True)]:
    print("BAL",lab,round(wsm(var,bi),3),round(wsm(var,bi,"iptw"),3))
# E-values
def ev(rr): rr=1/rr if rr<1 else rr; return rr+math.sqrt(rr*(rr-1))
for nm,o,hi in [("ICU",0.472,0.930),("AKI",0.516,0.926),("COMP",0.471,0.759)]:
    print("E",nm,round(ev(o),2),round(ev(hi),2))
# sex disaggregated corrected
a["comp_c"]=np.where(a.aki.notna(),((a.aki==1)|(a.icu_idx==1)|(a.d_inh_idx==1)).astype(float),np.nan)
for sx,l in [(1,"Female"),(0,"Male")]:
    s=a[a.female==sx]
    for gg,gl in [(0,"General"),(1,"Neuraxial")]:
        t=s[s.neur==gg]; kk=t.dropna(subset=["aki"])
        print(f"SEX {l} {gl} n={len(t)} ICU {int(t.icu_idx.sum())} ({100*t.icu_idx.mean():.1f}) AKI {int(kk.aki.sum())}/{len(kk)} ({100*kk.aki.mean():.1f}) COMP {int(kk.comp_c.sum())}/{len(kk)} ({100*kk.comp_c.mean():.1f}) LOSmed {t.los.median():.1f} [{t.los.quantile(.25):.1f}-{t.los.quantile(.75):.1f}]")
k=a.dropna(subset=["aki"])
for out,d in [("ICU",a),("AKI",k),("COMP",k)]:
    yy=d[{"ICU":"icu_idx","AKI":"aki","COMP":"comp_c"}[out]].values.astype(float)
    Xi=np.column_stack([XB(d),d.neur*d.female]); b,se,_=wlogit(Xi,yy); print("sex interaction",out,round(pnorm(b[-1]/se[-1]),3))
# AKI baseline creatinine adjusted + testing
print("ICU timing within index: <=24h",int(a.icu_early.sum()),"total",int(a.icu_idx.sum()))
# ICU c-stat for check
Xi=np.column_stack([X,a.neur*a.knee]); b,se=gamma_fit(Xi,a.los.values); print("LOS interaction p (unclustered, no BMI)",round(pnorm(b[-1]/se[-1]),3))
ab=a.dropna(subset=["bmi"]); Xi=np.column_stack([XB(ab,True),ab.neur*ab.knee]); b,se=gamma_fit(Xi,ab.los.values,ab.subject_id); print("LOS interaction p (BMI, clustered)",round(pnorm(b[-1]/se[-1]),3))
yy=ab.icu_old.values if "icu_old" in ab else None
Xi=np.column_stack([X,a.neur*a.knee]); b,se,_=wlogit(Xi,a.icu.values.astype(float)); print("OLD ICU interaction p",round(pnorm(b[-1]/se[-1]),3))
k=a.dropna(subset=["aki"]); b,se,_=wlogit(np.column_stack([XB(k),k.base]),k.comp_c.values); print("COMP + baseline cr",fmt(b,se))
