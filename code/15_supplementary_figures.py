# INSPIRE arthroplasty anesthesia analysis (v1.02) - 15_supplementary_figures.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd, os
from matplotlib.backends.backend_pdf import PdfPages
src=open(os.path.join(HERE,"04_aki_composite_sensitivity.py")).read().split("# propensity / IPTW on full cohort")[0]
exec(src)
Xp=np.column_stack([np.ones(len(a)),a.age,a.female,a.asa3,a.emop,a.knee]).astype(float)
_,_,ps=wlogit(Xp,a.neur.values.astype(float)); prev=a.neur.mean()
iw=np.where(a.neur==1,prev/ps,(1-prev)/(1-ps)); cut=np.quantile(iw,.99); tw=np.clip(iw,None,cut)
OUT=WORK+"/"
cap=a.los.quantile(.97)
with PdfPages(OUT+"Supplementary_Figures_S1-S3.pdf") as pdf:
    fig,ax=plt.subplots(figsize=(7.5,5))
    bins=np.linspace(0,cap,40)
    for g,l,c in [(0,"General (n = 269)","#c55a11"),(1,"Neuraxial (n = 2271)","#1f4e79")]:
        ax.hist(a[a.neur==g].los.clip(upper=cap),bins=bins,density=True,alpha=.55,color=c,label=l)
    ax.set_xlabel("Hospital length of stay, days (capped at the 97th percentile for display)"); ax.set_ylabel("Density"); ax.legend(frameon=False)
    ax.set_title("Figure S1. Distribution of hospital length of stay by anesthesia type",fontsize=11); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.5,5)); bins=np.linspace(0,1,41)
    for g,l,c in [(0,"General","#c55a11"),(1,"Neuraxial","#1f4e79")]:
        ax.hist(ps[a.neur.values==g],bins=bins,density=True,alpha=.55,color=c,label=l)
    ax.set_xlabel("Propensity score for neuraxial anesthesia"); ax.set_ylabel("Density"); ax.legend(frameon=False)
    ax.set_title("Figure S2. Propensity-score distributions before weighting",fontsize=11); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(9,4.5),sharey=True)
    ax[0].hist(iw,bins=40,color="#7f7f7f"); ax[0].set_title("Before truncation",fontsize=10); ax[0].set_xlabel("Stabilized weight"); ax[0].set_ylabel("Operations")
    ax[1].hist(tw,bins=40,color="#1f4e79"); ax[1].set_title(f"After truncation at 99th percentile ({cut:.2f})",fontsize=10); ax[1].set_xlabel("Stabilized weight")
    fig.suptitle("Figure S3. Stabilized IPTW weights before and after truncation",fontsize=11); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
print("cap",round(cap,1),"cut",round(cut,2),"ESS",round(tw.sum()**2/(tw**2).sum()),"ps range",round(ps.min(),3),round(ps.max(),3))
