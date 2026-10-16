# INSPIRE arthroplasty anesthesia analysis (v1.02) - 14_main_figures.py
import os
DATA = os.environ.get("INSPIRE_DATA_DIR", "inspire_data")      # folder holding INSPIRE v1.4.2 *.csv.gz
WORK = os.environ.get("INSPIRE_WORK_DIR", "inspire_work")      # local, patient-level intermediates (never commit)
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "code"
os.makedirs(WORK, exist_ok=True)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, os
from matplotlib.patches import FancyBboxPatch
OUT=os.path.join(WORK,"figures")+"/"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
# Figure 1 flow
fig,ax=plt.subplots(figsize=(7,8.2)); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,12)
def box(x,y,w,h,t,fc="#eef3f8"):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.02,rounding_size=0.15",fc=fc,ec="#2f4f6f",lw=1.2)); ax.text(x,y,t,ha="center",va="center",fontsize=9.5)
def arr(x1,y1,x2,y2): ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="->",lw=1.2,color="#2f4f6f"))
box(3.6,11.2,5.6,0.9,"Orthopedic operations in INSPIRE v1.4.2\n(2011–2020), n = 17,433")
box(8.0,10.0,3.6,0.9,"Excluded: anesthesia other than\ngeneral or neuraxial, n = 211",fc="#f8eeee")
box(3.6,8.8,5.6,0.9,"General or neuraxial anesthesia\nn = 17,222")
box(8.0,7.6,3.6,0.9,"Excluded: not hip or knee\narthroplasty, n = 14,666",fc="#f8eeee")
box(3.6,6.4,5.6,0.9,"Hip or knee arthroplasty\nn = 2556")
box(8.0,5.2,3.6,0.9,"Excluded: invalid/missing LOS\nor covariates, n = 16",fc="#f8eeee")
box(3.6,4.0,5.6,0.9,"Primary analytic cohort\nn = 2540",fc="#e3efe3")
box(1.9,2.0,3.2,0.9,"General anesthesia\nn = 269")
box(5.3,2.0,3.2,0.9,"Neuraxial anesthesia\nn = 2271")
box(3.6,0.5,6.4,0.8,"AKI and composite analyses: n = 2525 (15 without baseline\nor postoperative creatinine excluded)",fc="#fdf6e3")
for y1,y2 in [(10.75,9.25),(8.35,6.85),(5.95,4.45)]: arr(3.6,y1,3.6,y2)
for y in [10.0,7.6,5.2]: arr(3.6,y,6.2,y)
arr(3.6,3.55,1.9,2.45); arr(3.6,3.55,5.3,2.45); arr(3.6,1.55,3.6,0.9)
fig.savefig(OUT+"Figure_1_Cohort_Flow_600dpi.png",dpi=600,bbox_inches="tight"); r=fig.canvas.get_renderer(); bb=[t.get_window_extent(r) for t in ax.texts]; print("flow overlaps:",[(i,j) for i in range(len(bb)) for j in range(i+1,len(bb)) if bb[i].overlaps(bb[j])]); pb=[p.get_window_extent(r) for p in ax.patches]; print("box overlaps:",[(i,j) for i in range(len(pb)) for j in range(i+1,len(pb)) if pb[i].overlaps(pb[j])]); plt.close(fig)
# Figure 2 forest
rows=[("ICU admission","Adjusted",0.472,0.239,0.930),("ICU admission","IPTW",0.370,0.179,0.764),("ICU admission","Comorbidity-adjusted",0.518,0.264,1.013),
("Acute kidney injury","Adjusted",0.516,0.287,0.926),("Acute kidney injury","IPTW",0.517,0.274,0.972),("Acute kidney injury","Comorbidity-adjusted",0.691,0.385,1.242),
("Composite adverse event","Adjusted",0.471,0.292,0.759),("Composite adverse event","IPTW",0.429,0.257,0.715),("Composite adverse event","Comorbidity-adjusted",0.569,0.353,0.917),
("Length of stay (ratio of means)","Adjusted",0.892,0.827,0.963),("Length of stay (ratio of means)","IPTW",0.900,0.842,0.961),("Length of stay (ratio of means)","Comorbidity-adjusted",0.912,0.847,0.982)]
fig,ax=plt.subplots(figsize=(9,6.2)); n=len(rows); ys=np.arange(n)[::-1]
labels=[]
for y,(o,m,e,l,h) in zip(ys,rows):
    c={"Adjusted":"#1f4e79","IPTW":"#c55a11","Comorbidity-adjusted":"#548235"}[m]; mk={"Adjusted":"s","IPTW":"o","Comorbidity-adjusted":"D"}[m]
    ax.plot([l,h],[y,y],color=c,lw=1.6); ax.plot(e,y,mk,color=c,ms=6)
    ax.text(2.6,y,f"{e:.3f} ({l:.3f}–{h:.3f})",va="center",fontsize=8.5)
    labels.append(f"{o}, {m if m=='IPTW' else m.lower()}")
ax.axvline(1,color="grey",ls="--",lw=1); ax.set_xscale("log"); ax.set_xlim(0.15,6.0)
ax.set_xticks([0.2,0.5,1,2]); ax.set_xticklabels(["0.2","0.5","1","2"]); ax.minorticks_off()
ax.set_yticks(ys); ax.set_yticklabels(labels,fontsize=8.5); ax.tick_params(axis="y",length=0)
ax.set_xlabel("Neuraxial vs. general anesthesia: odds ratio (length of stay: ratio of means), log scale")
for sp in ["top","right","left"]: ax.spines[sp].set_visible(False)
ax.text(2.6,n-0.4,"Estimate (95% CI)",fontsize=8.5,fontweight="bold")
ax.text(0.95,-1.2,"Favors neuraxial",ha="right",fontsize=8.5); ax.text(1.05,-1.2,"Favors general",ha="left",fontsize=8.5)
ax.set_ylim(-1.6,n-0.2)
fig.savefig(OUT+"Figure_2_Forest_Plot_600dpi.png",dpi=600,bbox_inches="tight"); 
r=fig.canvas.get_renderer(); bbs=[t.get_window_extent(r) for t in ax.texts]
ov=[(i,j) for i in range(len(bbs)) for j in range(i+1,len(bbs)) if bbs[i].overlaps(bbs[j])]
lines=[l.get_window_extent(r) for l in ax.lines if len(l.get_xdata())==2 and l.get_xdata()[0]!=l.get_xdata()[1]]
ovl=[i for i,b in enumerate(bbs) if any(b.overlaps(L) for L in lines)]
print("forest text overlaps:",ov,"text-line overlaps:",ovl)
plt.close(fig)
# Figure 3 balance
cov=["Age","Female","Emergency","Knee","ASA ≥3"]; bef=[0.152,0.201,0.242,0.339,0.402]; aft=[0.049,0.040,0.017,0.033,0.013]
fig,ax=plt.subplots(figsize=(6.4,3.4)); y=np.arange(len(cov))
ax.plot(bef,y,"o",color="#c55a11",label="Before weighting"); ax.plot(aft,y,"s",color="#1f4e79",label="After IPTW")
ax.axvline(0.10,ls="--",color="grey"); ax.set_yticks(y); ax.set_yticklabels(cov); ax.set_xlabel("Absolute standardized mean difference")
ax.set_xlim(0,0.45); ax.legend(frameon=False,loc="lower right"); [ax.spines[s].set_visible(False) for s in ["top","right"]]
fig.savefig(OUT+"Figure_3_Balance_Plot_600dpi.png",dpi=600,bbox_inches="tight"); plt.close(fig)
# Graphical abstract (>=1100x560 px)
fig=plt.figure(figsize=(11,5.6)); ax=fig.add_axes([0,0,1,1]); ax.axis("off"); ax.set_xlim(0,11); ax.set_ylim(0,5.6)
ax.text(5.5,5.15,"Neuraxial vs. general anesthesia and postoperative adverse events after hip and knee arthroplasty",ha="center",fontsize=13,fontweight="bold",color="#1f3b57")
ax.add_patch(FancyBboxPatch((0.3,0.6),3.0,3.9,boxstyle="round,pad=0.05",fc="#eef3f8",ec="#2f4f6f"))
ax.text(1.8,4.1,"INSPIRE cohort",ha="center",fontsize=12,fontweight="bold"); ax.text(1.8,3.5,"2540 hip or knee\narthroplasty operations\nSeoul, 2011–2020",ha="center",va="top",fontsize=11)
ax.text(1.8,1.75,"General: 269\nNeuraxial: 2271",ha="center",va="top",fontsize=11,color="#333")
sub=fig.add_axes([0.37,0.17,0.40,0.62]); lab=["ICU admission","Acute kidney\ninjury","Composite\nadverse event"]
g=[7.1,6.7,12.7]; nn=[1.7,2.4,3.9]; x=np.arange(3); w=0.36
sub.bar(x-w/2,g,w,color="#c55a11",label="General"); sub.bar(x+w/2,nn,w,color="#1f4e79",label="Neuraxial")
for i in range(3): sub.text(x[i]-w/2,g[i]+0.3,f"{g[i]}%",ha="center",fontsize=10); sub.text(x[i]+w/2,nn[i]+0.3,f"{nn[i]}%",ha="center",fontsize=10)
sub.set_xticks(x); sub.set_xticklabels(lab,fontsize=10); sub.set_ylabel("Operations with event (%)"); sub.set_ylim(0,15); sub.legend(frameon=False)
[sub.spines[s].set_visible(False) for s in ["top","right"]]
ax.add_patch(FancyBboxPatch((8.7,0.6),2.0,3.9,boxstyle="round,pad=0.05",fc="#e3efe3",ec="#2f4f6f"))
ax.text(9.7,4.1,"Adjusted OR",ha="center",fontsize=12,fontweight="bold")
ax.text(9.7,3.6,"Primary (+comorbidity)\n\nICU 0.47 (0.52)\nAKI 0.52 (0.69)\nComposite 0.47 (0.57)\nLOS ratio 0.89 (0.91)",ha="center",va="top",fontsize=10)
ax.text(5.5,0.2,"Observational associations; anesthesia assignment was nonrandom.",ha="center",fontsize=10,style="italic",color="#444")
fig.savefig(OUT+"Graphical_Abstract.png",dpi=200); r=fig.canvas.get_renderer(); bb=[t.get_window_extent(r) for t in ax.texts]+[t.get_window_extent(r) for t in sub.texts]; print("GA text overlaps:",[(i,j) for i in range(len(bb)) for j in range(i+1,len(bb)) if bb[i].overlaps(bb[j])]); sb=sub.get_tightbbox(r); fb=fig.bbox; print("GA texts outside fig:",[t.get_text()[:20] for t in ax.texts if not (t.get_window_extent(r).x0>=fb.x0 and t.get_window_extent(r).x1<=fb.x1)]); pb=ax.patches[1].get_window_extent(r); print("GA side text inside box:",[(t.get_window_extent(r).x0>=pb.x0 and t.get_window_extent(r).x1<=pb.x1) for t in ax.texts if t.get_position()[0]==9.7]); print("GA sub vs side boxes:",[i for i,p in enumerate(ax.patches) if p.get_window_extent(r).overlaps(sb)]); plt.close(fig)
print("ok")
