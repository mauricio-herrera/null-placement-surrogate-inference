#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

OUT=Path('/mnt/data/null_transport_manuscript_revised')
OUT.mkdir(parents=True,exist_ok=True)

BLUE='#0072B2'; ORANGE='#D55E00'; GRAY='#777777'; LIGHT='#E6E6E6'; BLACK='#222222'
plt.rcParams.update({'font.size':10,'axes.titlesize':10,'axes.labelsize':10,'legend.fontsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'font.family':'DejaVu Serif'})

# --- Figure 4: Gaussian controls and soft-threshold ablation ---
df=pd.read_csv('/mnt/data/operator_ablation_softpath_B249/operator_ablation_softpath_trajectory_results.csv')

def pooled(cond):
    g=df[(df.condition==cond)&(df.statistic=='ferro_segers')].copy()
    return {
        'oracle':g.rej_oracle.mean()*100,'native':g.rej_native.mean()*100,'index':g.rej_index.mean()*100,
        'io':((g.rej_index)&(~g.rej_native)).mean()*100,
        'no':((g.rej_native)&(~g.rej_index)).mean()*100,
        'sd_io':np.median(g.null_sd_index/g.null_sd_oracle),
        'sd_no':np.median(g.null_sd_native/g.null_sd_oracle),
        'spec':g.index_spectral_error_median.median()
    }

fig,axs=plt.subplots(1,3,figsize=(11.4,3.55))
# panel a
conds=['mean_raw','mean_cubic','count_lam0_cubic','count_lam4_cubic']
labels=['mean\nraw','mean\ncubic','count\n$\\lambda=0$','count\n$\\lambda=4$']
vals=[pooled(c) for c in conds]
x=np.arange(len(conds)); w=.22
axs[0].bar(x-w,[v['oracle'] for v in vals],width=w,color=GRAY,label='oracle')
axs[0].bar(x,[v['native'] for v in vals],width=w,color=BLUE,label='native')
axs[0].bar(x+w,[v['index'] for v in vals],width=w,color=ORANGE,label='index IAAFT')
axs[0].axhline(4.8,ls='--',lw=1,color=BLACK,alpha=.55)
axs[0].set_xticks(x); axs[0].set_xticklabels(labels)
axs[0].set_ylabel('rejection rate (%)'); axs[0].set_ylim(0,27)
axs[0].legend(frameon=False,loc='upper left')
axs[0].text(-.14,1.03,'(a)',transform=axs[0].transAxes,fontweight='bold',fontsize=12)

# panel b: lambda=4 soft path, categorical smooth->hard
steps=[2,1,.5,.25,.1,0]
slabs=['2','1','0.5','0.25','0.1','hard']
vs=[]
for t in steps:
    c='count_lam4_cubic' if t==0 else f'soft_lam4_tau{t:g}_cubic'
    vs.append(pooled(c))
xx=np.arange(len(steps))
axs[1].plot(xx,[v['oracle'] for v in vs],marker='o',color=GRAY,label='oracle')
axs[1].plot(xx,[v['native'] for v in vs],marker='s',color=BLUE,label='native')
axs[1].plot(xx,[v['index'] for v in vs],marker='o',color=ORANGE,label='index IAAFT')
axs[1].axhline(4.8,ls='--',lw=1,color=BLACK,alpha=.55)
axs[1].set_xticks(xx); axs[1].set_xticklabels(slabs)
axs[1].set_xlabel('soft-threshold scale $\\tau$')
axs[1].set_ylabel('rejection rate (%)'); axs[1].set_ylim(0,27)
axs[1].text(-.14,1.03,'(b)',transform=axs[1].transAxes,fontweight='bold',fontsize=12)
axs[1].annotate('hard threshold',xy=(5,vs[-1]['index']),xytext=(3.25,23.5),arrowprops=dict(arrowstyle='->',lw=1,color=BLACK),fontsize=9)

# panel c null dispersion ratios
axs[2].plot(xx,[v['sd_no'] for v in vs],marker='s',color=BLUE,label='native / oracle')
axs[2].plot(xx,[v['sd_io'] for v in vs],marker='o',color=ORANGE,label='index / oracle')
axs[2].axhline(1,ls='--',lw=1,color=BLACK,alpha=.55)
axs[2].set_xticks(xx); axs[2].set_xticklabels(slabs)
axs[2].set_xlabel('soft-threshold scale $\\tau$')
axs[2].set_ylabel('median null-SD ratio'); axs[2].set_ylim(.45,1.08)
axs[2].legend(frameon=False,loc='lower left')
axs[2].text(-.14,1.03,'(c)',transform=axs[2].transAxes,fontweight='bold',fontsize=12)
axs[2].annotate(f"{vs[-1]['sd_io']:.2f}",xy=(5,vs[-1]['sd_io']),xytext=(4.25,.60),arrowprops=dict(arrowstyle='->',lw=1,color=ORANGE),color=ORANGE)

fig.tight_layout(pad=1.0,w_pad=1.7)
fig.savefig(OUT/'fig4_gaussian_operator_ablation.pdf',bbox_inches='tight')
fig.savefig(OUT/'fig4_gaussian_operator_ablation.png',dpi=220,bbox_inches='tight')
plt.close(fig)

# --- Figure 5: within-fiber lumpability diagnostic ---
wf=pd.read_csv('/mnt/data/within_fiber_b249_full/within_fiber_trajectory_results.csv')
features=['ferro_segers','lag1','tdiff_raw','lowfreq_cubic']
flabels=['Ferro--Segers','lag-1 corr.','$T_{\\rm diff}$','low-freq. share']
fig,axs=plt.subplots(1,3,figsize=(11.4,3.55))
# a median distance ratio
x=np.arange(len(features)); w=.34
for j,(lam,col,off) in enumerate([(0,BLUE,-w/2),(4,ORANGE,w/2)]):
    med=[]
    for f in features:
        g=wf[(wf.lambda_value==lam)&(wf.feature==f)]
        med.append(g.distance_ratio.median())
    axs[0].bar(x+off,med,width=w,color=col,label=f'$\\lambda={lam}$')
axs[0].axhline(1,ls='--',lw=1,color=BLACK,alpha=.6)
axs[0].set_xticks(x); axs[0].set_xticklabels(flabels,rotation=18,ha='right')
axs[0].set_ylabel('median fiber / same-$x$ distance')
axs[0].set_ylim(0,17)
axs[0].legend(frameon=False)
axs[0].text(-.14,1.03,'(a)',transform=axs[0].transAxes,fontweight='bold',fontsize=12)
# b fractions fiber > baseline
for j,(lam,col,off) in enumerate([(0,BLUE,-w/2),(4,ORANGE,w/2)]):
    frac=[]
    for f in features:
        g=wf[(wf.lambda_value==lam)&(wf.feature==f)]
        frac.append((g.fiber_wass_scaled>g.same_wass_scaled).mean()*100)
    axs[1].bar(x+off,frac,width=w,color=col,label=f'$\\lambda={lam}$')
axs[1].axhline(50,ls='--',lw=1,color=BLACK,alpha=.55)
axs[1].set_xticks(x); axs[1].set_xticklabels(flabels,rotation=18,ha='right')
axs[1].set_ylabel('fiber distance exceeds baseline (%)')
axs[1].set_ylim(0,105)
axs[1].text(-.14,1.03,'(b)',transform=axs[1].transAxes,fontweight='bold',fontsize=12)
# c lambda=4 spectrum difference vs transported Tdiff distance
sg=wf[(wf.lambda_value==4)&(wf.feature=='tdiff_raw')]
rho,p=spearmanr(sg.native_target_spectral_distance,sg.fiber_wass_scaled)
axs[2].scatter(sg.native_target_spectral_distance,sg.fiber_wass_scaled,s=20,alpha=.72,color=ORANGE,edgecolor='none')
# robust simple least-squares visual only
coef=np.polyfit(sg.native_target_spectral_distance,sg.fiber_wass_scaled,1)
grid=np.linspace(sg.native_target_spectral_distance.min(),sg.native_target_spectral_distance.max(),100)
axs[2].plot(grid,np.polyval(coef,grid),color=BLACK,lw=1)
axs[2].set_xlabel('hidden target-spectrum distance')
axs[2].set_ylabel('within-fiber $T_{\\rm diff}$ distance')
axs[2].text(.04,.92,f'Spearman $\\rho={rho:.2f}$',transform=axs[2].transAxes)
axs[2].text(-.14,1.03,'(c)',transform=axs[2].transAxes,fontweight='bold',fontsize=12)
fig.tight_layout(pad=1.0,w_pad=1.6)
fig.savefig(OUT/'fig5_within_fiber_lumpability.pdf',bbox_inches='tight')
fig.savefig(OUT/'fig5_within_fiber_lumpability.png',dpi=220,bbox_inches='tight')
plt.close(fig)

# supplementary near-zero figure from pooled summary
nz=pd.read_csv('/mnt/data/null_transport_mechanism/nearzero_pooled_summary.csv')
# If tau col parsed numeric, sort descending soft to hard 0
nz=nz.sort_values('tau',ascending=False)
fig,axs=plt.subplots(1,2,figsize=(7.6,3.1))
xx=np.arange(len(nz)); labs=[('hard' if t==0 else f'{t:g}') for t in nz.tau]
axs[0].plot(xx,nz.oracle_reject_rate*100,marker='o',color=GRAY,label='oracle')
axs[0].plot(xx,nz.native_reject_rate*100,marker='s',color=BLUE,label='native')
axs[0].plot(xx,nz.index_reject_rate*100,marker='o',color=ORANGE,label='index IAAFT')
axs[0].axhline(5,ls='--',lw=1,color=BLACK,alpha=.5)
axs[0].set_xticks(xx); axs[0].set_xticklabels(labs,rotation=35)
axs[0].set_xlabel('$\\tau$ (near-zero diagnostic)'); axs[0].set_ylabel('rejection rate (%)')
axs[0].legend(frameon=False,fontsize=8); axs[0].text(-.15,1.03,'(a)',transform=axs[0].transAxes,fontweight='bold')
axs[1].plot(xx,nz.median_sd_ratio_native_oracle,marker='s',color=BLUE,label='native/oracle')
axs[1].plot(xx,nz.median_sd_ratio_index_oracle,marker='o',color=ORANGE,label='index/oracle')
axs[1].axhline(1,ls='--',lw=1,color=BLACK,alpha=.5)
axs[1].set_xticks(xx); axs[1].set_xticklabels(labs,rotation=35)
axs[1].set_xlabel('$\\tau$ (near-zero diagnostic)'); axs[1].set_ylabel('median null-SD ratio')
axs[1].text(-.15,1.03,'(b)',transform=axs[1].transAxes,fontweight='bold')
fig.tight_layout()
fig.savefig(OUT/'figS_nearzero_path.pdf',bbox_inches='tight')
fig.savefig(OUT/'figS_nearzero_path.png',dpi=220,bbox_inches='tight')
plt.close(fig)
print('wrote figures to',OUT)
