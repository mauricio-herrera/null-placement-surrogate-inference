import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

IN='/mnt/data/SINCA_PM25_AOAS_RESULTS_V1_FINAL/SINCA_station_results_V1.csv'
OUT=Path('/mnt/data/final_manuscript_work/NULL_TRANSPORT_REVISED_PACKAGE_V1')
df=pd.read_csv(IN)
fixed=df[df.branch=='fixed'].copy().sort_values('station_id')
eq=df[df.branch=='equalized'].copy().sort_values('station_id')

plt.rcParams.update({
    'font.family':'serif',
    'mathtext.fontset':'dejavuserif',
    'font.size':9.5,
    'axes.labelsize':10.5,
    'axes.titlesize':10.5,
    'legend.fontsize':9,
    'xtick.labelsize':9,
    'ytick.labelsize':9,
    'figure.dpi':160,
    'savefig.bbox':'tight',
    'pdf.fonttype':42,
    'ps.fonttype':42,
})

# Main 3-panel applied figure
fig, axes = plt.subplots(1,3,figsize=(13.4,4.05))
for ax,g,title in [
    (axes[0],fixed,r'(a) Fixed 24-h standard level ($50\,\mu$g m$^{-3}$)'),
    (axes[1],eq,r'(b) Probability-equalized thresholds')]:
    ax.scatter(g.p_native_fs90,g.p_index_fs90,s=31,alpha=.9)
    lim=(0.0017,1.05)
    ax.plot(lim,lim,'--',lw=1,color='0.35')
    ax.axhline(.05,ls=':',lw=1,color='0.45')
    ax.axvline(.05,ls=':',lw=1,color='0.45')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.set_xlabel('native transported-null $p$')
    ax.set_ylabel('index-level IAAFT $p$')
    ax.set_title(title,loc='left',fontweight='bold')
    io=int(np.sum(g.paired_outcome=='index_only'))
    no=int(np.sum(g.paired_outcome=='native_only'))
    bo=int(np.sum(g.paired_outcome=='both'))
    ne=int(np.sum(g.paired_outcome=='neither'))
    ax.text(.03,.97,f'index-only={io}\nnative-only={no}\nboth={bo}\nneither={ne}',
            transform=ax.transAxes,va='top',ha='left',fontsize=8.7,
            bbox=dict(boxstyle='round,pad=.3',fc='white',ec='0.75',alpha=.9))

# Panel c: center-gap intervention
ax=axes[2]
fix=fixed.set_index('station_id')
equ=eq.set_index('station_id').loc[fix.index]
for sid in fix.index:
    ax.plot([0,1],[fix.loc[sid,'center_gap_fs90'],equ.loc[sid,'center_gap_fs90']],
            color='0.75',lw=.8,zorder=1)
ax.scatter(np.zeros(len(fix)),fix.center_gap_fs90,s=25,zorder=2,label='fixed')
ax.scatter(np.ones(len(equ)),equ.center_gap_fs90,s=25,zorder=2,label='equalized')
med0=float(fix.center_gap_fs90.median()); med1=float(equ.center_gap_fs90.median())
ax.plot([0,1],[med0,med1],color='black',lw=2.2,zorder=3)
ax.scatter([0,1],[med0,med1],color='black',s=40,zorder=4)
ax.set_xticks([0,1],['fixed 50','equalized'])
ax.set_xlim(-.25,1.25)
ax.set_ylabel(r'null-center gap $G_s$')
ax.set_title('(c) Seasonal equalization intervention',loc='left',fontweight='bold')
ax.text(.03,.97,f'median: {med0:.3f} $\\rightarrow$ {med1:.3f}\ndecrease: 31/31 stations',
        transform=ax.transAxes,va='top',ha='left',fontsize=8.7,
        bbox=dict(boxstyle='round,pad=.3',fc='white',ec='0.75',alpha=.9))
fig.tight_layout(w_pad=2.0)
fig.savefig(OUT/'fig6_sinca_application.pdf')
fig.savefig(OUT/'fig6_sinca_application.png',dpi=220)
plt.close(fig)

# Supplement: prespecified width mechanism check
fig,ax=plt.subplots(figsize=(5.8,4.45))
y=np.log(fixed.width_ratio_fs90.to_numpy(float))
x=fixed.neff_frozen.to_numpy(float)
ax.scatter(x,y,s=34)
ax.axhline(0,ls='--',lw=1,color='0.4')
ax.set_xlabel(r'$N_{\mathrm{eff},p}$')
ax.set_ylabel(r'$\log\{\mathrm{SD}_{index}/\mathrm{SD}_{native}\}$')
ax.set_title('SINCA prespecified null-width mechanism check',fontweight='bold')
ax.text(.03,.97,r'$\rho_S=-0.379$'+'\nregion bootstrap 95% CI\n[-0.624, 0.136]',
        transform=ax.transAxes,va='top',ha='left',fontsize=9,
        bbox=dict(boxstyle='round,pad=.3',fc='white',ec='0.75',alpha=.9))
fig.tight_layout()
fig.savefig(OUT/'figS2_sinca_width_mechanism.pdf')
fig.savefig(OUT/'figS2_sinca_width_mechanism.png',dpi=220)
plt.close(fig)
print('created figures')
