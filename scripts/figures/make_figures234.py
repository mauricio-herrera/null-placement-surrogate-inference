"""Figures 2-4 for Paper A (PRE style, v2)."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import style_pre as sp

sp.apply()
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "data" / "figure_data"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

ATTAIN_B249 = 100 * 12 / 250  # 4.8% attainable at B=249, strict p<0.05

# ================= FIGURE 2 : seasonal concentration =================
d2 = pd.read_csv(DATA / "figure2_stage3_concentration.csv").sort_values("lambda")
x = d2["median_neff_p"].to_numpy()
mis = 100 * d2["mechanistic_misattribution_rate"].to_numpy()
lo = 100 * d2["ci_low"].to_numpy()
hi = 100 * d2["ci_high"].to_numpy()
nat = 100 * d2["native_rejection_rate"].to_numpy()
ann = 100 * d2["annual_rejection_rate"].to_numpy()
lam = d2["lambda"].to_numpy()

fig, ax = plt.subplots(figsize=(3.40, 2.95))
fig.subplots_adjust(left=0.155, right=0.965, bottom=0.155, top=0.845)

# attainable Monte Carlo level
ax.axhline(ATTAIN_B249, color="0.55", lw=0.7, ls=(0, (4, 3)), zorder=1)
ax.text(11.15, ATTAIN_B249 - 1.25, "native-test attainable size (4.8%)",
        fontsize=6.6, color=sp.C_NATIVE, ha="left")

# native rejection: calibrated, flat
ax.plot(x, nat, marker="s", ms=3.6, color=sp.C_NATIVE, lw=1.0, zorder=3,
        label="native rejection")
# index-only discordance: the effect
ax.fill_between(x, lo, hi, color=sp.C_INDEX, alpha=0.18, lw=0, zorder=2)
ax.plot(x, mis, marker="o", color=sp.C_INDEX, lw=1.2, zorder=4,
        label="index-only discordance")

ax.invert_xaxis()
ax.set_xlabel(r"effective seasonal probability support $N_{\mathrm{eff},p}$")
ax.set_ylabel("rate (%)")
ax.set_ylim(0, 19.5)
ax.set_yticks([0, 5, 10, 15])
ax.set_xlim(12.55, 4.35)

# secondary top axis: lambda
sec = ax.secondary_xaxis("top")
sec.set_xticks(x)
sec.set_xticklabels([f"{v:g}" for v in lam], fontsize=7.0)
sec.set_xlabel(r"seasonal concentration $\lambda$", fontsize=7.8)
sec.tick_params(direction="in", length=3.0, width=0.6)
ax.tick_params(top=False)

ax.annotate("", xy=(0.97, 0.085), xycoords="axes fraction",
            xytext=(0.62, 0.085), textcoords="axes fraction",
            arrowprops=dict(arrowstyle="-|>", color="0.35", lw=0.7))
ax.text(0.795, 0.115, "concentration increases", transform=ax.transAxes,
        fontsize=6.8, color="0.35", ha="center")
ax.legend(loc="upper left", handlelength=1.6, borderaxespad=0.4)

fig.savefig(OUT / "fig2_seasonal_concentration.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig2_seasonal_concentration.png", dpi=600, bbox_inches="tight")
plt.close(fig)
print("fig2 ok")

# ================= FIGURE 3 : probability equalization =================
d3 = pd.read_csv(DATA / "figure3_stage5_mitigation.csv")
order = ["absolute_lambda4", "oracle_percentile", "estimated_percentile_15y",
         "estimated_percentile_30y", "estimated_percentile_60y"]
d3 = d3.set_index("scenario").loc[order].reset_index()
labels = ["concentrated\nprofile\n($\\lambda=4$)", "oracle\nuniform",
          "15-y\npercentile", "30-y\npercentile", "60-y\npercentile"]
X = np.arange(5)
mis3 = 100 * d3["mechanistic_misattribution_rate"].to_numpy()
lo3 = 100 * d3["misattrib_ci_low"].to_numpy()
hi3 = 100 * d3["misattrib_ci_high"].to_numpy()
neff3 = d3["median_true_p_neff"].to_numpy()

fig, (axT, axB) = plt.subplots(
    2, 1, figsize=(3.40, 3.55), sharex=True,
    gridspec_kw=dict(height_ratios=[2.1, 1.0], hspace=0.10))
fig.subplots_adjust(left=0.155, right=0.965, bottom=0.185, top=0.965)

# group shading: equalized constructions
for ax in (axT, axB):
    ax.axvspan(0.55, 4.55, color="0.94", zorder=0)

cols3 = [sp.C_INDEX] + [sp.C_NATIVE] * 4
axT.errorbar(X, mis3, yerr=[mis3 - lo3, hi3 - mis3], fmt="none",
             ecolor="0.35", elinewidth=0.9, capsize=2.4, zorder=3)
axT.scatter(X, mis3, s=34, c=cols3, zorder=4, edgecolor="white", lw=0.6)
axT.set_ylabel("index-only\ndiscordance (%)")
axT.set_ylim(0, 20.0)
axT.set_yticks([0, 5, 10, 15, 20])

# reduction arrow
axT.annotate("", xy=(1.0, mis3[1] + 1.1), xytext=(1.0, mis3[0] - 0.4),
             arrowprops=dict(arrowstyle="-|>", color="0.25", lw=0.9))
axT.text(1.17, 0.5 * (mis3[0] + mis3[1]) + 0.6, r"$-75\%$",
         fontsize=8.6, color="0.15", fontweight="bold")
axT.text(2.55, 17.6, "equalized event probability", fontsize=7.0,
         color="0.35", ha="center")

# bottom: threshold-probability geometry
axB.bar(X, neff3, width=0.55, color=["0.65"] + ["0.45"] * 4, edgecolor="none")
axB.axhline(12, color=sp.C_NATIVE, lw=0.8, ls=(0, (4, 3)))
axB.text(2.5, 13.55, "uniform limit $N_{\\mathrm{eff},p}=12$",
         fontsize=6.6, color=sp.C_NATIVE, ha="center", va="center",
         bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.85))
for xi, v in zip(X, neff3):
    axB.text(xi, v + 0.45, f"{v:.1f}", ha="center", fontsize=6.6, color="0.25")
axB.set_ylim(0, 15.4)
axB.set_yticks([0, 6, 12])
axB.set_ylabel("$N_{\\mathrm{eff},p}$")
axB.set_xticks(X)
axB.set_xticklabels(labels, fontsize=6.9)
axB.tick_params(axis="x", length=0)

axT.text(0.012, 0.965, "(a)", transform=axT.transAxes, fontsize=9.5,
         fontweight="bold", va="top")
axB.text(0.012, 0.93, "(b)", transform=axB.transAxes, fontsize=9.5,
         fontweight="bold", va="top")

fig.savefig(OUT / "fig3_probability_equalization.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig3_probability_equalization.png", dpi=600, bbox_inches="tight")
plt.close(fig)
print("fig3 ok")

# ================= FIGURE 4 : mechanism-specific power =================
d4 = pd.read_csv(DATA / "figure4_stage4_power.csv")
fig, ax = plt.subplots(figsize=(3.40, 2.95))
fig.subplots_adjust(left=0.155, right=0.965, bottom=0.155, top=0.965)

# reference lines
ax.axhline(50, color="0.60", lw=0.7, ls=(0, (4, 3)))
ax.text(1.045, 51.3, "50% power", fontsize=6.8, color="0.45", ha="right")
ax.axhline(4.0, color="0.60", lw=0.7, ls=(0, (1.5, 2.2)))
ax.text(1.045, 1.6, "null gate rate (4.0%)", fontsize=6.8, color="0.45", ha="right")

specs = [("history_feedback", "o", sp.C_INDEX, "history feedback"),
         ("persistent_regime_alignment", "s", sp.C_NATIVE, "persistent regime")]
for mech, mk, col, lab in specs:
    dd = d4[d4["curve_mechanism"] == mech].sort_values("mechanism_level")
    xx = dd["mechanism_level"].to_numpy()
    yy = 100 * dd["full_gate_power"].to_numpy()
    lo = 100 * dd["full_gate_power_ci_low"].to_numpy()
    hi = 100 * dd["full_gate_power_ci_high"].to_numpy()
    ax.fill_between(xx, lo, hi, color=col, alpha=0.15, lw=0)
    ax.plot(xx, yy, marker=mk, color=col, lw=1.2, label=lab, ms=4.2)

ax.annotate(r"$30.7\%\ [23.8,38.5]$",
            xy=(1.0, 30.7), xytext=(0.585, 41.5), fontsize=7.4, color=sp.C_INDEX,
            arrowprops=dict(arrowstyle="-", color=sp.C_INDEX, lw=0.6,
                            shrinkA=2, shrinkB=3))
ax.set_xlabel(r"injection strength $\kappa$")
ax.set_ylabel("full-gate detection rate (%)")
ax.set_xlim(-0.04, 1.06)
ax.set_ylim(0, 58)
ax.legend(loc="upper left", handlelength=1.7, borderaxespad=0.4)

fig.savefig(OUT / "fig4_mechanism_specific_power.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig4_mechanism_specific_power.png", dpi=600, bbox_inches="tight")
plt.close(fig)
print("fig4 ok")
