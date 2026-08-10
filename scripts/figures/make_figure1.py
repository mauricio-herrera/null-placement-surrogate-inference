"""Figure 1 - null placement and paired discordance (PRE style, v2)."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import style_pre as sp

sp.apply()
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "data" / "figure_data"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

d = pd.read_csv(DATA / "figure1_stage3_paired_outcomes.csv")
lk = {(r.index_resolution, r.native_resolution): int(r["count"]) for _, r in d.iterrows()}
n_io, n_no = lk[("reject", "no_reject")], lk[("no_reject", "reject")]
n_bb, n_nn = lk[("reject", "reject")], lk[("no_reject", "no_reject")]

fig = plt.figure(figsize=(7.05, 4.30))

# ---------------- panel (a): the two inferential paths ----------------
axA = fig.add_axes([0.005, 0.52, 0.99, 0.46])
axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
axA.text(0.006, 0.985, "(a)", fontsize=10, fontweight="bold", va="top")

BH = 0.205  # box height


def box(ax, x, y, w, text, fc="white", ec="0.15", lw=0.9, fs=8.0, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, BH,
                                boxstyle="round,pad=0.006,rounding_size=0.013",
                                facecolor=fc, edgecolor=ec, lw=lw))
    ax.text(x + w / 2, y + BH / 2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal")


def arr(ax, x1, x2, y, color="0.15"):
    ax.add_patch(FancyArrowPatch((x1, y + BH / 2), (x2, y + BH / 2),
                                 arrowstyle="-|>", mutation_scale=9,
                                 lw=0.9, color=color, shrinkA=0, shrinkB=0))


yN, yI = 0.575, 0.10
# faint path backgrounds
axA.add_patch(Rectangle((0.215, yN - 0.055), 0.775, BH + 0.20,
                        facecolor=sp.C_FILL_N, edgecolor="none", alpha=0.45, zorder=0))
axA.add_patch(Rectangle((0.215, yI - 0.055), 0.775, BH + 0.20,
                        facecolor=sp.C_FILL_I, edgecolor="none", alpha=0.45, zorder=0))

axA.text(0.222, yN + BH + 0.075, "Native-resolution contract",
         fontsize=8.4, fontweight="bold", color=sp.C_NATIVE, va="bottom")
axA.text(0.222, yI + BH + 0.075, "Index-resolution contract",
         fontsize=8.4, fontweight="bold", color=sp.C_INDEX, va="bottom")

# shared source node, vertically centered between the two rows
ymid = 0.5 * (yN + yI)
box(axA, 0.015, ymid, 0.135, "monthly\nprocess $X$", fc="white", ec="0.15")
# branch connectors
for ytgt, rad in ((yN, -0.25), (yI, 0.25)):
    axA.add_patch(FancyArrowPatch((0.150, ymid + BH / 2), (0.225, ytgt + BH / 2),
                                  arrowstyle="-|>", mutation_scale=9, lw=0.9,
                                  color="0.15",
                                  connectionstyle=f"arc3,rad={rad}"))

# native path: surrogate FIRST (filled = where randomization acts)
box(axA, 0.228, yN, 0.165, "surrogate\n$S_X$", fc=sp.C_NATIVE, ec=sp.C_NATIVE, fs=8.0, bold=True)
axA.texts[-1].set_color("white")
arr(axA, 0.393, 0.425, yN)
box(axA, 0.428, yN, 0.205, "threshold + aggregate\n$\\mathcal{O}$")
arr(axA, 0.633, 0.663, yN)
box(axA, 0.666, yN, 0.150, "detrend\n$D$")
arr(axA, 0.816, 0.846, yN)
box(axA, 0.849, yN, 0.135, r"$\hat\theta$, $p_{\rm nat}$")

# index path: surrogate LAST
box(axA, 0.228, yI, 0.205, "threshold + aggregate\n$\\mathcal{O}$")
arr(axA, 0.433, 0.463, yI)
box(axA, 0.466, yI, 0.150, "detrend\n$D$")
arr(axA, 0.616, 0.646, yI)
box(axA, 0.649, yI, 0.165, "surrogate\n$S_Y$", fc=sp.C_INDEX, ec=sp.C_INDEX, fs=8.0, bold=True)
axA.texts[-1].set_color("white")
arr(axA, 0.814, 0.844, yI)
box(axA, 0.847, yI, 0.137, r"$\hat\theta$, $p_{\rm idx}$")

# ---------------- panel (b): paired 2x2 outcomes ----------------
axB = fig.add_axes([0.075, 0.075, 0.46, 0.375])
axB.set_xlim(0, 2); axB.set_ylim(0, 2)
axB.set_xticks([]); axB.set_yticks([])
[s.set_visible(False) for s in axB.spines.values()]
fig.text(0.012, 0.455, "(b)", fontsize=10, fontweight="bold", va="top")

cells = {  # (row: index outcome 0=reject top, col: native outcome)
    (0, 0): (n_io, sp.C_FILL_I, sp.C_INDEX, "index-only"),
    (0, 1): (n_bb, sp.C_NEUT, "0.25", "both reject"),
    (1, 0): (n_nn, sp.C_NEUT, "0.25", "neither"),
    (1, 1): (n_no, sp.C_FILL_N, sp.C_NATIVE, "native-only"),
}
for (r, c), (n, fc, tc, lab) in cells.items():
    x, y = c, 1 - r
    emph = lab in ("index-only", "native-only")
    axB.add_patch(Rectangle((x + 0.02, y + 0.02), 0.96, 0.96,
                            facecolor=fc, edgecolor="0.25", lw=0.8))
    axB.text(x + 0.5, y + 0.56, f"{n:,}", ha="center", va="center",
             fontsize=13.5 if emph else 10.5,
             fontweight="bold" if emph else "normal",
             color=tc if emph else "0.15")
    axB.text(x + 0.5, y + 0.24, lab, ha="center", va="center",
             fontsize=7.3, color=tc if emph else "0.35",
             style="italic" if not emph else "normal")

axB.text(0.5, 2.07, "native: no reject", ha="center", fontsize=7.8)
axB.text(1.5, 2.07, "native: reject", ha="center", fontsize=7.8)
axB.text(-0.06, 1.5, "index:\nreject", ha="right", va="center", fontsize=7.8)
axB.text(-0.06, 0.5, "index:\nno reject", ha="right", va="center", fontsize=7.8)
axB.text(1.0, -0.16, "13,500 short-memory trajectories, $B=249$",
         ha="center", va="top", fontsize=7.6, color="0.3")

# ---------------- panel (c): the asymmetry, log scale ----------------
axC = fig.add_axes([0.665, 0.115, 0.315, 0.30])
fig.text(0.575, 0.455, "(c)", fontsize=10, fontweight="bold", va="top")
ypos = [1, 0]
vals = [n_io, n_no]
cols = [sp.C_INDEX, sp.C_NATIVE]
labs = ["index-only", "native-only"]
axC.barh(ypos, vals, height=0.58, color=cols, edgecolor="none")
axC.set_xscale("log")
axC.set_xlim(1, 4000)
axC.set_yticks(ypos)
axC.set_yticklabels(labs, fontsize=7.8)
for y, v, c in zip(ypos, vals, cols):
    axC.text(v * 1.25, y, f"{v:,}", va="center", fontsize=8.6,
             fontweight="bold", color=c)
axC.set_xlabel("discordant trajectories (log scale)", fontsize=7.8)
axC.tick_params(axis="y", length=0)
axC.spines["left"].set_visible(False)
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)
axC.tick_params(axis="x", top=False, labeltop=False, bottom=True, labelbottom=True)
axC.tick_params(axis="y", right=False)
axC.xaxis.set_ticks_position("bottom")
axC.text(60, 0.5, r"$68{:}1$", fontsize=11.5, ha="center", va="center",
         color="0.15",
         bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="0.55", lw=0.7))

fig.savefig(OUT / "fig1_null_placement_discordance.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig1_null_placement_discordance.png", dpi=600, bbox_inches="tight")
plt.close(fig)
print("fig1 ok")
