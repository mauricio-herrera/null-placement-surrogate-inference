"""Shared PRE house style for Paper A figures."""
import matplotlib as mpl

# Okabe-Ito colorblind-safe palette
C_INDEX  = "#D55E00"   # vermillion: index-resolution contract
C_NATIVE = "#0072B2"   # blue: native-resolution contract
C_GRAY   = "#555555"
C_LIGHT  = "#BBBBBB"
C_FILL_I = "#F6DCC8"   # light vermillion fill
C_FILL_N = "#CFE3F0"   # light blue fill
C_NEUT   = "#EFEFEF"

def apply():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 8.0,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9.0,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.size": 1.7,
        "ytick.minor.size": 1.7,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "lines.linewidth": 1.1,
        "lines.markersize": 4.5,
        "legend.frameon": False,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
    })
