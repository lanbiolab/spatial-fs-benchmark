"""Create the exploratory batch-mixing diagnostic; requires pandas and matplotlib."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

root = Path(__file__).resolve().parent
d = pd.read_csv(root / "results/control_component_diagnostics.tsv", sep="\t")
d = d[d.component.eq("BatchMixingMean")]
datasets = ["DLPFC", "MouseBrainSerialSections", "STOmics0212", "STOmics0218", "STOmics0224", "E8p5Embryo", "E9p5Embryo"]
pretty = ["DLPFC", "Mouse brain", "STOmics 0212", "STOmics 0218", "STOmics 0224", "E8.5 embryo", "E9.5 embryo"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.2), sharey=True, constrained_layout=True)
for ax, integrator, title in zip(axes, ["scvi", "cellcharter"], ["scVI", "CellCharter"]):
    g = d[d.integration_method.eq(integrator)]
    for j, dataset in enumerate(datasets):
        sub = g[g.dataset.eq(dataset)].set_index("control")
        assert len(sub) == 2
        ax.scatter(sub.loc["random", "competitive_median"], j, marker="|", s=200, color="#555555", linewidths=2, zorder=2)
        ax.scatter(sub.loc["all_features", "control_value"], j, marker="s", s=33, color="#d88734", zorder=3)
        ax.scatter(sub.loc["random", "control_value"], j, marker="o", s=39, color="#275a99", zorder=3)
    ax.set_title(title, weight="bold", loc="left")
    ax.set_xlim(-0.02, 1.02)
    ax.set_xticks([0, .25, .5, .75, 1])
    ax.grid(axis="x", color="#dddddd", linewidth=.6)
    ax.set_axisbelow(True)
    ax.set_yticks(range(len(datasets)), pretty)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Batch-mixing component (frozen scaled score)")
axes[0].invert_yaxis()
handles = [Line2D([0],[0],marker="|",linestyle="",color="#555555",markersize=14,markeredgewidth=2,label="Median of 29 competitive methods"),
           Line2D([0],[0],marker="o",linestyle="",color="#275a99",markersize=6,label="Random (500 genes)"),
           Line2D([0],[0],marker="s",linestyle="",color="#d88734",markersize=6,label="All features")]
fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=8)
fig.suptitle("Batch mixing: controls versus the competitive panel", fontweight="bold", x=.05, ha="left")
fig.savefig(root / "results/batch_mixing_control_diagnostic.png", dpi=200, bbox_inches="tight")
