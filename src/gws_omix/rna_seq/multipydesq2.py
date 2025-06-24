#!/usr/bin/env python3
# PyDESeq2 multi-contrast wrapper – PCA interactive + named resources
# v6.1 (2025-08-06)

import os
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import scipy.cluster.hierarchy as sch

from gws_core import (
    ConfigParams, ConfigSpecs, File, FloatParam, InputSpec, InputSpecs,
    OutputSpec, OutputSpecs, PlotlyResource, ResourceSet,
    ShellProxy, StrParam, TableImporter, Task, TaskInputs, TaskOutputs,
    task_decorator,
)
from gws_omix.base_env.pydesq2_env_task import Pydesq2ShellProxyHelper


# ────────────────────────────────────────────────────────────────────
@task_decorator(
    "pyDESeq2MultiContrast",
    human_name="pyDESeq2 multi-contrast",
    short_description="DESeq2 contrasts + interactive PCA / volcano / heat-maps",
)
class Pydesq2Multi(Task):
    input_specs: Final[InputSpecs] = InputSpecs({
        "count_table_file": InputSpec(File, human_name="Counts CSV"),
        "metadata_file":    InputSpec(File, human_name="Metadata TSV"),
        "gtf_file":         InputSpec(File, human_name="Gene annotation GTF"),
    })
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "tables":        OutputSpec(ResourceSet, human_name="DE + summary tables"),
        "pca_plots":     OutputSpec(ResourceSet, human_name="Interactive PCA"),
        "heatmap_plots": OutputSpec(ResourceSet, human_name="Interactive heat-maps"),
        "volcano_plots": OutputSpec(ResourceSet, human_name="Interactive volcano"),
    })
    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "genes_colname":        StrParam(),
        "control_condition":    StrParam(),
        "pvalue_value":         FloatParam(default_value=0.05, min_value=0.0),
        "log2FoldChange_value": FloatParam(default_value=0.5),
    })
    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_pydeseq2_multicontrast.py",
    )

    # 1) PCA with two dropdown filters (Group, Timepoint) -------------
    @staticmethod
    def _pca_dropdown(meta: pd.DataFrame, props: pd.DataFrame) -> go.Figure:
        # ----- size ↔ timepoint --------------------------------------
        if "Timepoint" in meta.columns:
            tp_sizes = {tp: 8 + 4*i for i, tp in enumerate(sorted(meta["Timepoint"].unique()))}
            meta["_size"] = meta["Timepoint"].map(tp_sizes)
        else:
            tp_sizes, meta["_size"] = {"all": 12}, 12

        # ----- symbol ↔ group ---------------------------------------
        if "Group" in meta.columns:
            symbols = ["circle", "square", "diamond", "triangle-up", "x", "cross"]
            grp_symbols = {g: symbols[i % len(symbols)]
                           for i, g in enumerate(sorted(meta["Group"].unique()))}
            meta["_symbol"] = meta["Group"].map(grp_symbols)
        else:
            grp_symbols, meta["_symbol"] = {"all": "circle"}, "circle"

        # ----- colour ↔ condition -----------------------------------
        colour_seq = px.colors.qualitative.Plotly
        cond_colours = {c: colour_seq[i % len(colour_seq)]
                        for i, c in enumerate(sorted(meta["Condition"].unique()))}

        # ----- build traces  (one per sample) -----------------------
        fig = go.Figure()
        for _, r in meta.iterrows():
            fig.add_trace(go.Scatter(
                x=[r["PC1"]], y=[r["PC2"]], mode="markers",
                marker=dict(color=cond_colours[r["Condition"]],
                            symbol=r["_symbol"], size=r["_size"],
                            line=dict(width=0.5, color="#333")),
                name=r["Condition"], legendgroup=r["Condition"], showlegend=False,
                customdata=[[r.get("Sample", r.name),
                             r.get("Group", ""), r.get("Timepoint", "")]],
                hovertemplate=("Sample: %{customdata[0]}<br>"
                               "Group: %{customdata[1]}<br>"
                               "TP: %{customdata[2]}<br>"
                               "PC1=%{x:.2f} · PC2=%{y:.2f}<extra></extra>")
            ))

        # dummy per condition for legend
        for cond, col in cond_colours.items():
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(color=col, size=10),
                name=cond, legendgroup=cond, showlegend=True,
                legendgrouptitle_text="Condition"))

        # dropdown menus --------------------------------------------
        grp_vals = ["all"] + list(grp_symbols.keys())
        tp_vals  = ["all"] + list(tp_sizes.keys())

        def _mask(sel_grp, sel_tp):
            m = [(sel_grp in ("all", g) and sel_tp in ("all", t))
                 for g, t in zip(meta.get("Group", ["all"]*len(meta)),
                                 meta.get("Timepoint", ["all"]*len(meta)))]
            m += [True] * len(cond_colours)   # keep condition dummies visible
            return m

        updatemenus = [
            dict(buttons=[dict(label=v, method="update",
                               args=[{"visible": _mask(v, "all")}])
                          for v in grp_vals],
                 x=1.01, y=1.15, direction="down",
                 pad={"r": 10}, showactive=True),
            dict(buttons=[dict(label=v, method="update",
                               args=[{"visible": _mask("all", v)}])
                          for v in tp_vals],
                 x=1.18, y=1.15, direction="down",
                 pad={"r": 10}, showactive=True),
        ]

        fig.update_layout(
            updatemenus=updatemenus,
            title=f"PCA — PC1 {props.iloc[0,0]*100:.1f}% · PC2 {props.iloc[0,1]*100:.1f}%",
            xaxis_title="PC1", yaxis_title="PC2",
            legend_title_text="Condition",
            margin=dict(l=60, r=200, t=100, b=50),
        )
        return fig

    @staticmethod
    def plot_pca(meta_csv: str, prop_csv: str) -> PlotlyResource:
        meta  = pd.read_csv(meta_csv)
        props = pd.read_csv(prop_csv)
        res   = PlotlyResource(Pydesq2Multi._pca_dropdown(meta, props))
        res.name = "Interactive PCA"
        return res

    # 2) heat-map ----------------------------------------------------
    @staticmethod
    def plot_heatmap(csv_path: str) -> PlotlyResource:
        df = pd.read_csv(csv_path, index_col=0)
        r_link = sch.linkage(df.values, method="average", metric="euclidean")
        c_link = sch.linkage(df.values.T, method="average", metric="euclidean")
        df = df.iloc[sch.dendrogram(r_link, no_plot=True)["leaves"],
                     sch.dendrogram(c_link, no_plot=True)["leaves"]]

        z = (df.T - df.mean(axis=1)).T
        z = (z.T / z.std(axis=1, ddof=0).replace(0, 1)).T
        vmax = z.abs().to_numpy().max()
        h, w = max(600, 20*z.shape[0]), max(900, 50*z.shape[1])

        fig = go.Figure(go.Heatmap(
            z=z.values, x=z.columns, y=z.index,
            colorscale="RdBu_r", zmin=-vmax, zmax=vmax, zmid=0,
            colorbar=dict(title="z-score (VST)", thickness=18,
                          lenmode="pixels", len=200),
            hovertemplate="<b>%{y}</b><br>%{x}<br>z=%{z:.2f}<extra></extra>"))
        fig.update_xaxes(tickangle=45, tickfont=dict(size=10))
        fig.update_yaxes(tickfont=dict(size=10))
        fig.update_layout(title=Path(csv_path).stem,
                          width=w, height=h,
                          margin=dict(l=160, r=40, t=60, b=120))

        res = PlotlyResource(fig)
        res.name = f"Interactive Heatmap {Path(csv_path).stem}"
        return res

    # 3) volcano -----------------------------------------------------
    @staticmethod
    def plot_volcano(path: str, gene_col: str,
                     *, padj_thr: float, fc_thr: float) -> PlotlyResource:
        df   = pd.read_csv(path)
        pcol = "padj" if "padj" in df.columns else "pvalue"
        eps  = df.loc[df[pcol] > 0, pcol].min() or 1e-300
        df["p_plot"] = -np.log10(df[pcol].replace(0, eps))
        df["signif"] = (df[pcol] < padj_thr) & (df.log2FoldChange.abs() > fc_thr)

        fig = px.scatter(
            df, x="log2FoldChange", y="p_plot", color="signif",
            hover_data=[gene_col, pcol, "log2FoldChange", "gene_name",
                        *[c for c in ("Timepoint", "Group") if c in df.columns]],
            color_discrete_map={True: "#D62728", False: "#A6ACBC"},
            labels={"log2FoldChange": "Log₂ FC",
                    "p_plot": f"-log₁₀({pcol})",
                    "signif": f"{pcol}<{padj_thr} & |log2FC|>{fc_thr}"},
            title=Path(path).stem)
        fig.add_hline(y=-np.log10(padj_thr), line_dash="dot", line_color="grey")
        fig.add_vline(x= fc_thr,  line_dash="dot", line_color="grey")
        fig.add_vline(x=-fc_thr, line_dash="dot", line_color="grey")

        res = PlotlyResource(fig)
        res.name = f"Interactive Volcano {Path(path).stem}"
        return res

    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        counts = inputs["count_table_file"].path
        meta   = inputs["metadata_file"].path
        gtf    = inputs["gtf_file"].path
        gcol, ctrl = params["genes_colname"], params["control_condition"]
        p_thr, fc_thr = params["pvalue_value"], params["log2FoldChange_value"]

        shell: ShellProxy = Pydesq2ShellProxyHelper.create_proxy(self.message_dispatcher)
        cmd = f"python3 {self.python_file_path} {counts} {meta} {gcol} {ctrl} {p_thr} {fc_thr} {gtf}"
        if shell.run(cmd, shell_mode=True) != 0:
            raise RuntimeError("DESeq2 analysis failed – check logs")

        work = shell.working_dir
        tables, pca, heat, volc = ResourceSet(), ResourceSet(), ResourceSet(), ResourceSet()

        pca.add_resource(self.plot_pca(f"{work}/pca_metadata.csv",
                                       f"{work}/pca_proportions.csv"),
                         "Interactive PCA")

        for fn in os.listdir(work):
            fp = os.path.join(work, fn)
            if fn.startswith("DE_") and fn.endswith(".csv"):
                tables.add_resource(
                    TableImporter.call(File(fp),
                                       {"delimiter": ",", "header": 0,
                                        "file_format": "csv", "index_column": 0}), fn)
                volc.add_resource(
                    self.plot_volcano(fp, gcol, padj_thr=p_thr, fc_thr=fc_thr),
                    f"Volcano {Path(fp).stem}")
            elif fn.startswith("Summary_") and fn.endswith(".csv"):
                tables.add_resource(
                    TableImporter.call(File(fp),
                                       {"delimiter": ",", "header": 0, "file_format": "csv"}), fn)
            elif fn.startswith("Heatmap_") and fn.endswith(".csv"):
                heat.add_resource(
                    self.plot_heatmap(fp),
                    f"Heatmap {Path(fp).stem}")

        return {
            "tables":        tables,
            "pca_plots":     pca,
            "heatmap_plots": heat,
            "volcano_plots": volc,
        }
