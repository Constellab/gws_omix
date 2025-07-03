# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com


import os
from pathlib import Path
from typing import Final, List, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scipy.cluster.hierarchy as sch

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    FloatParam,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    PlotlyResource,
    ResourceSet,
    ShellProxy,
    StrParam,
    TableImporter,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)
from gws_omix.base_env.pydesq2_env_task import Pydesq2ShellProxyHelper


# ------------------------------------------------------------------
@task_decorator(
    "pyDESeq2MultiContrast",
    human_name="pyDESeq2 multi-contrast",
    short_description="DESeq2 contrasts + interactive PCA / volcano / heat-maps",
)
class Pydesq2Multi(Task):

    # ---------- I/O -------------------------------------------------
    input_specs: Final[InputSpecs] = InputSpecs({
        "count_table_file": InputSpec(File, human_name="Counts CSV"),
        "metadata_file":    InputSpec(File, human_name="Metadata TSV"),
        "gtf_file":         InputSpec(File, human_name="Gene annotation GTF"),
    })
    output_specs: Final[OutputSpecs] = OutputSpecs({
        "tables":        OutputSpec(ResourceSet, human_name="DE & summaries"),
        "pca_plots":     OutputSpec(ResourceSet, human_name="Interactive PCA"),
        "heatmap_plots": OutputSpec(ResourceSet, human_name="Interactive heat-maps"),
        "volcano_plots": OutputSpec(ResourceSet, human_name="Interactive volcano"),
    })

    # ---------- user parameters ------------------------------------
    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "samples_colname":      StrParam(short_description=
            "Name of the **metadata** column holding sample IDs "
            "(must match count-table headers)."),
        "genes_colname":        StrParam(short_description=
            "Name of the **count-table** column containing gene identifiers."),
        "condition_column":     StrParam(default_value="Condition", short_description=
            "Metadata column that defines treatments / experimental conditions."),
        "control_condition":    StrParam(short_description=
            "Reference (control) level inside *condition_column*."),
        "timepoint_column":     StrParam(default_value="", short_description=
            "Optional time-point column (per-time-point contrasts if provided)."),
        "group_column":         StrParam(default_value="", short_description=
            "Optional second grouping column (e.g. sex, strain)."),
        "extra_covariates":     StrParam(default_value="", short_description=
            "Comma-separated list of additional covariate columns (added as main effects)."),
        "pvalue_value":         FloatParam(default_value=0.05, min_value=0.0,
            short_description="P-value/FDR threshold for filtered outputs."),
        "log2FoldChange_value": FloatParam(default_value=0.5, short_description=
            "Absolute |log₂FC| cutoff for filtered outputs."),
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_pydeseq2_multicontrast.py",
    )

    # ---------- helper: PCA dropdown -------------------------------
    @staticmethod
    def _pca_dropdown(
        meta: pd.DataFrame, props: pd.DataFrame,
        cond: str, grp: str | None, tp: str | None,
    ) -> go.Figure:

        g_vals: List[str]; t_vals: List[str]

        # symbols per-group
        if grp and grp in meta.columns:
            symbols = ["circle", "square", "diamond", "triangle-up", "x", "cross"]
            grp_syms = {g: symbols[i % len(symbols)]
                        for i, g in enumerate(sorted(meta[grp].unique()))}
            meta["_symbol"] = meta[grp].map(grp_syms)
            g_vals = meta[grp].astype(str).tolist()
        else:
            grp_syms = {"all": "circle"}
            meta["_symbol"] = "circle"
            g_vals = ["all"] * len(meta)

        # marker sizes per-time-point
        if tp and tp in meta.columns:
            tp_sizes = {t: 8 + 4 * i
                        for i, t in enumerate(sorted(meta[tp].unique()))}
            meta["_size"] = meta[tp].map(tp_sizes)
            t_vals = meta[tp].astype(str).tolist()
        else:
            tp_sizes = {"all": 12}
            meta["_size"] = 12
            t_vals = ["all"] * len(meta)

        pal = px.colors.qualitative.Plotly
        cond_cols = {c: pal[i % len(pal)]
                     for i, c in enumerate(sorted(meta[cond].unique()))}

        fig = go.Figure()
        hover_tpl = (
            "Sample: %{customdata[0]}<br>"
            f"{grp or 'Group'}: %{{customdata[1]}}<br>"
            f"{tp  or 'Time'}:  %{{customdata[2]}}<br>"
            "PC1=%{x:.2f} • PC2=%{y:.2f}<extra></extra>"
        )

        for idx, r in meta.iterrows():
            fig.add_trace(
                go.Scatter(
                    x=[r["PC1"]], y=[r["PC2"]],
                    mode="markers",
                    marker=dict(
                        color=cond_cols[r[cond]],
                        symbol=r["_symbol"],
                        size=r["_size"],
                        line=dict(width=0.5, color="#333"),
                    ),
                    customdata=[[idx, g_vals[r.name], t_vals[r.name]]],
                    hovertemplate=hover_tpl, showlegend=False,
                )
            )

        # legend
        for c, col in cond_cols.items():
            fig.add_trace(go.Scatter(x=[None], y=[None],
                                     mode="markers",
                                     marker=dict(color=col, size=10),
                                     name=c))

        grp_vals, tp_vals = ["all"] + list(grp_syms.keys()), ["all"] + list(tp_sizes.keys())

        def _mask(sel_g: str, sel_t: str) -> list[bool]:
            base = [(sel_g in ("all", g_vals[i])) and (sel_t in ("all", t_vals[i]))
                    for i in range(len(meta))]
            return base + [True] * len(cond_cols)

        fig.update_layout(
            updatemenus=[
                dict(buttons=[dict(label=g, method="update",
                                   args=[{"visible": _mask(g, "all")}])
                              for g in grp_vals],
                     direction="down", x=1.01, y=1.15, showactive=True),
                dict(buttons=[dict(label=t, method="update",
                                   args=[{"visible": _mask("all", t)}])
                              for t in tp_vals],
                     direction="down", x=1.18, y=1.15, showactive=True),
            ],
            title=f"PCA — PC1 {props.iloc[0,0]*100:.1f}% • "
                  f"PC2 {props.iloc[0,1]*100:.1f}%",
            xaxis_title="PC1", yaxis_title="PC2",
            margin=dict(l=60, r=200, t=100, b=50),
        )
        return fig

    # ---- static → Plotly converters (unchanged) -------------------
    @staticmethod
    def _heatmap(csv_path: str) -> PlotlyResource:
        df = pd.read_csv(csv_path, index_col=0)
        rl = sch.linkage(df.values, method="average", metric="euclidean")
        cl = sch.linkage(df.values.T, method="average", metric="euclidean")
        df = df.iloc[
            sch.dendrogram(rl, no_plot=True)["leaves"],
            sch.dendrogram(cl, no_plot=True)["leaves"],
        ]
        z = (df.T - df.mean(axis=1)).T
        z = (z.T / z.std(axis=1, ddof=0).replace(0, 1)).T
        vmax = z.abs().to_numpy().max()
        fig = go.Figure(
            go.Heatmap(
                z=z.values,
                x=z.columns,
                y=z.index,
                colorscale="RdBu_r",
                zmin=-vmax, zmax=vmax, zmid=0,
                hovertemplate="<b>%{y}</b><br>%{x}<br>z=%{z:.2f}<extra></extra>",
            )
        )
        fig.update_xaxes(tickangle=45)
        fig.update_yaxes(autorange="reversed")
        pr = PlotlyResource(fig)
        pr.name = f"Heatmap {Path(csv_path).stem}"
        return pr

    @staticmethod
    def _volcano(path: str, gene_col: str, padj: float, fc: float) -> PlotlyResource:
        d = pd.read_csv(path)
        pcol = "padj" if "padj" in d.columns else "pvalue"
        eps = d.loc[d[pcol] > 0, pcol].min() or 1e-300
        d["p_plot"] = -np.log10(d[pcol].replace(0, eps))
        d["signif"] = (d[pcol] < padj) & (d.log2FoldChange.abs() > fc)
        fig = px.scatter(
            d, x="log2FoldChange", y="p_plot", color="signif",
            hover_data=[gene_col, pcol, "log2FoldChange", "gene_name"],
            color_discrete_map={True: "#D62728", False: "#A6ACBC"},
            labels={
                "log2FoldChange": "Log₂ FC",
                "p_plot": f"-log₁₀({pcol})",
                "signif": f"{pcol}<{padj} & |log2FC|>{fc}",
            },
            title=Path(path).stem,
        )
        fig.add_hline(y=-np.log10(padj), line_dash="dot", line_color="grey")
        fig.add_vline(x=fc,  line_dash="dot", line_color="grey")
        fig.add_vline(x=-fc, line_dash="dot", line_color="grey")
        pr = PlotlyResource(fig)
        pr.name = f"Volcano {Path(path).stem}"
        return pr

    # ---------- run ------------------------------------------------
    def run(self, p: ConfigParams, ins: TaskInputs) -> TaskOutputs:

        counts = ins["count_table_file"].path
        meta   = ins["metadata_file"].path
        gtf    = ins["gtf_file"].path

        samp   = p["samples_colname"]
        gene   = p["genes_colname"]
        cond   = p["condition_column"]
        ctrl   = p["control_condition"]
        tp     = p["timepoint_column"] or "None"
        grp    = p["group_column"]     or "None"
        extras = p["extra_covariates"]
        pthr   = p["pvalue_value"]
        fcthr  = p["log2FoldChange_value"]

        shell: ShellProxy = Pydesq2ShellProxyHelper.create_proxy(self.message_dispatcher)
        cmd = (
            "python3 {py} {cts} {meta} "
            "{samp} {gene} {cond} {ctrl} "
            "{tp} {grp} '{extras}' "
            "{pthr} {fcthr} {gtf}"
        ).format(
            py=self.python_file_path, cts=counts, meta=meta,
            samp=samp, gene=gene, cond=cond, ctrl=ctrl,
            tp=tp, grp=grp, extras=extras,
            pthr=pthr, fcthr=fcthr, gtf=gtf,
        )
        if shell.run(cmd, shell_mode=True) != 0:
            raise RuntimeError("DESeq2 analysis failed – check logs")

        work = shell.working_dir
        tables, pca, heat, volc = ResourceSet(), ResourceSet(), ResourceSet(), ResourceSet()

        pca.add_resource(
            PlotlyResource(
                self._pca_dropdown(
                    pd.read_csv(f"{work}/pca_metadata.csv"),
                    pd.read_csv(f"{work}/pca_proportions.csv"),
                    cond,
                    None if grp == "None" else grp,
                    None if tp  == "None" else tp,
                )
            ),
            "Interactive PCA",
        )

        for fn in os.listdir(work):
            fp = os.path.join(work, fn)

            if fn.startswith("DE_") and fn.endswith(".csv"):
                tables.add_resource(
                    TableImporter.call(
                        File(fp),
                        {"delimiter": ",", "header": 0,
                         "file_format": "csv", "index_column": 0},
                    ), fn,
                )
                volc.add_resource(
                    self._volcano(fp, gene, pthr, fcthr),
                    f"Volcano {Path(fp).stem}",
                )
            elif fn.startswith("Summary_") and fn.endswith(".csv"):
                tables.add_resource(
                    TableImporter.call(
                        File(fp),
                        {"delimiter": ",", "header": 0, "file_format": "csv"},
                    ), fn,
                )
            elif fn.startswith("Heatmap_") and fn.endswith(".csv"):
                heat.add_resource(self._heatmap(fp), f"Heatmap {Path(fp).stem}")

        return {
            "tables":        tables,
            "pca_plots":     pca,
            "heatmap_plots": heat,
            "volcano_plots": volc,
        }
