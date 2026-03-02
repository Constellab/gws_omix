# LICENSE
# This software is the exclusive property of Gencovery SAS.
# The use and distribution of this software is prohibited without the prior consent of Gencovery SAS.
# About us: https://gencovery.com

import os
import shlex
from pathlib import Path
from typing import Final

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


def _split_csv(s: str) -> list[str]:
    return [z.strip() for z in str(s or "").split(",") if z.strip()]


@task_decorator(
    "pyDESeq2MultiContrast",
    human_name="pyDESeq2 multi-contrast",
    short_description="DESeq2 contrasts + interactive PCA / volcano / heat-maps",
)
class Pydesq2Multi(Task):
    """
    pyDESeq2 multi-contrast differential expression task.

    Runs DESeq2-style (PyDESeq2) differential expression for each treatment vs the chosen control,
    plus an optional pooled “ALL vs CTRL” comparison. Covariates found in the metadata (timepoint,
    group, and any user-provided extra covariates) are included in the design when available.

    Outputs per contrast:
    - FULL (unfiltered) results table and a FILTERED table (padj if available, else pvalue; + |log2FC| cutoff)
    - interactive volcano plot (log2FC vs −log10(padj/pvalue))
    - VST-based heatmap (top genes by |log2FC|, up/down balanced)

    Also outputs a global VST-based PCA of all samples.
    full documentation is available in community : https://constellab.community/bricks/gws_omix/latest/doc/use-cases/rna-seq-splice-aware-aligner-and-quasi-mapper
    """

    input_specs: Final[InputSpecs] = InputSpecs({
        "count_table_file": InputSpec(File, human_name="Counts CSV"),
        "metadata_file":    InputSpec(File, human_name="Metadata TSV"),
        "gtf_file":         InputSpec(File, human_name="Gene annotation GTF"),
    })

    output_specs: Final[OutputSpecs] = OutputSpecs({
        "tables_full":    OutputSpec(ResourceSet, human_name="DE results (FULL, unfiltered)"),
        "tables_filtered": OutputSpec(ResourceSet, human_name="DE results (FILTERED)"),
        "pca_plots":      OutputSpec(ResourceSet, human_name="Interactive PCA"),
        "heatmap_plots":  OutputSpec(ResourceSet, human_name="Interactive heat-maps"),
        "volcano_plots":  OutputSpec(ResourceSet, human_name="Interactive volcano"),
    })

    config_specs: Final[ConfigSpecs] = ConfigSpecs({
        "samples_colname": StrParam(short_description=
            "Name of the **metadata** column holding sample IDs (must match count-table headers)."),
        "genes_colname": StrParam(short_description=
            "Name of the **count-table** column containing feature identifiers (must be UNIQUE)."),
        "condition_column": StrParam(default_value="Condition", short_description=
            "Metadata column that defines treatments / experimental conditions."),
        "control_condition": StrParam(short_description=
            "Reference (control) level inside *condition_column*."),
        "timepoint_column": StrParam(default_value="", short_description=
            "Optional time-point column (per-time-point contrasts if provided)."),
        "group_column": StrParam(default_value="", short_description=
            "Optional second grouping column (e.g. sex, strain)."),
        "extra_covariates": StrParam(default_value="", short_description=
            "Comma-separated list of additional covariate columns (added as main effects)."),
        "padj_value": FloatParam(default_value=0.05, min_value=0.0,
            short_description="Threshold used for FILTERED outputs (uses padj if available else pvalue)."),
        "log2FoldChange_value": FloatParam(default_value=0.5, short_description=
            "Absolute |log₂FC| cutoff for FILTERED outputs."),

        "annotation_columns": StrParam(default_value="", short_description=
            "Comma-separated extra columns for outputs/plots (e.g. gene_name,gene_biotype). "
            "If not present in the counts table, they will be fetched from the GTF when possible."),
    })

    python_file_path: Final[str] = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_pydeseq2_multicontrast.py",
    )

    @staticmethod
    def _pca_dropdown(meta: pd.DataFrame, props: pd.DataFrame,
                      cond: str, grp: str | None, tp: str | None) -> go.Figure:

        if grp and grp in meta.columns:
            symbols = ["circle", "square", "diamond", "triangle-up", "x", "cross"]
            grp_syms = {g: symbols[i % len(symbols)]
                        for i, g in enumerate(sorted(meta[grp].astype(str).unique()))}
            meta["_symbol"] = meta[grp].astype(str).map(grp_syms)
            g_vals = meta[grp].astype(str).tolist()
        else:
            grp_syms = {"all": "circle"}
            meta["_symbol"] = "circle"
            g_vals = ["all"] * len(meta)

        if tp and tp in meta.columns:
            tp_sizes = {t: 8 + 4 * i
                        for i, t in enumerate(sorted(meta[tp].astype(str).unique()))}
            meta["_size"] = meta[tp].astype(str).map(tp_sizes)
            t_vals = meta[tp].astype(str).tolist()
        else:
            tp_sizes = {"all": 12}
            meta["_size"] = 12
            t_vals = ["all"] * len(meta)

        pal = px.colors.qualitative.Plotly
        cond_cols = {c: pal[i % len(pal)]
                     for i, c in enumerate(sorted(meta[cond].astype(str).unique()))}

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
                        color=cond_cols[str(r[cond])],
                        symbol=r["_symbol"],
                        size=r["_size"],
                        line=dict(width=0.5, color="#333"),
                    ),
                    customdata=[[idx, g_vals[r.name], t_vals[r.name]]],
                    hovertemplate=hover_tpl,
                    showlegend=False,
                )
            )

        for c, col in cond_cols.items():
            fig.add_trace(go.Scatter(x=[None], y=[None],
                                     mode="markers",
                                     marker=dict(color=col, size=10),
                                     name=c))

        grp_vals = ["all"] + list(grp_syms.keys())
        tp_vals = ["all"] + list(tp_sizes.keys())

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
            title=f"PCA — PC1 {props.iloc[0,0]*100:.1f}% • PC2 {props.iloc[0,1]*100:.1f}%",
            xaxis_title="PC1", yaxis_title="PC2",
            margin=dict(l=60, r=200, t=100, b=50),
        )
        return fig

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
    def _volcano(path_full_csv: str, id_col: str, padj_thr: float, fc_thr: float, extra_hover: list[str]) -> PlotlyResource:
        """
        Volcano interactif basé sur le CSV FULL (non filtré),
        coloration par seuils (padj si dispo sinon pvalue).
        """
        d = pd.read_csv(path_full_csv)
        pcol = "padj" if "padj" in d.columns else ("pvalue" if "pvalue" in d.columns else None)
        if pcol is None or "log2FoldChange" not in d.columns:
            pr = PlotlyResource(go.Figure())
            pr.name = f"Volcano {Path(path_full_csv).stem}"
            return pr

        eps = d.loc[d[pcol] > 0, pcol].min()
        eps = eps if (eps is not None and not np.isnan(eps)) else 1e-300

        d["p_plot"] = -np.log10(d[pcol].replace(0, eps))
        d["signif"] = (d[pcol] < float(padj_thr)) & (d["log2FoldChange"].abs() > float(fc_thr))

        hover = []
        for c in [id_col, *extra_hover, pcol, "log2FoldChange", "baseMean"]:
            if c and c in d.columns and c not in hover:
                hover.append(c)

        fig = px.scatter(
            d, x="log2FoldChange", y="p_plot", color="signif",
            hover_data=hover,
            color_discrete_map={True: "#D62728", False: "#A6ACBC"},
            labels={
                "log2FoldChange": "Log₂ FC",
                "p_plot": f"-log₁₀({pcol})",
                "signif": f"{pcol}<{padj_thr} & |log2FC|>{fc_thr}",
            },
            title=Path(path_full_csv).stem.replace("DE_FULL_", ""),
        )
        fig.add_hline(y=-np.log10(float(padj_thr)), line_dash="dot", line_color="grey")
        fig.add_vline(x=float(fc_thr),  line_dash="dot", line_color="grey")
        fig.add_vline(x=-float(fc_thr), line_dash="dot", line_color="grey")

        pr = PlotlyResource(fig)
        pr.name = f"Volcano {Path(path_full_csv).stem}"
        return pr

    def run(self, p: ConfigParams, ins: TaskInputs) -> TaskOutputs:
        counts = ins["count_table_file"].path
        meta   = ins["metadata_file"].path
        gtf    = ins["gtf_file"].path

        samp   = p["samples_colname"]
        feat   = p["genes_colname"]
        cond   = p["condition_column"]
        ctrl   = p["control_condition"]
        tp     = p["timepoint_column"] or "None"
        grp    = p["group_column"]     or "None"
        extras = p["extra_covariates"] or ""
        pthr   = float(p["padj_value"])
        fcthr  = float(p["log2FoldChange_value"])

        annot_cols_str = p.get("annotation_columns", "") or ""
        annot_cols = _split_csv(annot_cols_str)

        shell: ShellProxy = Pydesq2ShellProxyHelper.create_proxy(self.message_dispatcher)

        cmd = " ".join([
            "python3",
            shlex.quote(self.python_file_path),
            shlex.quote(counts),
            shlex.quote(meta),
            shlex.quote(samp),
            shlex.quote(feat),
            shlex.quote(cond),
            shlex.quote(ctrl),
            shlex.quote(tp),
            shlex.quote(grp),
            shlex.quote(extras),
            shlex.quote(str(pthr)),
            shlex.quote(str(fcthr)),
            shlex.quote(gtf),
            shlex.quote(annot_cols_str),
        ])

        if shell.run(cmd, shell_mode=True) != 0:
            raise RuntimeError("DESeq2 analysis failed – check logs")

        work = shell.working_dir

        tables_full = ResourceSet()
        tables_filt = ResourceSet()
        pca = ResourceSet()
        heat = ResourceSet()
        volc = ResourceSet()

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

            # FULL results
            if fn.startswith("DE_FULL_") and fn.endswith(".csv"):
                tables_full.add_resource(
                    TableImporter.call(
                        File(fp),
                        {"delimiter": ",", "header": 0, "file_format": "csv", "index_column": -1},
                    ),
                    fn,
                )
                volc.add_resource(
                    self._volcano(fp, feat, pthr, fcthr, annot_cols),
                    f"Volcano {Path(fp).stem}",
                )

            # FILTERED results
            elif fn.startswith("DE_FILT_") and fn.endswith(".csv"):
                tables_filt.add_resource(
                    TableImporter.call(
                        File(fp),
                        {"delimiter": ",", "header": 0, "file_format": "csv", "index_column": -1},
                    ),
                    fn,
                )

            elif fn.startswith("Summary_") and fn.endswith(".csv"):
                # summaries = filtrés (par construction dans le script)
                tables_filt.add_resource(
                    TableImporter.call(
                        File(fp),
                        {"delimiter": ",", "header": 0, "file_format": "csv"},
                    ),
                    fn,
                )

            elif fn.startswith("Heatmap_") and fn.endswith(".csv"):
                heat.add_resource(self._heatmap(fp), f"Heatmap {Path(fp).stem}")

        return {
            "tables_full":     tables_full,
            "tables_filtered": tables_filt,
            "pca_plots":       pca,
            "heatmap_plots":   heat,
            "volcano_plots":   volc,
        }