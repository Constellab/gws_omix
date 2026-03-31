# gsea.py
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
import plotly.graph_objects as go

from gws_core import (
    ConfigParams,
    ConfigSpecs,
    File,
    IntParam,
    InputSpec,
    InputSpecs,
    OutputSpec,
    OutputSpecs,
    PlotlyResource,
    ResourceSet,
    ShellProxy,
    StrParam,
    Table,
    TableExporter,
    TableImporter,
    Task,
    TaskInputs,
    TaskOutputs,
    task_decorator,
)

from .gsea_env import GseaShellProxyHelper


def _grid_from_long_gsea(
    long_df: pd.DataFrame,
    title: str,
    genes_per_page: int = 25,
    sets_per_page: int = 20,
) -> PlotlyResource:
    """
    Interactive grid: Gene sets × Leading-edge genes
    - Color: padj
    - Size: |stat|
    """
    if long_df is None or long_df.empty:
        fig = go.Figure()
        fig.update_layout(title=title, width=1400, height=420)
        return PlotlyResource(fig)

    d = long_df.copy()
    d["set_id"] = d["set_id"].astype(str)
    d["gene_label"] = d["gene_label"].astype(str)
    d["padj"] = pd.to_numeric(d.get("padj", np.nan), errors="coerce")
    d["stat"] = pd.to_numeric(d.get("stat", np.nan), errors="coerce")
    d["NES"] = pd.to_numeric(d.get("NES", np.nan), errors="coerce")

    order = (
        d.groupby("set_id")
        .agg(
            best_padj=("padj", "min"),
            best_abs_nes=("NES", lambda x: np.nanmax(np.abs(pd.to_numeric(x, errors="coerce")))),
        )
        .reset_index()
        .sort_values(["best_padj", "best_abs_nes"], ascending=[True, False])
    )
    sets = order["set_id"].tolist()
    y_positions = np.arange(len(sets), dtype=float)
    ymap = dict(zip(sets, y_positions))
    d["_y"] = d["set_id"].map(ymap)

    genes = sorted(d["gene_label"].dropna().unique().tolist())
    pages_g = [genes[i : i + genes_per_page] for i in range(0, len(genes), genes_per_page)] or [genes]

    sp = max(1, int(sets_per_page))
    windows = [(s, min(s + sp, len(sets))) for s in range(0, len(sets), sp)]
    ypad = 0.7

    def y_layout(win):
        s, e = win
        return dict(
            tickmode="array",
            tickvals=y_positions[s:e].tolist(),
            ticktext=sets[s:e],
            autorange=False,
            range=[(e - 0.5) + ypad, (s - 0.5) - ypad],
            showgrid=True,
            gridcolor="#f7f7f7",
            tickfont=dict(size=11),
        )

    p = d["padj"].replace([np.inf, -np.inf], np.nan)
    max95 = np.nanpercentile(p, 95) if np.isfinite(p).any() else 1.0
    d["_p"] = np.clip(p, 0, max95)

    abs_stat = d["stat"].abs().replace([np.inf, -np.inf], np.nan)
    s95 = np.nanpercentile(abs_stat, 95) if np.isfinite(abs_stat).any() else 1.0
    scale = (abs_stat / (s95 if s95 > 0 else 1.0)).fillna(0.0)
    d["_sz"] = (8.0 + 10.0 * np.clip(scale, 0, 1)).astype(float)

    fig = go.Figure()
    for gi, page_genes in enumerate(pages_g, start=1):
        sub = d[d["gene_label"].isin(page_genes)].copy()
        xmap = {g: i for i, g in enumerate(page_genes)}
        sub["_x"] = sub["gene_label"].map(xmap)

        fig.add_trace(
            go.Scatter(
                x=sub["_x"],
                y=sub["_y"],
                mode="markers",
                marker=dict(
                    size=sub["_sz"],
                    color=sub["_p"],
                    colorscale="Viridis",
                    reversescale=True,
                    colorbar=dict(title="p.adjust", thickness=12, len=0.85),
                    line=dict(width=0.3, color="#444"),
                    opacity=0.9,
                ),
                visible=(gi == 1),
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Gene: %{customdata[1]}<br>"
                    "NES: %{customdata[2]}<br>"
                    "p.adjust: %{customdata[3]}<br>"
                    "stat: %{customdata[4]}<extra></extra>"
                ),
                customdata=np.stack(
                    [
                        sub["set_id"].values,
                        sub["gene_label"].values,
                        pd.to_numeric(sub["NES"], errors="coerce").round(3).astype(object).values,
                        pd.to_numeric(sub["padj"], errors="coerce")
                        .map(lambda v: f"{v:.2e}" if np.isfinite(v) else "NA")
                        .values,
                        pd.to_numeric(sub["stat"], errors="coerce").round(3).astype(object).values,
                    ],
                    axis=1,
                ),
                name=f"Genes page {gi}",
            )
        )

    buttons_genes = []
    for gi, page_genes in enumerate(pages_g, start=1):
        vis = [False] * len(pages_g)
        vis[gi - 1] = True
        buttons_genes.append(
            dict(
                label=f"Genes {gi}/{len(pages_g)}",
                method="update",
                args=[
                    {"visible": vis},
                    {
                        "xaxis": {
                            "tickmode": "array",
                            "tickvals": list(range(len(page_genes))),
                            "ticktext": page_genes,
                            "tickangle": 45,
                            "tickfont": {"size": 10},
                        }
                    },
                ],
            )
        )

    buttons_sets = []
    for wi, win in enumerate(windows, start=1):
        buttons_sets.append(dict(label=f"Sets {wi}/{len(windows)}", method="relayout", args=[{"yaxis": y_layout(win)}]))

    fig.update_layout(
        title=title,
        width=1600,
        height=min(8200, 520 + 30 * max(1, len(sets))),
        xaxis=dict(
            title="<b>Leading-edge genes</b>",
            tickmode="array",
            tickvals=list(range(len(pages_g[0]))),
            ticktext=pages_g[0],
            tickangle=45,
            showgrid=True,
            gridcolor="#f0f0f0",
        ),
        yaxis=y_layout(windows[0]),
        updatemenus=[
            dict(type="dropdown", x=0.80, y=1.12, showactive=True, buttons=buttons_sets),
            dict(type="dropdown", x=0.90, y=1.12, showactive=True, buttons=buttons_genes),
        ],
        margin=dict(l=340, r=90, t=120, b=180),
    )
    return PlotlyResource(fig)


@task_decorator(
    "GSEAEnrichmentGMT",
    human_name="Gene set enrichment analysis (GSEA)",
    short_description="Pre-ranked GSEA using a signed ranking statistic (stat column) and a single GMT gene set file."
)
class GSEAEnrichmentGMT(Task):

    input_specs: Final[InputSpecs] = InputSpecs(
        {
            "de_table_file": InputSpec((Table, File), human_name="DE/Full results (Table or CSV file)"),
            "gmt_file": InputSpec(File, human_name="Gene sets (GMT file)"),
        }
    )

    output_specs: Final[OutputSpecs] = OutputSpecs(
        {
            "results": OutputSpec(ResourceSet, human_name="GSEA results (table + plots + grid)"),
            "enrichplots": OutputSpec(ResourceSet, human_name="GSEA enrichment plots (PNG)"),
        }
    )

    config_specs: Final[ConfigSpecs] = ConfigSpecs(
        {
            "genes_colname": StrParam(
                default_value="gene",
                short_description="Column containing gene identifiers. Default: 'gene'.",
            ),
            "stat_colname": StrParam(
                default_value="stat",
                short_description="Column containing the signed ranking statistic. Default: 'stat'.",
            ),
            "min_size": IntParam(default_value=15, min_value=1, short_description="Minimum gene set size after overlap."),
            "max_size": IntParam(default_value=500, min_value=15, short_description="Maximum gene set size."),
            "nperm": IntParam(default_value=1000, min_value=1000, short_description="Number of permutations."),
            "processes": IntParam(default_value=4, min_value=1, short_description="Parallel processes used by gseapy."),
            "n_enrichplots": IntParam(default_value=10, min_value=0, short_description="Number of enrichment plots (PNG) to generate."),
            "grid_genes_per_page": IntParam(default_value=25, min_value=1, short_description="Genes per page in the interactive grid."),
            "grid_sets_per_page": IntParam(default_value=20, min_value=1, short_description="Gene sets per page in the interactive grid."),
        }
    )

    python_file_path: Final[str] = os.path.join(os.path.abspath(os.path.dirname(__file__)), "_gsea.py")

    def run(self, p: ConfigParams, ins: TaskInputs) -> TaskOutputs:
        de_res = ins["de_table_file"]
        if isinstance(de_res, Table):
            de_file: File = TableExporter.call(de_res, {"file_format": "csv"})
            de_csv = de_file.path
        elif isinstance(de_res, File):
            de_csv = de_res.path
        else:
            raise RuntimeError("Unsupported type for 'de_table_file' (expected Table or File).")

        gmt: File = ins["gmt_file"]

        genes_col = str(p.get("genes_colname", "gene") or "gene").strip()
        stat_col = str(p.get("stat_colname", "stat") or "stat").strip()

        shell: ShellProxy = GseaShellProxyHelper.create_proxy(self.message_dispatcher)
        work = Path(shell.working_dir)

        csv_out = work / "csv"
        fig_out = work / "figs"
        csv_out.mkdir(parents=True, exist_ok=True)
        fig_out.mkdir(parents=True, exist_ok=True)

        script_path = Path(os.path.dirname(os.path.realpath(__file__))) / "_gsea.py"

        cmd = [
            "python3",
            str(script_path),
            "--infile", de_csv,
            "--genes_colname", genes_col,
            "--stat_colname", stat_col,
            "--gmt_file", gmt.path,
            "--outdir", str(work),
            "--csv_dir", str(csv_out),
            "--fig_dir", str(fig_out),
            "--min_size", str(int(p["min_size"])),
            "--max_size", str(int(p["max_size"])),
            "--nperm", str(int(p["nperm"])),
            "--processes", str(int(p["processes"])),
            "--n_enrichplots", str(int(p["n_enrichplots"])),
        ]
        cmd_str = " ".join(shlex.quote(str(x)) for x in cmd)

        if shell.run(cmd_str, shell_mode=True) != 0:
            return {"results": ResourceSet(), "enrichplots": ResourceSet()}

        # ── Base ResourceSet ─────────────────────────────────────────
        rs = ResourceSet()
        rs.name = "GSEA results"

        # 1) Results table
        table_fp = csv_out / "GSEA_results.csv"
        if table_fp.exists() and table_fp.stat().st_size > 0:
            rs.add_resource(
                TableImporter.call(File(str(table_fp)), {"delimiter": ",", "header": 0, "file_format": "csv"}),
                "GSEA_results.csv",
            )

        long_fp = csv_out / "GSEA_matrix_long.csv"
        long_df = pd.DataFrame()
        if long_fp.exists() and long_fp.stat().st_size > 0:
            try:
                long_tbl = TableImporter.call(File(str(long_fp)), {"delimiter": ",", "header": 0, "file_format": "csv"})
                rs.add_resource(long_tbl, "GSEA_matrix_long.csv")
            except Exception:
                # Fallback: keep it as a File so it's not lost
                rs.add_resource(File(str(long_fp)), "GSEA_matrix_long.csv")

            # For grid creation (need dataframe)
            try:
                long_df = pd.read_csv(long_fp)
            except Exception:
                long_df = pd.DataFrame()

        # 3) Add all PNGs except enrichment plots into base output
        for fp in sorted(fig_out.glob("GSEA_*.png")):
            if not fp.exists() or fp.stat().st_size == 0:
                continue
            if fp.name.startswith("GSEA_enrichment_"):
                continue
            rs.add_resource(File(str(fp)), fp.stem)

        # 4) Grid stays in base output
        try:
            pr = _grid_from_long_gsea(
                long_df,
                title="Gene sets × leading-edge genes",
                genes_per_page=int(p["grid_genes_per_page"]),
                sets_per_page=int(p["grid_sets_per_page"]),
            )
            rs.add_resource(pr, "GSEA_grid")
        except Exception:
            fig = go.Figure()
            fig.update_layout(title="Gene sets × leading-edge genes", width=1200, height=400)
            rs.add_resource(PlotlyResource(fig), "GSEA_grid")

        # ── Enrichment plots ResourceSet ─────────────────────────────
        rs_enrich = ResourceSet()
        rs_enrich.name = "GSEA enrichment plots"

        for fp in sorted(fig_out.glob("GSEA_enrichment_*.png")):
            if fp.exists() and fp.stat().st_size > 0:
                rs_enrich.add_resource(File(str(fp)), fp.stem)

        return {"results": rs, "enrichplots": rs_enrich}