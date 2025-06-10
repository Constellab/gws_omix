import os, pandas as pd, numpy as np, plotly.express as px
import plotly.graph_objects as go, scipy.cluster.hierarchy as sch

from gws_core import (
    ConfigParams, ConfigSpecs, File,
    FloatParam, InputSpec, InputSpecs,
    OutputSpec, OutputSpecs,
    PlotlyResource, ResourceSet,
    ShellProxy, StrParam, TableImporter,
    Task, TaskInputs, TaskOutputs,
    task_decorator
)
from gws_omix.base_env.pydesq2_env_task import Pydesq2ShellProxyHelper

@task_decorator(
    "pyDESeq2MultiContrast",
    human_name="pyDESeq2 multi-contrast",
    short_description="Compute differential analysis using pyDESeq2 python package: CTRL vs each treatment + pooled all"
)
class Pydesq2Multi(Task):
    """
- MultiDE – automated multi-contrast bulk RNA-seq DEA
    This Python script loops over every treatment level in the metadata file and runs **one differential-expression test per treatment vs control**, then (if ≥ 2 treatments are present) a **pooled “ALL vs CTRL”** comparison.
    The design matrix is chosen automatically:
    • `~ Batch + Condition` when a **`Batch`** column exists, otherwise `~ Condition`.<br>

- **Normalisation workflow:**<br>
    • Before differential expression analysis: Median-of-ratios normalization (dds.deseq2()) is used to correct for differences in sequencing depth across samples.
    • After differential expression analysis: VST (Variance Stabilizing Transformation) is applied for exploratory analyses such as PCA, clustering, and heatmaps.<br>

- **Statistical test:**<br>
     each contrast is analysed with the **Wald test** implemented in PyDESeq2’s `DeseqStats`.<br>

- **The Wald test is used to compare two conditions.**<br>
    The null hypothesis of the Wald test states that for each gene, there is no differential expression between two sample groups (e.g., treated vs. control).<br>
    If the p-value is small (e.g., p < 0.05), the null hypothesis is rejected, suggesting there is only a 5% chance that the observed difference occurred by random chance.<br>
    However, when testing many genes, a number of non-differentially expressed genes may still appear significant due to random chance (false positives).<br>

- **Visual outputs:**<br>
    • `pydesq2_results_table.csv` — DE results sorted by `log2FoldChange` (raw Wald values).<br>
    • **Volcano plot** `Volcano_<contrast>` (red = sig, grey = non-sig).<br>
    • **Heat-map** `Heatmap_<contrast>` for the top-50 DE genes (VST counts).<br>
    • **Global PCA**: `.<br>

- **Example of metadata file:**<br>
    Sample	forward-absolute-filepath	reverse-absolute-filepath	Condition	Replicate	Batch
    SRR13978645	SRR13978645_1.fastq.gz	SRR13978645_2.fastq.gz	CTRL	R1	Batch1
    SRR13978644	SRR13978644_1.fastq.gz	SRR13978644_2.fastq.gz	CTRL	R2	Batch1
    SRR13978643	SRR13978643_1.fastq.gz	SRR13978643_2.fastq.gz	CTRL	R3	Batch2
    SRR13978642	SRR13978642_1.fastq.gz	SRR13978642_2.fastq.gz	SPRC1	R1	Batch1
    SRR13978641	SRR13978641_1.fastq.gz	SRR13978641_2.fastq.gz	SPRC2	R1	Batch1
    SRR13978640	SRR13978640_1.fastq.gz	SRR13978640_2.fastq.gz	SPRC2	R2	Batch2
    """
    input_specs = InputSpecs({
        'count_table_file': InputSpec(File, human_name="counts CSV"),
        'metadata_file':    InputSpec(File, human_name="metadata TSV"),
    })
    output_specs = OutputSpecs({
        "tables":        OutputSpec(ResourceSet, human_name="DE tables"),
        "pca_plots":     OutputSpec(ResourceSet, human_name="Interactive PCA"),
        "heatmap_plots": OutputSpec(ResourceSet, human_name="Interactive heatmaps"),
        "volcano_plots": OutputSpec(ResourceSet, human_name="Interactive volcano plots"),
    })
    config_specs = ConfigSpecs({
        "genes_colname":        StrParam(short_description="gene ID column"),
        "control_condition":    StrParam(short_description="CTRL level"),
        "pvalue_value":         FloatParam(default_value=0.05, min_value=0.0),
        "log2FoldChange_value": FloatParam(default_value=0.5),
    })

    python_file_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "_pydeseq2_multicontrast.py"
    )

    @staticmethod
    def plot_pca(meta_csv, prop_csv):
        df = pd.read_csv(meta_csv); pr = pd.read_csv(prop_csv)
        fig = px.scatter(
            df, x="PC1", y="PC2", color="Condition", text="Sample",
            labels={
                "PC1": f"PC1 ({pr.iloc[0,0]*100:.1f}%)",
                "PC2": f"PC2 ({pr.iloc[0,1]*100:.1f}%)",
            }, title="PCA"
        )
        return PlotlyResource(fig)

    @staticmethod
    def plot_heatmap(path: str) -> PlotlyResource:
        df = pd.read_csv(path, index_col=0)
        rl = sch.linkage(df.values,   method="average", metric="euclidean")
        cl = sch.linkage(df.values.T, method="average", metric="euclidean")
        df = df.iloc[sch.dendrogram(rl,no_plot=True)["leaves"],
                    sch.dendrogram(cl,no_plot=True)["leaves"]]
        fig = go.Figure(go.Heatmap(
            z=df.values,
            x=df.columns,
            y=df.index,
            colorscale="RdYlBu",
            showscale=False
        ))
        base = os.path.splitext(os.path.basename(path))[0]
        fig.update_layout(width=900, height=1200, title=base)
        res = PlotlyResource(fig)
        res.name = f"Interactive Heatmap {base}"
        return res

    @staticmethod
    def plot_volcano(path: str,gene_col: str,*,padj_thr: float,fc_thr: float) -> PlotlyResource:
        # ───────────────────────────────────────────────────────────────
        # 1) Load & prepare data
        # ───────────────────────────────────────────────────────────────
        df = pd.read_csv(path)

        # Prefer ‘padj’ if available
        p_col = "padj" if "padj" in df.columns else "pvalue"

        # Prevent –log10(0) = inf
        eps = df.loc[df[p_col] > 0, p_col].min() or 1e-300
        df["padj_plot"] = -np.log10(df[p_col].replace(0, eps))

        # Decide significance
        df["signif"] = (df[p_col] < padj_thr) & (df["log2FoldChange"].abs() > fc_thr)

        sig_df    = df[df["signif"]]      # True
        nonsig_df = df[~df["signif"]]     # False

        # ───────────────────────────────────────────────────────────────
        # 2) Define friendly labels
        # ───────────────────────────────────────────────────────────────
        labels = {
            "log2FoldChange": "Log₂ Fold Change",
            "padj_plot":      "-log₁₀(adjusted p-value)",
            "signif":         f"{p_col}<{padj_thr} & |log2FC|>{fc_thr}",
        }

        # ───────────────────────────────────────────────────────────────
        # 3) Build figure
        # ───────────────────────────────────────────────────────────────
        fig = go.Figure()

        # False first → True stays on top
        fig.add_trace(go.Scatter(
            x=nonsig_df["log2FoldChange"],
            y=nonsig_df["padj_plot"],
            mode="markers",
            name="False",
            marker=dict(
                color="#A6ACBC",
                size=6,
                opacity=0.55,
                symbol="circle"
            ),
            text=nonsig_df[gene_col],
            hovertemplate=(
                "%{text}"
                "<br>log2FC=%{x:.2f}"
                "<br>-log10(p)=%{y:.2f}"
                "<extra></extra>"
            )
        ))

        fig.add_trace(go.Scatter(
            x=sig_df["log2FoldChange"],
            y=sig_df["padj_plot"],
            mode="markers",
            name="True",
            marker=dict(
                color="#D62728",
                size=8,
                symbol="circle"
            ),
            text=sig_df[gene_col],
            hovertemplate=(
                "%{text}"
                "<br>log2FC=%{x:.2f}"
                "<br>-log10(p)=%{y:.2f}"
                "<extra></extra>"
            )
        ))

        # Threshold lines
        fig.add_hline(
            y=-np.log10(padj_thr),
            line_dash="dot",
            line_color="grey",
            annotation_text=f"{p_col} = {padj_thr}",
            annotation_position="top left"
        )
        fig.add_vline(x= fc_thr, line_dash="dot", line_color="grey")
        fig.add_vline(x=-fc_thr, line_dash="dot", line_color="grey")

        # Layout with the new labels
        fig.update_layout(
            title=os.path.splitext(os.path.basename(path))[0],
            xaxis_title=labels["log2FoldChange"],
            yaxis_title=labels["padj_plot"],
            legend_title_text=labels["signif"],
            width=800,
            height=600,
            margin=dict(l=70, r=20, t=60, b=70)
        )

        res = PlotlyResource(fig)
        res.name = f"Interactive Volcano {os.path.splitext(os.path.basename(path))[0]}"
        return res


    def run(self, params: ConfigParams, inputs: TaskInputs) -> TaskOutputs:
        counts = inputs["count_table_file"].path
        meta   = inputs["metadata_file"].path
        gcol   = params["genes_colname"]
        ctrl   = params["control_condition"]
        p_thr  = params["pvalue_value"]
        fc_thr = params["log2FoldChange_value"]

        shell: ShellProxy = Pydesq2ShellProxyHelper.create_proxy(self.message_dispatcher)
        cmd = f"python3 {self.python_file_path} {counts} {meta} {gcol} {ctrl} {p_thr} {fc_thr}"
        if shell.run(cmd, shell_mode=True) != 0:
            raise RuntimeError("DESeq2 analysis failed")

        work = shell.working_dir
        tables = ResourceSet(); pca_set = ResourceSet()
        heat_set = ResourceSet(); volc_set = ResourceSet()

        pca_set.add_resource(
            self.plot_pca(f"{work}/pca_metadata.csv", f"{work}/pca_proportions.csv"),
            "PCA"
        )

        for fn in os.listdir(work):
            fp = os.path.join(work, fn)
            if fn.startswith("DE_") and fn.endswith(".csv"):
                tables.add_resource(
                    TableImporter.call(File(fp),
                                       {"delimiter":",","header":0,"file_format":"csv","index_column":0}),
                    fn
                )
                volc_set.add_resource(
                    self.plot_volcano(fp, gcol, padj_thr=p_thr, fc_thr=fc_thr),
                    fn.replace("DE_","Volcano_").replace(".csv","")
                )
            elif fn.startswith("Heatmap_") and fn.endswith(".csv"):
                heat_set.add_resource(
                    self.plot_heatmap(fp),
                    fn.replace(".csv","")
                )

        return {
            "tables":        tables,
            "pca_plots":     pca_set,
            "heatmap_plots": heat_set,
            "volcano_plots": volc_set,
        }
