import sys, warnings
from pathlib import Path
from typing import List, Optional, Union
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds  import DeseqStats

def _parse_float(x: Union[str, float, None]) -> Optional[float]:
    if x is None or (isinstance(x, str) and x.lower() == "none"):
        return None
    return float(x)

def _check_samples(count_idx: pd.Index, meta_idx: pd.Index) -> None:
    miss = set(meta_idx) - set(count_idx)
    extra = set(count_idx) - set(meta_idx)
    if miss:
        raise ValueError(f"Missing samples in counts: {sorted(miss)}")
    if extra:
        print(f"⚠️ Extra samples ignored: {sorted(extra)}", file=sys.stderr)

def _volcano(df: pd.DataFrame, label: str,
             pthr: Optional[float], fcthr: Optional[float]) -> None:
    # prepare
    eps = df.loc[df.pvalue > 0, "pvalue"].min() or 1e-300
    df = df.assign(neglog10p = -np.log10(df.pvalue.replace(0, eps)))
    # color
    sig = (df.pvalue < pthr) & (df.log2FoldChange.abs() > fcthr)
    plt.figure(figsize=(5,4))
    plt.scatter(df.log2FoldChange[~sig], df.neglog10p[~sig],
                s=8, alpha=0.4, label="non-signif", color="grey")
    plt.scatter(df.log2FoldChange[sig], df.neglog10p[sig],
                s=8, alpha=0.6, label="signif", color="red")
    # thresholds
    if pthr:
        plt.axhline(-np.log10(pthr), ls="--", c="grey")
    if fcthr:
        plt.axvline( fcthr, ls="--", c="grey")
        plt.axvline(-fcthr, ls="--", c="grey")
    plt.xlabel("log₂ fold change")
    plt.ylabel("P-value (−Log₁₀)")
    plt.title(f"Volcano {label}")
    plt.legend(frameon=False, fontsize="small")
    plt.tight_layout()
    plt.savefig(f"Volcano_{label}.png", dpi=300)
    plt.close()

def _heatmap(ad: sc.AnnData, genes: List[str], label: str, figsize=(12,12)) -> None:
    mat = pd.DataFrame(ad[:, genes].X.T, index=genes, columns=ad.obs_names)
    mat.to_csv(f"Heatmap_{label}.csv")
    cg = sns.clustermap(mat, z_score=0, cmap="RdYlBu_r", figsize=figsize)
    cg.ax_heatmap.set_title(f"Heatmap {label}")
    plt.savefig(f"Heatmap_{label}.png", dpi=300)
    plt.close()

class MultiDE:
    def __init__(self, counts_fp, meta_fp, gene_col, ctrl, pthr, fcthr):
        self.counts_fp = Path(counts_fp)
        self.meta_fp   = Path(meta_fp)
        self.gene_col  = gene_col
        self.ctrl      = ctrl
        self.pthr, self.fcthr = pthr, fcthr
        self.dds = None

    def prepare(self):
        counts = (pd.read_csv(self.counts_fp, sep=None, engine="python")
                    .set_index(self.gene_col)
                    .apply(pd.to_numeric, errors="coerce")
                    .fillna(0).astype(int))
        counts = counts[counts.sum(axis=1)>0].T

        meta = pd.read_csv(self.meta_fp, sep="\t", comment="#").set_index("Sample")
        if "Condition" not in meta.columns:
            raise ValueError("metadata must have a 'Condition' column")
        meta["Condition"] = meta["Condition"].astype("category")

        _check_samples(counts.index, meta.index)
        counts = counts.loc[meta.index]

        design = "~ Batch + Condition" if "Batch" in meta.columns else "~ Condition"
        self.dds = DeseqDataSet(counts=counts, metadata=meta, design=design)
        self.dds.deseq2()

    def _contrast(self, grp, ref, label, dds=None):
        dds = dds or self.dds
        stats = DeseqStats(dds, contrast=("Condition", grp, ref))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try: stats.lfc_shrink()
            except: pass
        stats.summary(save=False)
        if hasattr(stats, "get_results_df"):
            df = stats.get_results_df()
        else:
            r = stats.results_df
            df = r() if callable(r) else r
        if self.pthr  is not None:
            df = df[df.pvalue < self.pthr]
        if self.fcthr is not None:
            df = df[df.log2FoldChange.abs() > self.fcthr]
        df.sort_values("log2FoldChange", ascending=False)\
          .to_csv(f"DE_{label}.csv")
        _volcano(df, label, self.pthr, self.fcthr)
        return df

    def pooled_contrast(self):
        meta2 = self.dds.obs.copy()
        collapsed = [x if x==self.ctrl else "ALL" for x in meta2["Condition"].astype(str)]
        meta2["Condition"] = pd.Categorical(collapsed, categories=[self.ctrl, "ALL"])
        dds2 = DeseqDataSet(
            counts   = pd.DataFrame(self.dds.X, index=self.dds.obs_names, columns=self.dds.var_names),
            metadata = meta2,
            design   = "~ Condition"
        )
        try:
            dds2.deseq2()
        except Exception as e:
            print(f"pooled_contrast skipped: {e}", file=sys.stderr)
            return None
        return self._contrast("ALL", self.ctrl, f"ALL_vs_{self.ctrl}", dds2)

    def _vst_pca_heat(self, dfs: List[pd.DataFrame], labels: List[str]):
        self.dds.vst_fit(use_design=False)
        vst = self.dds.vst_transform()
        ad  = sc.AnnData(vst, obs=self.dds.obs.copy(), var=pd.DataFrame(index=self.dds.var_names))

        # PCA
        sc.pp.pca(ad, n_comps=min(30, ad.n_obs-1, ad.n_vars-1))
        pd.DataFrame(ad.obsm["X_pca"][:,:2], index=ad.obs_names, columns=["PC1","PC2"])\
          .join(ad.obs).to_csv("pca_metadata.csv")
        v1, v2 = ad.uns["pca"]["variance_ratio"][:2]
        pd.DataFrame({"PC1 Proportion":[v1],"PC2 Proportion":[v2]})\
          .to_csv("pca_proportions.csv", index=False)

        for df, lbl in zip(dfs, labels):
            if df is None:
                continue
            top = df.nlargest(50, "log2FoldChange").index.tolist()
            if lbl.startswith("ALL_vs_"):
                ad_sub = ad
            else:
                treat = lbl.split("_vs_")[0]
                ad_sub = ad[ad.obs["Condition"].isin([self.ctrl, treat]), :]
            _heatmap(ad_sub, top, lbl)

    def run(self):
        self.prepare()
        levels = [l for l in self.dds.obs["Condition"].cat.categories if l != self.ctrl]

        dfs, lbls = [], []
        for t in levels:
            dfs.append(self._contrast(t, self.ctrl, f"{t}_vs_{self.ctrl}"))
            lbls.append(f"{t}_vs_{self.ctrl}")

        if len(levels) > 1:
            pooled = self.pooled_contrast()
            if pooled is not None:
                dfs.append(pooled)
                lbls.append(f"ALL_vs_{self.ctrl}")

        self._vst_pca_heat(dfs, lbls)

if __name__ == "__main__":
    _, counts, meta, gcol, ctrl, p, fc = sys.argv
    MultiDE(counts, meta, gcol, ctrl, _parse_float(p), _parse_float(fc)).run()
