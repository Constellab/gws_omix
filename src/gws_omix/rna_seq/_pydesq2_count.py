import sys
from pathlib import Path
from typing import Optional, Union
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def _parse_float(x: Union[str, float, None]) -> Optional[float]:
    if x is None or (isinstance(x, str) and x.lower() == "none"):
        return None
    return float(x)


def _check_samples(count_idx: pd.Index, meta_idx: pd.Index) -> None:
    missing = set(meta_idx) - set(count_idx)
    extra   = set(count_idx) - set(meta_idx)
    if missing:
        raise ValueError(f"Samples missing in counts: {sorted(missing)}")
    if extra:
        print(f"Extra samples in counts ignored: {sorted(extra)}", file=sys.stderr)


class DeseqAnalyzer:
    def __init__(
        self,
        counts_fp: str | Path,
        meta_fp:   str | Path,
        gene_col:  str,
        ctrl:      str,
        treat:     str,
        pval_thr:  Optional[float] = None,
        fc_thr:    Optional[float] = None,
    ) -> None:
        self.counts_fp = Path(counts_fp)
        self.meta_fp   = Path(meta_fp)
        self.gene_col  = gene_col
        self.ctrl, self.treat = ctrl, treat
        self.pval_thr, self.fc_thr = pval_thr, fc_thr

        self.dds  : DeseqDataSet | None = None
        self.vst  : sc.AnnData   | None = None
        self.sigs : pd.DataFrame | None = None

    # 1 ─ load
    def load(self) -> None:
        counts = (
            pd.read_csv(self.counts_fp, sep=None, engine="python")
              .set_index(self.gene_col)
              .apply(pd.to_numeric, errors="coerce")
              .fillna(0).astype(int)
        )
        counts = counts[counts.sum(axis=1) > 0].T

        meta = pd.read_csv(self.meta_fp, sep="\t", comment="#").set_index("Sample")
        if "Condition" not in meta.columns:
            raise ValueError("metadata.tsv needs a 'Condition' column")
        meta["Condition"] = meta["Condition"].astype("category")

        _check_samples(counts.index, meta.index)
        counts = counts.loc[meta.index]

        self.dds = DeseqDataSet(counts=counts, metadata=meta,
                                design_factors=["Condition"])

    # 2 ─ DE statistics
    def de_stats(self) -> None:
        self.dds.deseq2()
        stats = DeseqStats(self.dds,
                        contrast=("Condition", self.treat, self.ctrl))

        # shrinkage
        try:
            stats.lfc_shrink()
        except TypeError:
            coef = f"Condition_{self.treat}_vs_{self.ctrl}"
            try:
                stats.lfc_shrink(coeff=coef)
            except Exception:
                pass
        except Exception:
            pass

        stats.summary(save=False)
        if hasattr(stats, "get_results_df"):
            res = stats.get_results_df()
        elif hasattr(stats, "results_df"):
            attr = stats.results_df
            res  = attr() if callable(attr) else attr
        else:
            raise RuntimeError("PyDESeq2 table not found (results_df / get_results_df)")

        # filtres utilisateur
        if self.pval_thr is not None:
            res = res[res.pvalue < self.pval_thr]
        if self.fc_thr is not None:
            res = res[res.log2FoldChange.abs() > self.fc_thr]

        self.sigs = res.sort_values("log2FoldChange", ascending=False)
        self.sigs.to_csv("pydesq2_results_table.csv")



    # 3 ─ VST (obligatoire)
    def make_vst(self) -> None:
        # nouvelle API (≥ 0.5)
        if hasattr(self.dds, "vst_fit") and hasattr(self.dds, "vst_transform"):
            self.dds.vst_fit(use_design=False)
            vst_counts = self.dds.vst_transform()
            ad = sc.AnnData(vst_counts,
                            obs=self.dds.obs.copy(),
                            var=pd.DataFrame(index=self.dds.var_names))
            self.vst = ad
            return

        # ancienne API (0.3–0.4)
        if hasattr(self.dds, "get_vst"):
            self.vst = self.dds.get_vst()
            return

        raise RuntimeError(
        )

    # 4 ─ PCA
    def pca(self) -> None:
        sc.pp.pca(self.vst, n_comps=min(30, self.vst.n_obs - 1, self.vst.n_vars - 1))
        (
            pd.DataFrame(self.vst.obsm["X_pca"][:, :2],
                         index=self.vst.obs_names, columns=["PC1", "PC2"])
            .join(self.vst.obs)
            .to_csv("pca_metadata.csv")
        )
        pc1, pc2 = self.vst.uns["pca"]["variance_ratio"][:2]
        pd.DataFrame({"PC1 Proportion":[pc1], "PC2 Proportion":[pc2]}
                     ).to_csv("pca_proportions.csv", index=False)

    # 5 ─ heatmap
    def heatmap(self) -> None:
        if self.sigs is None or self.sigs.empty:
            print("No DE genes – heatmap skipped", file=sys.stderr); return
        ad  = self.vst[:, self.sigs.index]
        mat = pd.DataFrame(ad.X.T, index=ad.var_names, columns=ad.obs_names)
        mat.to_csv("grapher.csv")
        sns.clustermap(mat, z_score=0, cmap="RdYlBu_r")
        plt.savefig("Heatmap.png", dpi=300); plt.close()

    # master
    def run(self) -> None:
        self.load(); self.de_stats(); self.make_vst(); self.pca(); self.heatmap()


# ─────────── CLI ───────────
if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit("Usage: python _pydesq2_count.py counts.csv metadata.tsv "
                 "gene_col ctrl treat pval_thr|None log2fc_thr|None")
    _, counts, meta, gcol, ctrl, treat, pval, fc = sys.argv
    DeseqAnalyzer(
        counts, meta, gcol, ctrl, treat,
        _parse_float(pval), _parse_float(fc)
    ).run()