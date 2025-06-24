from __future__ import annotations
import csv, sys, warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np, pandas as pd, scanpy as sc, seaborn as sns, matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds  import DeseqStats


def _parse_float(x: Union[str, float, None]) -> Optional[float]:
    return None if x is None or (isinstance(x, str) and x.lower() == "none") else float(x)


def _check_samples(count_idx: pd.Index, meta_idx: pd.Index) -> None:
    miss, extra = set(meta_idx) - set(count_idx), set(count_idx) - set(meta_idx)
    if miss:
        raise ValueError(f"Missing samples in counts: {sorted(miss)}")
    if extra:
        print(f"⚠️ Extra samples ignored: {sorted(extra)}", file=sys.stderr)


def _load_gene_names(gtf_path: Optional[Union[str, Path]]) -> Dict[str, str]:
    if gtf_path is None:
        return {}
    gtf_path = Path(gtf_path)
    if not gtf_path.is_file():
        raise FileNotFoundError(gtf_path)

    mapping: Dict[str, str] = {}
    with gtf_path.open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#") or row[2] != "gene":
                continue
            attrs = {
                k.strip(): v.strip('" ')
                for part in row[8].split(";") if part.strip()
                for k, v in [part.strip().split(" ", 1)]
            }
            gid = attrs.get("gene_id")
            gname = attrs.get("gene_name", gid)
            if gid:
                mapping[gid.split(".")[0]] = gname
    return mapping


# ───────────────────── volcano / heatmap ───────────────────
def _volcano(df: pd.DataFrame, label: str,
             pthr: Optional[float], fcthr: Optional[float]) -> None:
    eps = df.loc[df.pvalue > 0, "pvalue"].min() or 1e-300
    df  = df.assign(neglog10p=-np.log10(df.pvalue.replace(0, eps)))
    sig = (df.pvalue < pthr) & (df.log2FoldChange.abs() > fcthr)

    plt.figure(figsize=(5, 4))
    plt.scatter(df.log2FoldChange[~sig], df.neglog10p[~sig],
                s=8, alpha=0.45, color="#A6ACBC", label="False")
    plt.scatter(df.log2FoldChange[sig], df.neglog10p[sig],
                s=8, alpha=0.75, color="#D62728", label="True")

    if pthr:  plt.axhline(-np.log10(pthr), ls="--", c="grey")
    if fcthr:
        plt.axvline( fcthr, ls="--", c="grey")
        plt.axvline(-fcthr, ls="--", c="grey")

    plt.xlabel("log₂ fold change")
    plt.ylabel("P value (−log₁₀)")
    plt.title(f"Volcano {label}")
    plt.legend(frameon=False, fontsize="small")
    plt.tight_layout()
    plt.savefig(f"Volcano_{label}.png", dpi=300)
    plt.close()


def _heatmap(ad: sc.AnnData, genes: List[str], label: str,
             figsize=(12, 12)) -> None:
    """
    Construit une heatmap VST (z-score) des 50 gènes sélectionnés.
    • Ignore les gènes/samples contenant des NaN ou Inf
    • Saute la figure si < 3 gènes ou < 2 échantillons après nettoyage
    """
    mat = pd.DataFrame(ad[:, genes].X.T, index=genes, columns=ad.obs_names)

    # ── nettoyage pour éviter NaN / ±Inf ─────────────────────────
    mat = mat.replace([np.inf, -np.inf], np.nan).dropna(how="any")

    # ── skip si matrice trop petite  ------------------------------
    if mat.shape[0] < 3 or mat.shape[1] < 2:
        return  # pas assez d’information → on ne génère pas la heat-map

    # ── sauvegarde CSV nettoyé (utile pour debug) -----------------
    mat.to_csv(f"Heatmap_{label}.csv")

    cg = sns.clustermap(
        mat, z_score=0, cmap="RdYlBu_r", figsize=figsize,
        cbar_kws={"label": "z-score (VST)"}
    )
    cg.ax_heatmap.set_title(f"Heatmap {label}")
    plt.savefig(f"Heatmap_{label}.png", dpi=300)
    plt.close()


#MultiDE class
class MultiDE:
    def __init__(self, counts_fp, meta_fp, gene_col, ctrl,
                 pthr: Optional[float], fcthr: Optional[float],
                 gtf_path: Optional[str] = None):
        self.counts_fp = Path(counts_fp)
        self.meta_fp   = Path(meta_fp)
        self.gene_col  = gene_col
        self.ctrl      = ctrl
        self.pthr, self.fcthr = pthr, fcthr
        self.dds: DeseqDataSet | None = None
        self.gene_names = _load_gene_names(gtf_path)

        self._summary_time, self._summary_group = [], []

    # ---------- 1 · load & fit --------------------------------------
    def prepare(self) -> None:
        counts = (pd.read_csv(self.counts_fp, sep=None, engine="python")
                    .set_index(self.gene_col)
                    .apply(pd.to_numeric, errors="coerce")
                    .fillna(0).astype(int))
        counts = counts[counts.sum(axis=1) > 0].T

        meta = pd.read_csv(self.meta_fp, sep="\t", comment="#").set_index("Sample")
        if "Condition" not in meta.columns:
            raise ValueError("Metadata must contain a 'Condition' column")
        meta["Condition"] = meta["Condition"].astype("category")

        _check_samples(counts.index, meta.index)
        counts = counts.loc[meta.index]

        extras = [c for c in ["Batch", "Timepoint", "Replicate", "Group"]
                  if c in meta.columns]
        for c in extras:
            meta[c] = meta[c].astype("category")

        self.dds = DeseqDataSet(
            counts=counts,
            metadata=meta,
            design="~ " + " + ".join(extras + ["Condition"]),
        )
        self.dds.deseq2()

    # ---------- 2 · single contrast --------------------------------
    def _contrast(self, group, ref, label,
                  dds: Optional[DeseqDataSet] = None,
                  ctx: Optional[dict] = None) -> pd.DataFrame:
        dds = dds or self.dds
        st  = DeseqStats(dds, contrast=("Condition", group, ref))
        st.summary(save=False)
        df = st.get_results_df() if hasattr(st, "get_results_df") else st.results_df

        if self.pthr is not None:
            df = df[df.pvalue < self.pthr]
        if self.fcthr is not None:
            df = df[df.log2FoldChange.abs() > self.fcthr]

        df = df.copy()
        df["gene_name"] = [self.gene_names.get(i.split(".")[0], "") for i in df.index]

        for col in ("Timepoint", "Group"):
            if ctx and col in ctx:
                df.insert(0, col, ctx[col])

        df.sort_values("log2FoldChange", ascending=False)\
          .to_csv(f"DE_{label}.csv")
        _volcano(df, label, self.pthr, self.fcthr)
        return df

    # ---------- 3 · pooled utilities -------------------------------
    def _pooled_dataset(self, mask: pd.Series) -> DeseqDataSet:
        counts_sub = pd.DataFrame(self.dds.X[mask.values, :],
                                  index=self.dds.obs_names[mask],
                                  columns=self.dds.var_names)
        meta_sub = self.dds.obs.loc[mask].copy()
        meta_sub["Condition"] = pd.Categorical(
            [self.ctrl if x == self.ctrl else "ALL"
             for x in meta_sub["Condition"].astype(str)],
            categories=[self.ctrl, "ALL"])
        return DeseqDataSet(counts=counts_sub, metadata=meta_sub,
                            design="~ Condition")

    def pooled_contrast_global(self):
        try:
            dds2 = self._pooled_dataset(self.dds.obs["Condition"].notna())
            dds2.deseq2()
            return self._contrast("ALL", self.ctrl, f"ALL_vs_{self.ctrl}", dds2)
        except Exception as e:
            print(f"⚠ pooled_contrast (global) skipped: {e}", file=sys.stderr)
            return None

    # ---------- 4 · VST → PCA & heatmaps ---------------------------
    def _vst_pca_heat(self, dfs: List[pd.DataFrame], labels: List[str]) -> None:
        self.dds.vst_fit(use_design=False)
        ad = sc.AnnData(self.dds.vst_transform(),
                        obs=self.dds.obs.copy(),
                        var=pd.DataFrame(index=self.dds.var_names))

        # global PCA
        sc.pp.pca(ad, n_comps=min(30, ad.n_obs - 1, ad.n_vars - 1))
        (pd.DataFrame(ad.obsm["X_pca"][:, :2],
                      index=ad.obs_names, columns=["PC1", "PC2"])
           .join(ad.obs).to_csv("pca_metadata.csv"))
        v1, v2 = ad.uns["pca"]["variance_ratio"][:2]
        pd.DataFrame({"PC1 Proportion": [v1], "PC2 Proportion": [v2]})\
          .to_csv("pca_proportions.csv", index=False)

        # per Timepoint & per Group PCA
        def _subset_pca(mask: pd.Series, tag: str):
            if mask.sum() < 3:
                return
            ad_sub = ad[mask, :]
            sc.pp.pca(ad_sub, n_comps=min(30, ad_sub.n_obs - 1, ad_sub.n_vars - 1))
            (pd.DataFrame(ad_sub.obsm["X_pca"][:, :2],
                          index=ad_sub.obs_names, columns=["PC1", "PC2"])
               .join(ad_sub.obs).to_csv(f"pca_metadata_{tag}.csv"))
            v1, v2 = ad_sub.uns["pca"]["variance_ratio"][:2]
            pd.DataFrame({"PC1 Proportion": [v1], "PC2 Proportion": [v2]})\
              .to_csv(f"pca_proportions_{tag}.csv", index=False)

        if "Timepoint" in ad.obs.columns:
            for tp in ad.obs["Timepoint"].cat.categories:
                _subset_pca(ad.obs["Timepoint"] == tp, f"T{tp}")

        if "Group" in ad.obs.columns:
            for gp in ad.obs["Group"].cat.categories:
                _subset_pca(ad.obs["Group"] == gp, f"G{gp}")

        # heatmaps for every DE result
        for df, lbl in zip(dfs, labels):
            if df is None or df.empty:
                continue
            genes = df.nlargest(50, "log2FoldChange").index.tolist()
            ad_sub = ad

            if "_T" in lbl:
                ad_sub = ad_sub[ad_sub.obs.get("Timepoint", "") == lbl.split("_T")[-1], :]
            if "_G" in lbl:
                ad_sub = ad_sub[ad_sub.obs.get("Group", "") == lbl.split("_G")[-1], :]

            if not lbl.startswith("ALL_vs_"):
                tr = lbl.split("_vs_")[0]
                ad_sub = ad_sub[ad_sub.obs["Condition"].isin([self.ctrl, tr]), :]

            _heatmap(ad_sub, genes, lbl)

    # 5 stratified contrasts (+pooled)
    def _run_stratified(self, column, summary_holder, tag, do_pooled=False):
        if column not in self.dds.obs.columns:
            return

        for lvl in self.dds.obs[column].cat.categories:
            mask = self.dds.obs[column] == lvl
            if mask.sum() < 2:
                continue

            counts_sub = pd.DataFrame(self.dds.X[mask.values, :],
                                      index=self.dds.obs_names[mask],
                                      columns=self.dds.var_names)
            meta_sub = self.dds.obs.loc[mask].copy()

            extras = [c for c in ["Batch", "Timepoint", "Replicate", "Group"]
                      if c in meta_sub.columns and c != column]
            design = ("~ " + " + ".join(extras + ["Condition"])
                      if extras else "~ Condition")

            dds_sub = DeseqDataSet(counts=counts_sub,
                                   metadata=meta_sub,
                                   design=design)
            try:
                dds_sub.deseq2()
            except Exception as e:
                warnings.warn(f"{column} level {lvl} failed: {e}")
                continue

            for tr in meta_sub["Condition"].cat.categories:
                if tr == self.ctrl:
                    continue
                lbl = f"{tr}_vs_{self.ctrl}_{tag}{lvl}"
                ctx = {column: lvl}
                df  = self._contrast(tr, self.ctrl, lbl, dds_sub, ctx)

                summary_holder.append(
                    (lvl, tr, len(df),
                     (df.log2FoldChange > 0).sum(),
                     (df.log2FoldChange < 0).sum())
                )

            # pooled ALL-vs-CTRL inside the level
            if do_pooled:
                try:
                    dds_pool = self._pooled_dataset(mask)
                    dds_pool.deseq2()
                    lbl = f"ALL_vs_{self.ctrl}_{tag}{lvl}"
                    df  = self._contrast("ALL", self.ctrl, lbl, dds_pool,
                                         ctx={column: lvl})
                    summary_holder.append(
                        (lvl, "ALL", len(df),
                         (df.log2FoldChange > 0).sum(),
                         (df.log2FoldChange < 0).sum())
                    )
                except Exception as e:
                    warnings.warn(f"{column} level {lvl}: pooled failed ({e})")

    # 6 · summary writer
    @staticmethod
    def _write_summary(rows: Sequence[tuple], out_csv: str, level: str) -> None:
        if rows:
            pd.DataFrame(rows,
                         columns=[level, "Treatment",
                                  "DEGs_total", "DEGs_up", "DEGs_down"])\
              .to_csv(out_csv, index=False)

    # 7 · orchestrator
    def run(self) -> None:
        self.prepare()

        dfs, labels = [], []
        for tr in self.dds.obs["Condition"].cat.categories:
            if tr == self.ctrl:
                continue
            dfs.append(self._contrast(tr, self.ctrl, f"{tr}_vs_{self.ctrl}"))
            labels.append(f"{tr}_vs_{self.ctrl}")

        pooled = self.pooled_contrast_global()
        if pooled is not None:
            dfs.append(pooled); labels.append(f"ALL_vs_{self.ctrl}")

        self._run_stratified("Timepoint", self._summary_time,  "T", do_pooled=True)
        self._run_stratified("Group",     self._summary_group, "G", do_pooled=True)

        self._write_summary(self._summary_time,  "Summary_Timepoint.csv", "Timepoint")
        self._write_summary(self._summary_group, "Summary_Group.csv",     "Group")

        self._vst_pca_heat(dfs, labels)


if __name__ == "__main__":
    _, counts, meta, gcol, ctrl, p_thr, fc_thr, *rest = sys.argv
    gtf = rest[0] if rest else None

    MultiDE(counts, meta, gcol, ctrl,
            _parse_float(p_thr), _parse_float(fc_thr),
            gtf_path=gtf).run()
