from __future__ import annotations
import csv, sys, warnings, re
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds  import DeseqStats


# ───────────────────────── utilities ──────────────────────────────
def _parse_float(x: str | float | None) -> Optional[float]:
    return None if x is None or str(x).lower() in {"", "none", "na"} else float(x)

split_comma = lambda s: [z.strip() for z in s.split(",") if z.strip()]

def _check_samples(count_idx: pd.Index, meta_idx: pd.Index) -> None:
    miss, extra = set(meta_idx) - set(count_idx), set(count_idx) - set(meta_idx)
    if miss:
        raise ValueError(f"Missing samples in counts: {sorted(miss)[:5]} …")
    if extra:
        print(f"⚠  Ignored samples (not in metadata): {sorted(extra)[:5]} …", file=sys.stderr)

def _load_gene_names(gtf: str | None) -> Dict[str, str]:
    if not gtf:
        return {}
    mapping: Dict[str, str] = {}
    with Path(gtf).open() as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if (not row) or row[0].startswith("#") or row[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\w+)\s+"([^"]+)"', row[8]))
            gid, gname = attrs.get("gene_id"), attrs.get("gene_name")
            if gid:
                mapping[gid.split(".")[0]] = gname or gid
    return mapping


# ─────────────────── volcano & heat-map functions ─────────────────
def _volcano(d: pd.DataFrame, label: str,
             pthr: Optional[float], fcthr: Optional[float]) -> None:
    eps = d.loc[d.pvalue > 0, "pvalue"].min() or 1e-300
    d = d.assign(
        neglog10p=-np.log10(d.pvalue.replace(0, eps)),
        signif=((d.pvalue < (pthr or 1)) &
                (d.log2FoldChange.abs() > (fcthr or 0)))
    )
    plt.figure(figsize=(5, 4))
    plt.scatter(d.log2FoldChange[~d.signif], d.neglog10p[~d.signif],
                s=8, alpha=.45, color="#A6ACBC")
    plt.scatter(d.log2FoldChange[d.signif], d.neglog10p[d.signif],
                s=8, alpha=.75, color="#D62728")
    if pthr:
        plt.axhline(-np.log10(pthr), ls="--", c="grey")
    if fcthr:
        plt.axvline( fcthr, ls="--", c="grey")
        plt.axvline(-fcthr, ls="--", c="grey")
    plt.xlabel("log₂ FC"); plt.ylabel("−log₁₀ p")
    plt.title(f"Volcano {label}")
    plt.tight_layout()
    plt.savefig(f"Volcano_{label}.png", dpi=300)
    plt.close()

def _heatmap(ad: sc.AnnData, genes: List[str], label: str) -> None:
    mat = (pd.DataFrame(ad[:, genes].X.T, index=genes, columns=ad.obs_names)
             .replace([np.inf, -np.inf], np.nan).dropna(how="any"))
    if mat.shape[0] < 3 or mat.shape[1] < 2:
        return
    mat.to_csv(f"Heatmap_{label}.csv")
    sns.clustermap(mat, z_score=0, cmap="RdYlBu_r",
                   cbar_kws={"label": "z-score (VST)"}) \
       .fig.savefig(f"Heatmap_{label}.png", dpi=300)
    plt.close()


# ─────────────────────────── MultiDE class ────────────────────────
class MultiDE:
    def __init__(
        self,
        counts_fp: str | Path,
        meta_fp:   str | Path,
        col_sample: str, col_gene: str,
        col_cond: str, ctrl_level: str,
        col_time: str | None, col_group: str | None,
        extras: list[str],
        pthr: Optional[float], fcthr: Optional[float],
        gtf_path: str | None = None,
    ):
        self.counts_fp, self.meta_fp = Path(counts_fp), Path(meta_fp)
        self.c = dict(Sample=col_sample, Gene=col_gene,
                      Cond=col_cond,  Time=col_time,  Group=col_group)
        self.ctrl       = ctrl_level
        self.extras     = extras
        self.pthr       = pthr
        self.fcthr      = fcthr
        self.gene_names = _load_gene_names(gtf_path)
        self.dds: DeseqDataSet | None = None
        self.sum_time:  list[tuple] = []
        self.sum_group: list[tuple] = []

    # ---------- 1 · data load & normalisation ----------------------
    def prepare(self):
        raw = pd.read_csv(self.counts_fp, sep=None, engine="python")
        if self.c["Gene"] not in raw.columns:
            raise ValueError(f"Column '{self.c['Gene']}' missing in count table")

        # keep numeric columns only
        count_mtx = (
            raw.drop(columns=[c for c in raw.columns
                              if c != self.c["Gene"]
                              and not pd.api.types.is_numeric_dtype(raw[c])])
               .set_index(self.c["Gene"])
               .apply(pd.to_numeric, errors="coerce")
               .fillna(0).astype(int)
        )
        count_mtx = count_mtx[count_mtx.sum(axis=1) > 0].T  # samples in rows

        meta = pd.read_csv(self.meta_fp, sep="\t").set_index(self.c["Sample"])
        if self.c["Cond"] not in meta.columns:
            raise ValueError(f"condition_column '{self.c['Cond']}' not in metadata")

        _check_samples(count_mtx.index, meta.index)
        count_mtx = count_mtx.loc[meta.index]
        meta[self.c["Cond"]] = meta[self.c["Cond"]].astype("category")

        # categorical columns
        for col in [self.c["Time"], self.c["Group"], *self.extras]:
            if col and col in meta.columns:
                meta[col] = meta[col].astype("category")

        # design (main effects only)
        terms = [c for c in [self.c["Time"], self.c["Group"], *self.extras]
                 if c and c in meta.columns] + [self.c["Cond"]]
        design = "~ " + " + ".join(terms)

        self.dds = DeseqDataSet(counts=count_mtx, metadata=meta, design=design)
        self.dds.deseq2()

    # ---------- 2 · single contrast --------------------------------
    def _contrast(self, tr: str, ref: str, label: str,
                  dds: DeseqDataSet | None = None,
                  ctx: dict[str, str] | None = None) -> pd.DataFrame:
        ds = dds or self.dds
        st = DeseqStats(ds, contrast=(self.c["Cond"], tr, ref))
        st.summary(save=False)
        df = st.get_results_df() if hasattr(st, "get_results_df") else st.results_df
        if self.pthr  is not None: df = df[df.pvalue < self.pthr]
        if self.fcthr is not None: df = df[df.log2FoldChange.abs() > self.fcthr]
        df = df.copy()
        df["gene_name"] = [self.gene_names.get(i.split(".")[0], "") for i in df.index]
        if ctx:
            for k, v in ctx.items():
                df.insert(0, k, v)
        df.sort_values("log2FoldChange", ascending=False).to_csv(f"DE_{label}.csv")
        _volcano(df, label, self.pthr, self.fcthr)
        return df

    # ---------- 3 · pooled “ALL vs CTRL” ---------------------------
    def _dds_pooled(self, mask: pd.Series) -> DeseqDataSet:
        cts  = pd.DataFrame(self.dds.X[mask.values, :],
                            index=self.dds.obs_names[mask],
                            columns=self.dds.var_names)
        meta = self.dds.obs.loc[mask].copy()
        meta[self.c["Cond"]] = pd.Categorical(
            [self.ctrl if x == self.ctrl else "ALL"
             for x in meta[self.c["Cond"]].astype(str)],
            categories=[self.ctrl, "ALL"])
        dds_p = DeseqDataSet(counts=cts, metadata=meta,
                             design=f"~ {self.c['Cond']}")
        dds_p.deseq2()
        return dds_p

    # ---------- 4 · VST / PCA / heat-maps --------------------------
    def _vst_pca_heat(self, dfs: list[pd.DataFrame], labels: list[str]) -> None:
        self.dds.vst_fit(use_design=False)
        ad = sc.AnnData(self.dds.vst_transform(),
                        obs=self.dds.obs.copy(),
                        var=pd.DataFrame(index=self.dds.var_names))
        sc.pp.pca(ad, n_comps=min(30, ad.n_obs - 1, ad.n_vars - 1))
        (pd.DataFrame(ad.obsm["X_pca"][:, :2],
                      index=ad.obs_names, columns=["PC1", "PC2"])
           .join(ad.obs).to_csv("pca_metadata.csv"))
        v1, v2 = ad.uns["pca"]["variance_ratio"][:2]
        pd.DataFrame({"PC1 Proportion": [v1], "PC2 Proportion": [v2]}) \
          .to_csv("pca_proportions.csv", index=False)

        for d, l in zip(dfs, labels):
            if not d.empty:
                _heatmap(ad, d.nlargest(50, "log2FoldChange").index.tolist(), l)

    # ---------- 5 · per-Time / per-Group contrasts -----------------
    def _stratified(self, col: str | None, store: list[tuple],
                    tag: str, pooled=False):
        if not col or col not in self.dds.obs.columns:
            return
        for lvl in self.dds.obs[col].cat.categories:
            mask = self.dds.obs[col] == lvl
            if mask.sum() < 2:
                continue
            extras = [e for e in self.extras if e != col and e in self.dds.obs.columns]
            design = "~ " + " + ".join([*extras, self.c["Cond"]])
            sub_dds = DeseqDataSet(
                counts=pd.DataFrame(self.dds.X[mask.values, :],
                                    index=self.dds.obs_names[mask],
                                    columns=self.dds.var_names),
                metadata=self.dds.obs.loc[mask], design=design,
            )
            try:
                sub_dds.deseq2()
            except Exception as e:
                warnings.warn(f"{col} {lvl}: {e}")
                continue

            for tr in sub_dds.obs[self.c["Cond"]].cat.categories:
                if tr == self.ctrl:
                    continue
                lbl = f"{tr}_vs_{self.ctrl}_{tag}{lvl}"
                df  = self._contrast(tr, self.ctrl, lbl, sub_dds, {col: lvl})
                store.append((lvl, tr, len(df),
                              (df.log2FoldChange > 0).sum(),
                              (df.log2FoldChange < 0).sum()))

            if pooled:
                try:
                    pool_dds = self._dds_pooled(mask)
                    lbl = f"ALL_vs_{self.ctrl}_{tag}{lvl}"
                    df  = self._contrast("ALL", self.ctrl, lbl, pool_dds, {col: lvl})
                    store.append((lvl, "ALL", len(df),
                                  (df.log2FoldChange > 0).sum(),
                                  (df.log2FoldChange < 0).sum()))
                except Exception as e:
                    warnings.warn(f"{col} {lvl} pooled: {e}")

    # ---------- 6 · write summary ----------------------------------
    @staticmethod
    def _write(rows: Sequence[tuple], csv_path: str, level_name: str):
        if rows:
            pd.DataFrame(rows, columns=[level_name, "Treatment",
                                        "DEGs_total", "DEGs_up", "DEGs_down"]) \
              .to_csv(csv_path, index=False)

    # ---------- 7 · orchestrate run --------------------------------
    def run(self):
        self.prepare()
        dfs, labels = [], []

        # 7.1 each treatment vs CTRL
        for tr in self.dds.obs[self.c["Cond"]].cat.categories:
            if tr == self.ctrl:
                continue
            dfs.append(self._contrast(tr, self.ctrl, f"{tr}_vs_{self.ctrl}"))
            labels.append(f"{tr}_vs_{self.ctrl}")

        # 7.2 pooled ALL vs CTRL
        if len(self.dds.obs[self.c["Cond"]].cat.categories) > 1:
            dds_pool = self._dds_pooled(self.dds.obs[self.c["Cond"]].notna())
            dfs.append(self._contrast("ALL", self.ctrl, f"ALL_vs_{self.ctrl}",
                                      dds_pool))
            labels.append(f"ALL_vs_{self.ctrl}")

        # 7.3 stratified contrasts
        self._stratified(self.c["Time"],  self.sum_time,  'T', True)
        self._stratified(self.c["Group"], self.sum_group, 'G', True)
        if self.sum_time:
            self._write(self.sum_time,  "Summary_Timepoint.csv", self.c["Time"])
        if self.sum_group:
            self._write(self.sum_group, "Summary_Group.csv",     self.c["Group"])

        # 7.4 VST / PCA / heat-maps
        self._vst_pca_heat(dfs, labels)


# ───────────────────────────── CLI ────────────────────────────────
if __name__ == "__main__":
    (
        _, counts, meta,
        col_sample, col_gene,
        col_cond,  ctrl,
        col_time,  col_group,
        extras,
        p_thr, fc_thr, *rest
    ) = sys.argv

    MultiDE(
        counts, meta,
        col_sample, col_gene,
        col_cond, ctrl,
        None if col_time  == "None" else col_time,
        None if col_group == "None" else col_group,
        split_comma(extras),
        _parse_float(p_thr), _parse_float(fc_thr),
        gtf_path=rest[0] if rest else None,
    ).run()
