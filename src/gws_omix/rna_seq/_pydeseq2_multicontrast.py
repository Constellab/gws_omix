from __future__ import annotations

import re
import sys
import warnings
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# ───────────────────────── utilities ──────────────────────────────
def _parse_float(x: str | float | None) -> float | None:
    return None if x is None or str(x).lower() in {"", "none", "na"} else float(x)

def _split_comma(s: str | None) -> list[str]:
    return [z.strip() for z in str(s or "").split(",") if z.strip()]

def _strip_version(x: str) -> str:
    return str(x).split(".", 1)[0]

def _parse_gtf_attr(attr_str: str) -> dict[str, str]:
    return dict(re.findall(r'(\w+)\s+"([^"]+)"', attr_str))

def _check_samples(count_idx: pd.Index, meta_idx: pd.Index) -> None:
    miss = [s for s in meta_idx if s not in count_idx]
    extra = [s for s in count_idx if s not in meta_idx]
    if miss:
        raise ValueError(f"Missing samples in counts: {miss[:10]}{' …' if len(miss) > 10 else ''}")
    if extra:
        print(f"⚠  Ignored samples (not in metadata): {extra[:10]}{' …' if len(extra) > 10 else ''}", file=sys.stderr)

def _norm_colname_for_match(x: str) -> str:
    b = Path(str(x)).name
    for ext in [".bam", ".sam", ".cram", ".txt", ".tsv", ".csv"]:
        if b.lower().endswith(ext):
            b = b[: -len(ext)]
            break
    suffixes = [
        ".sorted", ".Sorted", ".Aligned", ".aligned",
        ".Aligned.sortedByCoord.out", ".aligned.sortedByCoord.out",
        "_Aligned.sortedByCoord.out", "_aligned.sortedByCoord.out",
    ]
    for suf in suffixes:
        if b.endswith(suf):
            b = b[: -len(suf)]
            break
    return b

def _match_sample_columns(raw_cols: list[str], sample_ids: list[str]) -> tuple[dict[str, str], list[str]]:
    raw_set = set(raw_cols)
    norm_to_raw: dict[str, str] = {}
    for c in raw_cols:
        norm_to_raw.setdefault(_norm_colname_for_match(c), c)

    mapping: dict[str, str] = {}
    for sid in sample_ids:
        if sid in raw_set:
            mapping[sid] = sid
            continue
        if sid in norm_to_raw:
            mapping[sid] = norm_to_raw[sid]
            continue
        ns = _norm_colname_for_match(sid)
        if ns in norm_to_raw:
            mapping[sid] = norm_to_raw[ns]
            continue

    missing = [s for s in sample_ids if s not in mapping]
    return mapping, missing

def _normalize_key_name(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(s).strip().lower())


def _pick_pcol(df: pd.DataFrame) -> str | None:
    if "padj" in df.columns:
        return "padj"
    if "pvalue" in df.columns:
        return "pvalue"
    return None


# ─────────────────── volcano & heat-map functions ─────────────────
def _volcano(full_df: pd.DataFrame, label: str, pthr: float | None, fcthr: float | None) -> None:
    """
    Volcano basé sur les résultats FULL (non filtrés), coloration par seuils.
    """
    pcol = _pick_pcol(full_df)
    if pcol is None:
        return

    pos = full_df[pcol].dropna()
    pos = pos[pos > 0]
    eps = pos.min() if not pos.empty else 1e-300

    thr_p = 1.0 if pthr is None else float(pthr)
    thr_fc = 0.0 if fcthr is None else float(fcthr)

    dd = full_df.copy()
    dd["neglog10p"] = -np.log10(dd[pcol].replace(0, eps))
    dd["signif"] = (dd[pcol] < thr_p) & (dd["log2FoldChange"].abs() > thr_fc)

    plt.figure(figsize=(5, 4))
    plt.scatter(dd.loc[~dd["signif"], "log2FoldChange"], dd.loc[~dd["signif"], "neglog10p"],
                s=8, alpha=.45, color="#A6ACBC")
    plt.scatter(dd.loc[dd["signif"], "log2FoldChange"], dd.loc[dd["signif"], "neglog10p"],
                s=8, alpha=.75, color="#D62728")

    if pthr is not None:
        plt.axhline(-np.log10(thr_p), ls="--", c="grey")
    if fcthr is not None:
        plt.axvline(thr_fc, ls="--", c="grey")
        plt.axvline(-thr_fc, ls="--", c="grey")

    plt.xlabel("log₂ FC")
    plt.ylabel(f"−log₁₀({pcol})")
    plt.title(f"Volcano {label}")
    plt.tight_layout()
    plt.savefig(f"Volcano_{label}.png", dpi=300)
    plt.close()


def _select_heatmap_genes(full_df: pd.DataFrame, n: int = 50) -> list[str]:
    """
    Sélectionne des gènes pour heatmap :
    - priorise padj si dispo (sinon pvalue)
    - prend up/down équilibré selon |log2FC|
    """
    if full_df.empty or "log2FoldChange" not in full_df.columns:
        return []

    pcol = _pick_pcol(full_df)
    dd = full_df.copy()

    if pcol is not None:
        dd = dd.dropna(subset=[pcol])
        dd = dd.sort_values(pcol, ascending=True)
    else:
        dd = dd.copy()

    dd["absLFC"] = dd["log2FoldChange"].abs()
    dd = dd.sort_values(["absLFC"], ascending=False)

    # équilibrer up/down si possible
    up = dd[dd["log2FoldChange"] > 0]
    down = dd[dd["log2FoldChange"] < 0]

    half = max(1, n // 2)
    genes = []
    if len(up) > 0:
        genes += up.head(half).index.tolist()
    if len(down) > 0:
        genes += down.head(half).index.tolist()

    # compléter si manque
    if len(genes) < n:
        rest = [g for g in dd.index.tolist() if g not in set(genes)]
        genes += rest[: (n - len(genes))]

    return genes[:n]


def _heatmap(ad: sc.AnnData, genes: list[str], label: str) -> None:
    if not genes:
        return

    genes = [g for g in genes if g in ad.var_names]
    if len(genes) < 3:
        return

    mat = (pd.DataFrame(ad[:, genes].X.T, index=genes, columns=ad.obs_names)
             .replace([np.inf, -np.inf], np.nan)
             .dropna(how="any"))
    if mat.shape[0] < 3 or mat.shape[1] < 2:
        return

    mat.to_csv(f"Heatmap_{label}.csv")
    sns.clustermap(mat, z_score=0, cmap="RdYlBu_r",
                   cbar_kws={"label": "z-score (VST)"}) \
       .fig.savefig(f"Heatmap_{label}.png", dpi=300)
    plt.close()


# ───────────────────────── GTF auto-mapping ───────────────────────
def _choose_gtf_id_key(gtf_path: str, feature_ids_sample: set[str], prefer_name: str) -> str | None:
    prefer_norm = _normalize_key_name(prefer_name)

    scores: dict[str, int] = {}
    seen_keys: set[str] = set()

    with Path(gtf_path).open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = _parse_gtf_attr(parts[8])
            for k, v in attrs.items():
                seen_keys.add(k)
                vv = _strip_version(v)
                if vv in feature_ids_sample:
                    scores[k] = scores.get(k, 0) + 1

    if not seen_keys:
        return None

    name_matches = [k for k in seen_keys if _normalize_key_name(k) == prefer_norm]
    if name_matches:
        name_matches.sort(key=lambda k: scores.get(k, 0), reverse=True)
        return name_matches[0]

    best = max(scores.items(), key=lambda kv: kv[1])[0] if scores else None
    if best and scores.get(best, 0) > 0:
        return best

    return None


def _build_gtf_annotation_table(
    gtf_path: str,
    feature_key: str,
    requested_cols: list[str],
    feature_ids_all: pd.Series,
) -> pd.DataFrame:
    feat_ids = feature_ids_all.astype(str).map(_strip_version)
    feat_set_sample = set(feat_ids.sample(min(200, len(feat_ids)), random_state=0).tolist())

    gene_keys: set[str] = set()
    gene_has: set[str] = set()
    bridge_candidates: dict[str, int] = {}
    id_line_has_cols: set[str] = set()

    with Path(gtf_path).open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attrs = _parse_gtf_attr(parts[8])

            if ftype == "gene":
                gene_keys.update(attrs.keys())
                for c in requested_cols:
                    if c in attrs:
                        gene_has.add(c)

            if feature_key in attrs:
                for c in requested_cols:
                    if c in attrs:
                        id_line_has_cols.add(c)

                for k in attrs.keys():
                    if k != feature_key:
                        bridge_candidates[k] = bridge_candidates.get(k, 0) + 1

    direct_cols = [c for c in requested_cols if c in id_line_has_cols]
    remaining = [c for c in requested_cols if c not in direct_cols]

    bridge_key: str | None = None
    if remaining and gene_has:
        candidates = [k for k in bridge_candidates if k in gene_keys]
        if candidates:
            candidates.sort(key=lambda k: bridge_candidates.get(k, 0), reverse=True)
            bridge_key = candidates[0]

    direct_map: dict[str, dict[str, str]] = {}
    feat_to_bridge: dict[str, str] = {}
    bridge_to_cols: dict[str, dict[str, str]] = {}

    feat_all_set = set(feat_ids.tolist())

    with Path(gtf_path).open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attrs = _parse_gtf_attr(parts[8])

            if feature_key in attrs:
                fid = _strip_version(attrs[feature_key])
                if fid in feat_set_sample or fid in feat_all_set:
                    if direct_cols:
                        dd = direct_map.get(fid, {})
                        for c in direct_cols:
                            if c in attrs and c not in dd:
                                dd[c] = attrs[c]
                        if dd:
                            direct_map[fid] = dd
                    if bridge_key and bridge_key in attrs:
                        feat_to_bridge[fid] = _strip_version(attrs[bridge_key])

            if bridge_key and remaining and ftype == "gene" and bridge_key in attrs:
                bid = _strip_version(attrs[bridge_key])
                dd = bridge_to_cols.get(bid, {})
                for c in remaining:
                    if c in attrs and c not in dd:
                        dd[c] = attrs[c]
                if dd:
                    bridge_to_cols[bid] = dd

    out = pd.DataFrame({"__feature_id__": feat_ids})

    for c in direct_cols:
        out[c] = out["__feature_id__"].map(lambda k: direct_map.get(k, {}).get(c, ""))

    if bridge_key and remaining:
        bridge_series = out["__feature_id__"].map(lambda k: feat_to_bridge.get(k, ""))
        for c in remaining:
            out[c] = bridge_series.map(lambda bid: bridge_to_cols.get(bid, {}).get(c, ""))

    return out.drop_duplicates(subset=["__feature_id__"]).reset_index(drop=True)


# ─────────────────────────── MultiDE class ────────────────────────
class MultiDE:
    def __init__(
        self,
        counts_fp: str | Path,
        meta_fp: str | Path,
        col_sample: str,
        col_feature: str,
        col_cond: str,
        ctrl_level: str,
        col_time: str | None,
        col_group: str | None,
        extras: list[str],
        pthr: float | None,
        fcthr: float | None,
        gtf_path: str | None = None,
        annotation_columns: list[str] | None = None,
    ):
        self.counts_fp, self.meta_fp = Path(counts_fp), Path(meta_fp)
        self.c = dict(Sample=col_sample, Feature=col_feature,
                      Cond=col_cond, Time=col_time, Group=col_group)
        self.ctrl = ctrl_level
        self.extras = extras
        self.pthr = pthr
        self.fcthr = fcthr
        self.gtf_path = gtf_path
        self.annot_cols = annotation_columns or []

        self.dds: DeseqDataSet | None = None
        self.feature_annot: pd.DataFrame | None = None
        self.sum_time: list[tuple] = []
        self.sum_group: list[tuple] = []

    def prepare(self):
        raw = pd.read_csv(self.counts_fp, sep=None, engine="python")
        feat_col = self.c["Feature"]
        if feat_col not in raw.columns:
            raise ValueError(f"Column '{feat_col}' missing in count table")

        meta = pd.read_csv(self.meta_fp, sep="\t", comment="#").set_index(self.c["Sample"])
        if self.c["Cond"] not in meta.columns:
            raise ValueError(f"condition_column '{self.c['Cond']}' not in metadata")

        sample_ids = [str(s) for s in meta.index.tolist()]
        raw_cols = [str(c) for c in raw.columns.tolist()]
        colmap, missing = _match_sample_columns(raw_cols, sample_ids)
        if missing:
            raise ValueError(
                "Could not match these metadata samples to count-table headers: "
                f"{missing[:10]}{' …' if len(missing) > 10 else ''}"
            )

        counts_df = raw[[feat_col] + [colmap[s] for s in sample_ids]].copy()
        counts_df.rename(columns={colmap[s]: s for s in sample_ids}, inplace=True)

        counts_df[feat_col] = counts_df[feat_col].astype(str)
        if counts_df[feat_col].duplicated().any():
            dup_n = int(counts_df[feat_col].duplicated().sum())
            raise ValueError(
                f"'{feat_col}' contains duplicate identifiers ({dup_n}). "
                "DESeq2 requires unique feature IDs."
            )

        count_mtx = (
            counts_df.set_index(feat_col)
                     .apply(pd.to_numeric, errors="coerce")
                     .fillna(0)
                     .astype(int)
        )
        count_mtx = count_mtx.loc[count_mtx.sum(axis=1) > 0]
        count_mtx = count_mtx.T
        _check_samples(count_mtx.index, meta.index)
        count_mtx = count_mtx.loc[meta.index]

        meta[self.c["Cond"]] = meta[self.c["Cond"]].astype("category")
        for col in [self.c["Time"], self.c["Group"], *self.extras]:
            if col and col in meta.columns:
                meta[col] = meta[col].astype("category")

        terms = [c for c in [self.c["Time"], self.c["Group"], *self.extras]
                 if c and c in meta.columns] + [self.c["Cond"]]
        design = "~ " + " + ".join(terms)

        self.feature_annot = self._build_feature_annotation(raw, sample_ids)

        self.dds = DeseqDataSet(counts=count_mtx, metadata=meta, design=design)
        self.dds.deseq2()

    def _build_feature_annotation(self, raw: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame | None:
        feat_col = self.c["Feature"]

        raw_cols = [str(c) for c in raw.columns.tolist()]
        colmap, _ = _match_sample_columns(raw_cols, [str(s) for s in sample_ids])
        sample_cols_in_raw = set(colmap.values())

        base_cols = [c for c in raw.columns if c not in sample_cols_in_raw]
        if feat_col not in base_cols:
            base_cols = [feat_col] + [c for c in base_cols if c != feat_col]

        annot = raw[base_cols].copy()
        annot[feat_col] = annot[feat_col].astype(str)
        annot = annot.drop_duplicates(subset=[feat_col]).reset_index(drop=True)

        if not (self.gtf_path and self.annot_cols):
            return annot

        requested = [c for c in self.annot_cols if c]
        missing = [c for c in requested if c not in annot.columns]
        if not missing:
            return annot

        feat_ids_series = annot[feat_col].astype(str).map(_strip_version)
        feat_sample_set = set(feat_ids_series.sample(min(200, len(feat_ids_series)), random_state=0).tolist())

        gtf_key = _choose_gtf_id_key(self.gtf_path, feat_sample_set, prefer_name=feat_col)
        if not gtf_key:
            warnings.warn(
                "Could not identify a GTF attribute key matching your count-table feature IDs. "
                "Requested GTF annotation columns will be empty."
            )
            for c in missing:
                annot[c] = ""
            return annot

        gtf_annot = _build_gtf_annotation_table(
            gtf_path=self.gtf_path,
            feature_key=gtf_key,
            requested_cols=missing,
            feature_ids_all=annot[feat_col].astype(str),
        )

        left = annot.copy()
        left["__k__"] = left[feat_col].astype(str).map(_strip_version)
        right = gtf_annot.copy()
        right.rename(columns={"__feature_id__": "__k__"}, inplace=True)

        merged = pd.merge(left, right, on="__k__", how="left", suffixes=("", "_gtf"))
        merged.drop(columns=["__k__"], inplace=True, errors="ignore")

        for c in missing:
            if c not in merged.columns:
                merged[c] = ""

        return merged.reset_index(drop=True)

    def _filter_df(self, df_full: pd.DataFrame) -> pd.DataFrame:
        """
        Filtrage sur padj si dispo, sinon pvalue.
        Applique aussi |log2FC| si demandé.
        """
        df = df_full.copy()
        pcol = _pick_pcol(df)

        if self.pthr is not None and pcol is not None:
            df = df[df[pcol] < float(self.pthr)]

        if self.fcthr is not None and "log2FoldChange" in df.columns:
            df = df[df["log2FoldChange"].abs() > float(self.fcthr)]

        return df

    def _format_and_merge_annot(self, res: pd.DataFrame) -> pd.DataFrame:
        feat_col = self.c["Feature"]

        out = res.copy()
        id_values = [str(i) for i in out.index]

        if feat_col in out.columns:
            out[feat_col] = id_values
            cols = [feat_col] + [c for c in out.columns if c != feat_col]
            out = out[cols]
        else:
            out.insert(0, feat_col, id_values)

        out = out.reset_index(drop=True)

        if self.feature_annot is not None and not self.feature_annot.empty:
            out = pd.merge(out, self.feature_annot, on=feat_col, how="left")
            annot_cols = [c for c in self.feature_annot.columns if c != feat_col]
            rest = [c for c in out.columns if c not in ([feat_col] + annot_cols)]
            out = out[[feat_col] + annot_cols + rest]

        return out

    def _contrast(
        self,
        tr: str,
        ref: str,
        label: str,
        dds: DeseqDataSet | None = None,
        ctx: dict[str, str] | None = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        ds = dds or self.dds
        st = DeseqStats(ds, contrast=(self.c["Cond"], tr, ref))
        st.summary(save=False)

        df_full = st.get_results_df() if hasattr(st, "get_results_df") else st.results_df
        df_full = df_full.copy()

        if ctx:
            for k, v in ctx.items():
                if k not in df_full.columns:
                    df_full.insert(0, k, v)

        df_filt = self._filter_df(df_full)

        # Exports (FULL then FILT), with annotations
        full_out = self._format_and_merge_annot(df_full)
        filt_out = self._format_and_merge_annot(df_filt)

        full_out.sort_values("log2FoldChange", ascending=False).to_csv(f"DE_FULL_{label}.csv", index=False)
        filt_out.sort_values("log2FoldChange", ascending=False).to_csv(f"DE_FILT_{label}.csv", index=False)

        # Volcano from FULL
        _volcano(df_full, label, self.pthr, self.fcthr)

        return df_full, df_filt

    def _dds_pooled(self, mask: pd.Series) -> DeseqDataSet:
        cts = pd.DataFrame(self.dds.X[mask.values, :],
                           index=self.dds.obs_names[mask],
                           columns=self.dds.var_names)
        meta = self.dds.obs.loc[mask].copy()
        meta[self.c["Cond"]] = pd.Categorical(
            [self.ctrl if x == self.ctrl else "ALL" for x in meta[self.c["Cond"]].astype(str)],
            categories=[self.ctrl, "ALL"],
        )
        dds_p = DeseqDataSet(counts=cts, metadata=meta, design=f"~ {self.c['Cond']}")
        dds_p.deseq2()
        return dds_p

    def _vst_pca_heat(self, dfs_full: list[pd.DataFrame], labels: list[str]) -> None:
        self.dds.vst_fit(use_design=False)
        ad = sc.AnnData(
            self.dds.vst_transform(),
            obs=self.dds.obs.copy(),
            var=pd.DataFrame(index=self.dds.var_names),
        )

        sc.pp.pca(ad, n_comps=min(30, ad.n_obs - 1, ad.n_vars - 1))
        (pd.DataFrame(ad.obsm["X_pca"][:, :2],
                      index=ad.obs_names, columns=["PC1", "PC2"])
           .join(ad.obs).to_csv("pca_metadata.csv"))

        v1, v2 = ad.uns["pca"]["variance_ratio"][:2]
        pd.DataFrame({"PC1 Proportion": [v1], "PC2 Proportion": [v2]}).to_csv("pca_proportions.csv", index=False)

        # Heatmaps: choose genes from FULL (prioritizing padj / absLFC)
        for df_full, lbl in zip(dfs_full, labels):
            genes = _select_heatmap_genes(df_full, n=50)
            if genes:
                _heatmap(ad, genes, lbl)

    def _stratified(self, col: str | None, store: list[tuple], tag: str, pooled: bool = False):
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
                metadata=self.dds.obs.loc[mask],
                design=design,
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
                df_full, df_filt = self._contrast(tr, self.ctrl, lbl, sub_dds, {col: lvl})

                # summary basé sur FILT (ce que les gens attendent en général)
                store.append((
                    lvl, tr, len(df_filt),
                    int((df_filt.log2FoldChange > 0).sum()) if not df_filt.empty else 0,
                    int((df_filt.log2FoldChange < 0).sum()) if not df_filt.empty else 0,
                ))

            if pooled:
                try:
                    pool_dds = self._dds_pooled(mask)
                    lbl = f"ALL_vs_{self.ctrl}_{tag}{lvl}"
                    df_full, df_filt = self._contrast("ALL", self.ctrl, lbl, pool_dds, {col: lvl})
                    store.append((
                        lvl, "ALL", len(df_filt),
                        int((df_filt.log2FoldChange > 0).sum()) if not df_filt.empty else 0,
                        int((df_filt.log2FoldChange < 0).sum()) if not df_filt.empty else 0,
                    ))
                except Exception as e:
                    warnings.warn(f"{col} {lvl} pooled: {e}")

    @staticmethod
    def _write(rows: Sequence[tuple], csv_path: str, level_name: str):
        if rows:
            pd.DataFrame(rows, columns=[level_name, "Treatment", "DEGs_total", "DEGs_up", "DEGs_down"]) \
              .to_csv(csv_path, index=False)

    def run(self):
        self.prepare()

        dfs_full: list[pd.DataFrame] = []
        labels: list[str] = []

        for tr in self.dds.obs[self.c["Cond"]].cat.categories:
            if tr == self.ctrl:
                continue
            df_full, _df_filt = self._contrast(tr, self.ctrl, f"{tr}_vs_{self.ctrl}")
            dfs_full.append(df_full)
            labels.append(f"{tr}_vs_{self.ctrl}")

        if len(self.dds.obs[self.c["Cond"]].cat.categories) > 1:
            dds_pool = self._dds_pooled(self.dds.obs[self.c["Cond"]].notna())
            df_full, _df_filt = self._contrast("ALL", self.ctrl, "ALL_vs_{0}".format(self.ctrl), dds_pool)
            dfs_full.append(df_full)
            labels.append(f"ALL_vs_{self.ctrl}")

        self._stratified(self.c["Time"], self.sum_time, "T", True)
        self._stratified(self.c["Group"], self.sum_group, "G", True)
        if self.sum_time:
            self._write(self.sum_time, "Summary_Timepoint.csv", self.c["Time"])
        if self.sum_group:
            self._write(self.sum_group, "Summary_Group.csv", self.c["Group"])

        self._vst_pca_heat(dfs_full, labels)


# ───────────────────────────── CLI ────────────────────────────────
if __name__ == "__main__":
    (
        _, counts, meta,
        col_sample, col_feature,
        col_cond, ctrl,
        col_time, col_group,
        extras,
        p_thr, fc_thr, *rest
    ) = sys.argv

    gtf_path = rest[0] if len(rest) >= 1 and rest[0] not in {"None", ""} else None
    annot_cols = _split_comma(rest[1]) if len(rest) >= 2 else []

    MultiDE(
        counts_fp=counts,
        meta_fp=meta,
        col_sample=col_sample,
        col_feature=col_feature,
        col_cond=col_cond,
        ctrl_level=ctrl,
        col_time=None if col_time == "None" else col_time,
        col_group=None if col_group == "None" else col_group,
        extras=_split_comma(extras),
        pthr=_parse_float(p_thr),
        fcthr=_parse_float(fc_thr),
        gtf_path=gtf_path,
        annotation_columns=annot_cols,
    ).run()