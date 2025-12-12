from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path
from typing import List, Optional, Sequence, Dict, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# ───────────────────────── utilities ──────────────────────────────
def _parse_float(x: str | float | None) -> Optional[float]:
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
    # try to match "gene_ID" -> "gene_id", "TranscriptID" -> "transcriptid"
    return re.sub(r"[^a-z0-9]+", "_", str(s).strip().lower())


# ─────────────────── volcano & heat-map functions ─────────────────
def _volcano(d: pd.DataFrame, label: str,
             pthr: Optional[float], fcthr: Optional[float]) -> None:
    pcol = "padj" if "padj" in d.columns else ("pvalue" if "pvalue" in d.columns else None)
    if pcol is None:
        return
    pos = d[pcol].dropna()
    pos = pos[pos > 0]
    eps = pos.min() if not pos.empty else 1e-300

    dd = d.assign(
        neglog10p=-np.log10(d[pcol].replace(0, eps)),
        signif=((d[pcol] < (pthr or 1)) & (d.log2FoldChange.abs() > (fcthr or 0)))
    )
    plt.figure(figsize=(5, 4))
    plt.scatter(dd.log2FoldChange[~dd.signif], dd.neglog10p[~dd.signif],
                s=8, alpha=.45, color="#A6ACBC")
    plt.scatter(dd.log2FoldChange[dd.signif], dd.neglog10p[dd.signif],
                s=8, alpha=.75, color="#D62728")
    if pthr:
        plt.axhline(-np.log10(pthr), ls="--", c="grey")
    if fcthr:
        plt.axvline(fcthr, ls="--", c="grey")
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


# ───────────────────────── GTF auto-mapping ───────────────────────
def _choose_gtf_id_key(gtf_path: str, feature_ids_sample: set[str], prefer_name: str) -> Optional[str]:
    """
    Choose which GTF attribute key best matches the IDs used in the count table.
    Uses:
      1) name match (case-insensitive normalized)
      2) value overlap score with sample IDs
    """
    prefer_norm = _normalize_key_name(prefer_name)

    # Pass 1: collect key scores by matching values
    scores: Dict[str, int] = {}
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

    # Strong preference: if a key name matches genes_colname (normalized), use it
    name_matches = [k for k in seen_keys if _normalize_key_name(k) == prefer_norm]
    if name_matches:
        # if multiple, take the one with highest overlap
        name_matches.sort(key=lambda k: scores.get(k, 0), reverse=True)
        return name_matches[0]

    # Otherwise pick the highest overlap key
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
    """
    Build DF keyed by the count-table feature IDs (feature_key in GTF),
    adding requested_cols. Works for gene-level or transcript-level, and if
    requested columns are only available on gene features, it auto-selects a bridge key.
    """
    # normalize IDs
    feat_ids = feature_ids_all.astype(str).map(_strip_version)
    feat_set_sample = set(feat_ids.sample(min(200, len(feat_ids)), random_state=0).tolist())

    # Pass A: discover what columns exist where, and candidate bridge keys
    gene_keys: set[str] = set()
    gene_has: set[str] = set()
    bridge_candidates: Dict[str, int] = {}  # key -> count on ID-lines
    id_line_has_cols: set[str] = set()
    id_lines_seen = 0

    with Path(gtf_path).open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attrs = _parse_gtf_attr(parts[8])

            # gene feature lines
            if ftype == "gene":
                gene_keys.update(attrs.keys())
                for c in requested_cols:
                    if c in attrs:
                        gene_has.add(c)

            # lines where the chosen feature_key appears
            if feature_key in attrs:
                id_lines_seen += 1
                for c in requested_cols:
                    if c in attrs:
                        id_line_has_cols.add(c)

                # candidate bridges: any attribute key on this line that also exists on gene lines
                for k in attrs.keys():
                    if k != feature_key:
                        bridge_candidates[k] = bridge_candidates.get(k, 0) + 1

    # Decide what we can get directly from ID lines
    direct_cols = [c for c in requested_cols if c in id_line_has_cols]

    # Remaining columns might need a bridge through gene features
    remaining = [c for c in requested_cols if c not in direct_cols]

    # choose a bridge key only if needed and possible
    bridge_key: Optional[str] = None
    if remaining and gene_has:
        # candidate must be present on gene lines too, and reasonably present on ID lines
        candidates = [k for k in bridge_candidates.keys() if k in gene_keys]
        if candidates:
            # score = presence on ID-lines + ability to retrieve remaining cols from gene
            candidates.sort(key=lambda k: bridge_candidates.get(k, 0), reverse=True)
            bridge_key = candidates[0]

    # Pass B: actually build mappings
    # direct: feature_id -> {cols}
    direct_map: Dict[str, Dict[str, str]] = {}
    # bridge: feature_id -> bridge_id
    feat_to_bridge: Dict[str, str] = {}
    # gene: bridge_id -> {remaining cols}
    bridge_to_cols: Dict[str, Dict[str, str]] = {}

    with Path(gtf_path).open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attrs = _parse_gtf_attr(parts[8])

            # direct annotations on ID lines
            if feature_key in attrs:
                fid = _strip_version(attrs[feature_key])
                if fid in feat_set_sample or fid in set(feat_ids.tolist()):
                    if direct_cols:
                        dd = direct_map.get(fid, {})
                        for c in direct_cols:
                            if c in attrs and c not in dd:
                                dd[c] = attrs[c]
                        if dd:
                            direct_map[fid] = dd
                    if bridge_key and bridge_key in attrs:
                        feat_to_bridge[fid] = _strip_version(attrs[bridge_key])

            # gene lines for remaining cols via bridge key
            if bridge_key and remaining and ftype == "gene" and bridge_key in attrs:
                bid = _strip_version(attrs[bridge_key])
                dd = bridge_to_cols.get(bid, {})
                for c in remaining:
                    if c in attrs and c not in dd:
                        dd[c] = attrs[c]
                if dd:
                    bridge_to_cols[bid] = dd

    # assemble final DF
    out = pd.DataFrame({ "__feature_id__": feat_ids })

    # apply direct cols
    for c in direct_cols:
        out[c] = out["__feature_id__"].map(lambda k: direct_map.get(k, {}).get(c, ""))

    # apply bridged cols if any
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
        col_feature: str,   # genes_colname (user-defined column name in count table)
        col_cond: str,
        ctrl_level: str,
        col_time: str | None,
        col_group: str | None,
        extras: list[str],
        pthr: Optional[float],
        fcthr: Optional[float],
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

        # Build counts table feature + sample columns
        counts_df = raw[[feat_col] + [colmap[s] for s in sample_ids]].copy()
        counts_df.rename(columns={colmap[s]: s for s in sample_ids}, inplace=True)

        counts_df[feat_col] = counts_df[feat_col].astype(str)
        if counts_df[feat_col].duplicated().any():
            dup_n = int(counts_df[feat_col].duplicated().sum())
            raise ValueError(
                f"'{feat_col}' contains duplicate identifiers ({dup_n}). "
                "DESeq2 requires unique feature IDs."
            )

        # numeric matrix features x samples -> transpose
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

        # Build annotation table: keep all non-sample columns + optional GTF columns
        self.feature_annot = self._build_feature_annotation(raw, sample_ids)

        self.dds = DeseqDataSet(counts=count_mtx, metadata=meta, design=design)
        self.dds.deseq2()

    def _build_feature_annotation(self, raw: pd.DataFrame, sample_ids: list[str]) -> pd.DataFrame | None:
        feat_col = self.c["Feature"]

        # Identify which raw columns are samples (so we keep everything else)
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

        # add only missing requested columns from GTF
        requested = [c for c in self.annot_cols if c]
        missing = [c for c in requested if c not in annot.columns]
        if not missing:
            return annot

        # auto-detect which GTF attribute key matches the feature IDs used in counts
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
        # merge: join on stripped IDs
        left = annot.copy()
        left["__k__"] = left[feat_col].astype(str).map(_strip_version)
        right = gtf_annot.copy()
        right.rename(columns={"__feature_id__": "__k__"}, inplace=True)

        merged = pd.merge(left, right, on="__k__", how="left", suffixes=("", "_gtf"))
        merged.drop(columns=["__k__"], inplace=True, errors="ignore")

        # ensure all missing columns exist
        for c in missing:
            if c not in merged.columns:
                merged[c] = ""

        return merged.reset_index(drop=True)

    def _contrast(self, tr: str, ref: str, label: str,
                  dds: DeseqDataSet | None = None,
                  ctx: dict[str, str] | None = None) -> pd.DataFrame:
        ds = dds or self.dds
        st = DeseqStats(ds, contrast=(self.c["Cond"], tr, ref))
        st.summary(save=False)
        df = st.get_results_df() if hasattr(st, "get_results_df") else st.results_df

        if self.pthr is not None:
            df = df[df.pvalue < self.pthr]
        if self.fcthr is not None:
            df = df[df.log2FoldChange.abs() > self.fcthr]
        df = df.copy()

        if ctx:
            for k, v in ctx.items():
                if k not in df.columns:
                    df.insert(0, k, v)

        feat_col = self.c["Feature"]

        out = df.copy()
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

        out.sort_values("log2FoldChange", ascending=False).to_csv(f"DE_{label}.csv", index=False)
        _volcano(df, label, self.pthr, self.fcthr)
        return df

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
        pd.DataFrame({"PC1 Proportion": [v1], "PC2 Proportion": [v2]}).to_csv("pca_proportions.csv", index=False)

        for d, l in zip(dfs, labels):
            if not d.empty:
                _heatmap(ad, d.nlargest(50, "log2FoldChange").index.tolist(), l)

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
                df = self._contrast(tr, self.ctrl, lbl, sub_dds, {col: lvl})
                store.append((lvl, tr, len(df),
                              (df.log2FoldChange > 0).sum(),
                              (df.log2FoldChange < 0).sum()))

            if pooled:
                try:
                    pool_dds = self._dds_pooled(mask)
                    lbl = f"ALL_vs_{self.ctrl}_{tag}{lvl}"
                    df = self._contrast("ALL", self.ctrl, lbl, pool_dds, {col: lvl})
                    store.append((lvl, "ALL", len(df),
                                  (df.log2FoldChange > 0).sum(),
                                  (df.log2FoldChange < 0).sum()))
                except Exception as e:
                    warnings.warn(f"{col} {lvl} pooled: {e}")

    @staticmethod
    def _write(rows: Sequence[tuple], csv_path: str, level_name: str):
        if rows:
            pd.DataFrame(rows, columns=[level_name, "Treatment", "DEGs_total", "DEGs_up", "DEGs_down"]) \
              .to_csv(csv_path, index=False)

    def run(self):
        self.prepare()
        dfs, labels = [], []

        for tr in self.dds.obs[self.c["Cond"]].cat.categories:
            if tr == self.ctrl:
                continue
            dfs.append(self._contrast(tr, self.ctrl, f"{tr}_vs_{self.ctrl}"))
            labels.append(f"{tr}_vs_{self.ctrl}")

        if len(self.dds.obs[self.c["Cond"]].cat.categories) > 1:
            dds_pool = self._dds_pooled(self.dds.obs[self.c["Cond"]].notna())
            dfs.append(self._contrast("ALL", self.ctrl, f"ALL_vs_{self.ctrl}", dds_pool))
            labels.append(f"ALL_vs_{self.ctrl}")

        self._stratified(self.c["Time"], self.sum_time, 'T', True)
        self._stratified(self.c["Group"], self.sum_group, 'G', True)
        if self.sum_time:
            self._write(self.sum_time, "Summary_Timepoint.csv", self.c["Time"])
        if self.sum_group:
            self._write(self.sum_group, "Summary_Group.csv", self.c["Group"])

        self._vst_pca_heat(dfs, labels)


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

    # Optional:
    # rest[0] = gtf_path (or "None")
    # rest[1] = annotation_columns (comma-separated)
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
