#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GC Content Calculator — window=100 nt (classique) + highlight de la sélection du 1er plot dans la séquence.

- Lit une séquence DNA/RNA (FASTA ou texte).
- Détecte DNA/RNA automatiquement.
- Calcule GC% par fenêtre de 100 nt (pas=1), CpG/GpC par 100 nt, et prédit les îlots CpG (200 nt).
- Génère un HTML autonome avec 2 sous-graphes Plotly + un bloc "Sequence".
- NOUVEAU : quand on sélectionne une région sur le 1er sous-graphe (GC%), la séquence correspondante est surlignée.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict

import plotly.graph_objects as go
from plotly.subplots import make_subplots

IUPAC_DNA = set("ACGTUacgtuNn")  # non-IUPAC retirés

# ---------------------------
# IO helpers
# ---------------------------

def read_fasta_or_plain(path: Path) -> str:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return ""
    if lines[0].startswith(">"):
        seq = "".join(ln for ln in lines[1:] if not ln.startswith(">"))
    else:
        seq = "".join(lines)
    seq = "".join(ch for ch in seq if ch in IUPAC_DNA)
    return seq.upper()

def auto_type(seq: str) -> str:
    return "RNA" if "U" in seq and "T" not in seq else "DNA"

# ---------------------------
# Core calculations
# ---------------------------

def rolling_counts(seq: str, win: int) -> Tuple[List[int], List[int], List[int], List[int]]:
    """Return arrays (A, C, G, T_or_U) counts in each window."""
    n = len(seq)
    if n < win:
        s = seq
        return [s.count("A")], [s.count("C")], [s.count("G")], [s.count("T") + s.count("U")]
    A = [0]*(n - win + 1)
    C = [0]*(n - win + 1)
    G = [0]*(n - win + 1)
    TU = [0]*(n - win + 1)

    w = seq[:win]
    A0 = w.count("A"); C0 = w.count("C"); G0 = w.count("G"); TU0 = w.count("T") + w.count("U")
    A[0] = A0; C[0] = C0; G[0] = G0; TU[0] = TU0

    for i in range(1, n - win + 1):
        out = seq[i-1]
        inn = seq[i+win-1]
        if out == "A": A0 -= 1
        elif out == "C": C0 -= 1
        elif out == "G": G0 -= 1
        elif out in ("T","U"): TU0 -= 1
        if inn == "A": A0 += 1
        elif inn == "C": C0 += 1
        elif inn == "G": G0 += 1
        elif inn in ("T","U"): TU0 += 1
        A[i] = A0; C[i] = C0; G[i] = G0; TU[i] = TU0
    return A, C, G, TU

def rolling_dinuc_counts(seq: str, win: int, dinuc: str) -> List[int]:
    """Count overlapping dinucleotides within each window."""
    n = len(seq)
    if n < win:
        return [sum(1 for i in range(len(seq)-1) if seq[i:i+2] == dinuc)]
    hits = [0]*(n-1)
    for i in range(n-1):
        if seq[i:i+2] == dinuc:
            hits[i] = 1
    k = win - 1
    out = []
    cur = sum(hits[:k])
    out.append(cur)
    for i in range(1, (n-1) - k + 1):
        cur -= hits[i-1]
        cur += hits[i+k-1]
        out.append(cur)
    return out

@dataclass
class CpGIsland:
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive
    gc: float
    obs_exp: float

def predict_cpg_islands(seq: str) -> List[CpGIsland]:
    """Classic criteria: len>=200, GC%>=50, Obs/Exp CpG>=0.6"""
    n = len(seq)
    win = 200
    if n < win:
        return []
    A, C, G, TU = rolling_counts(seq, win)
    cpg = rolling_dinuc_counts(seq, win, "CG")

    islands: List[CpGIsland] = []
    for i in range(len(A)):
        c = C[i]; g = G[i]
        gc = 100.0 * (c + g) / win
        expected = (c * g) / win if win else 0.0
        obs = cpg[i]
        obs_exp = (obs / expected) if expected > 0 else 0.0
        if gc >= 50.0 and obs_exp >= 0.6:
            start = i
            end = i + win
            islands.append(CpGIsland(start, end, gc, obs_exp))
    if not islands:
        return []
    merged: List[CpGIsland] = []
    cur = islands[0]
    for isl in islands[1:]:
        if isl.start <= cur.end:
            cur = CpGIsland(cur.start, max(cur.end, isl.end), max(cur.gc, isl.gc), max(cur.obs_exp, isl.obs_exp))
        else:
            merged.append(cur)
            cur = isl
    merged.append(cur)
    return merged

# ---------------------------
# Plot/colors
# ---------------------------

def color_band(gc: float) -> str:
    # buckets: <20, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, >=80
    if gc < 20:   return "rgba(90,140,200,0.9)"
    if gc < 30:   return "rgba(120,170,210,0.9)"
    if gc < 40:   return "rgba(160,200,140,0.95)"
    if gc < 50:   return "rgba(210,220,120,0.95)"
    if gc < 60:   return "rgba(240,170,80,0.95)"
    if gc < 70:   return "rgba(230,120,60,0.95)"
    if gc < 80:   return "rgba(200,70,60,0.95)"
    return "rgba(140,30,30,0.95)"

# ---------------------------
# Sequence HTML
# ---------------------------

def make_sequence_html(seq: str, line=100, group=10) -> str:
    """Monospaced colored sequence block with coordinates every line.

    Chaque chunk (10 nt) reçoit data-start/data-end (1-based) pour pouvoir
    surligner lors d'une sélection sur le graphe.
    """
    out = []
    n = len(seq)
    out.append("""
<div id="seqPanel" style="font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size:12px;">
  <div style="margin-bottom:6px;color:#555">Full length: %d</div>
""" % n)
    pos1 = 1  # 1-based
    for row_start in range(0, n, line):
        row_end = min(n, row_start + line)
        row = seq[row_start:row_end]
        out.append(f'<div><span style="display:inline-block;width:50px;color:#666">{row_start+1}</span>')
        for g in range(0, len(row), group):
            chunk = row[g:g+group]
            c_len = len(chunk)
            start = pos1
            end = pos1 + c_len - 1
            # style par défaut (vert pâle), classe seq-chunk pour pouvoir toggler .sel
            out.append(
                '<span class="seq-chunk" '
                f'data-start="{start}" data-end="{end}" '
                'style="background:#d7ffd7;border-radius:3px;margin:1px;'
                'padding:2px 3px;display:inline-block">'
                f'{chunk}</span>'
            )
            pos1 += c_len
        out.append("</div>")
    out.append("</div>")
    return "\n".join(out)

# ---------------------------
# HTML assembly
# ---------------------------

def make_html(payload: Dict, out_html: Path) -> None:
    seq = payload["seq"]
    seq_type = payload["type"]
    win = payload["win"]              # 100
    x_gc = payload["x_gc"]            # window start (1-based for display)
    gc_pct = payload["gc_pct"]
    cpg = payload["cpg_per100"]
    gpc = payload["gpc_per100"]
    islands: List[CpGIsland] = payload["islands"]

    # Downsample (visuel identique)
    stride = max(1, len(x_gc) // 3000)
    x_gc_ds = x_gc[::stride]
    gc_pct_ds = gc_pct[::stride]
    band_colors = [color_band(v) for v in gc_pct_ds]

    fig = make_subplots(
        rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.07,
        subplot_titles=(f"GC% (window = {win} nt) · {seq_type}",
                        "CpG & GpC per 100 nt · CpG Islands")
    )

    # Row 1: GC% bar colorée
    fig.add_trace(
        go.Bar(
            x=x_gc_ds,
            y=gc_pct_ds,
            marker=dict(color=band_colors),
            hovertemplate="pos: %{x}–%{customdata}<br>GC: %{y:.2f}%<extra></extra>",
            customdata=[x + win - 1 for x in x_gc_ds],
            name="GC%"
        ),
        row=1, col=1
    )

    # Row 2: CpG / GpC
    stride2 = max(1, len(cpg) // 3000)
    xs2 = x_gc[::stride2]
    fig.add_trace(
        go.Scatter(
            x=xs2, y=cpg[::stride2],
            mode="lines",
            name="CpG / 100nt",
            hovertemplate="pos: %{x}–%{customdata}<br>CpG: %{y:.3g}<extra></extra>",
            customdata=[x + win - 1 for x in xs2]
        ),
        row=2, col=1
    )
    fig.add_trace(
        go.Scatter(
            x=xs2, y=gpc[::stride2],
            mode="lines",
            name="GpC / 100nt",
            hovertemplate="pos: %{x}–%{customdata}<br>GpC: %{y:.3g}<extra></extra>",
            customdata=[x + win - 1 for x in xs2]
        ),
        row=2, col=1
    )

    # Îlots CpG (barres noires)
    for isl in islands:
        fig.add_shape(
            type="rect",
            x0=isl.start+1, x1=isl.end, y0=-0.5, y1=0.5,
            line=dict(color="black", width=0),
            fillcolor="black",
            opacity=1.0,
            row=2, col=1
        )

    fig.update_yaxes(title_text="GC %", row=1, col=1, range=[0, 100])
    fig.update_yaxes(title_text="Count/100nt", row=2, col=1)
    fig.update_xaxes(title_text="Position (nt)", row=2, col=1)

    fig.update_layout(
        height=720,
        bargap=0,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0.01),
        margin=dict(l=60, r=30, t=60, b=60),
        dragmode="select"  # permettre la sélection à la souris
    )

    # Rendu avec un ID fixe pour brancher les handlers
    fig_html = fig.to_html(include_plotlyjs="cdn", full_html=False, div_id="fig")

    # Bloc séquence
    seq_panel = make_sequence_html(seq)

    # CSS léger pour le highlight
    css_extra = """
<style>
.seq-chunk.sel{
  background:#ffe08a !important;  /* surbrillance douce */
  box-shadow: inset 0 0 0 1px #cc9a00;
}
</style>
"""

    # JS pour relier sélection -> highlight
    js_link = """
<script>
(function(){
  const fig = document.getElementById('fig');
  const chunks = Array.from(document.querySelectorAll('#seqPanel .seq-chunk'));

  function highlightRange(x0, x1){
    if (x0 == null || x1 == null) { clear(); return; }
    const a = Math.min(x0, x1);
    const b = Math.max(x0, x1);
    for (const el of chunks){
      const s = +el.dataset.start;
      const e = +el.dataset.end;
      // chevauchement ?
      if (e >= a && s <= b) el.classList.add('sel');
      else el.classList.remove('sel');
    }
  }

  function clear(){
    for (const el of chunks){ el.classList.remove('sel'); }
  }

  // Sélection rectangulaire (dragmode='select')
  fig.on('plotly_selected', function(ev){
    if(!ev || !ev.points || !ev.points.length){ clear(); return; }
    // On s'intéresse aux points de la trace 0 (GC% bar)
    let minX = Infinity, maxX = -Infinity;
    for(const p of ev.points){
      if(p.curveNumber !== 0) continue; // GC% = trace 0
      const start = +p.x;
      const end = (p.customdata != null) ? +p.customdata : start;
      if (isFinite(start)) { minX = Math.min(minX, start); maxX = Math.max(maxX, start); }
      if (isFinite(end))   { minX = Math.min(minX, end);   maxX = Math.max(maxX, end); }
    }
    if(minX === Infinity || maxX === -Infinity) { clear(); return; }
    highlightRange(minX, maxX);
  });

  // Zoom/relayout (p. ex. l'utilisateur zoome au lieu de sélectionner)
  fig.on('plotly_relayout', function(ev){
    const x0 = ev['xaxis.range[0]'];
    const x1 = ev['xaxis.range[1]'];
    if (x0 != null && x1 != null) highlightRange(+x0, +x1);
  });

  // Double-clic pour nettoyer
  fig.on('plotly_doubleclick', function(){ clear(); });
})();
</script>
"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>GC Content Report</title>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
{css_extra}
</head>
<body>
  <div class="card" style="border:1px solid #e5e7eb;border-radius:12px;padding:16px;margin-bottom:16px;box-shadow:0 1px 2px rgba(0,0,0,0.04);">
    <h1 style="font-family:Inter,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;font-size:20px;margin:0 0 8px;">GC Content Calculator</h1>
    <div style="margin:6px 0 10px;">
      <span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:rgba(90,140,200,0.9);margin-right:6px;vertical-align:middle"></span>&lt;20
      <span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:rgba(120,170,210,0.9);margin:0 6px;vertical-align:middle"></span>20–29
      <span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:rgba(160,200,140,0.95);margin:0 6px;vertical-align:middle"></span>30–39
      <span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:rgba(210,220,120,0.95);margin:0 6px;vertical-align:middle"></span>40–49
      <span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:rgba(240,170,80,0.95);margin:0 6px;vertical-align:middle"></span>50–59
      <span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:rgba(230,120,60,0.95);margin:0 6px;vertical-align:middle"></span>60–69
      <span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:rgba(200,70,60,0.95);margin:0 6px;vertical-align:middle"></span>70–79
      <span style="display:inline-block;width:14px;height:14px;border-radius:3px;background:rgba(140,30,30,0.95);margin:0 6px;vertical-align:middle"></span>&ge;80
      <span style="display:inline-block;border-radius:999px;padding:4px 8px;background:#eef2ff;color:#1f2937;font-size:12px;margin-left:6px">Black bars: predicted CpG islands</span>
    </div>
    {fig_html}
  </div>

  <div class="card" style="border:1px solid #e5e7eb;border-radius:12px;padding:16px;margin-bottom:16px;box-shadow:0 1px 2px rgba(0,0,0,0.04);">
    <h2 style="font-family:Inter,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;font-size:16px;margin:0 0 8px;">Sequence</h2>
    {seq_panel}
  </div>

{js_link}
</body>
</html>
"""
    out_html.write_text(html, encoding="utf-8")

# ---------------------------
# Main
# ---------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="GC content (window 100) → HTML + sequence highlight")
    ap.add_argument("--in", dest="in_path", required=True, help="Input FASTA or plain sequence")
    ap.add_argument("--out", dest="out_dir", required=True, help="Output folder")
    ap.add_argument("--prefix", dest="prefix", default="gc_content", help="Output prefix")
    args = ap.parse_args()

    in_path = Path(args.in_path)
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix.strip() or "gc_content"

    seq = read_fasta_or_plain(in_path)
    if not seq:
        raise SystemExit("Empty or unparsable sequence.")

    seq_type = auto_type(seq)
    WIN = 100  # fenêtre fixe

    # Séries de base
    A, C, G, TU = rolling_counts(seq, WIN)
    n = len(seq)
    x_gc = list(range(1, n - WIN + 2)) if n >= WIN else [1]
    gc_pct = [ ((C[i] + G[i]) / WIN) * 100.0 for i in range(len(x_gc)) ]

    cpg = rolling_dinuc_counts(seq, WIN, "CG")
    gpc = rolling_dinuc_counts(seq, WIN, "GC")
    islands = predict_cpg_islands(seq)

    payload = dict(
        seq=seq,
        type=seq_type,
        win=WIN,
        x_gc=x_gc,
        gc_pct=gc_pct,
        cpg_per100=cpg,
        gpc_per100=gpc,
        islands=islands,
    )

    html_path = out_dir / f"{prefix}.gc_content.html"
    make_html(payload, html_path)
    print(f"[OK] HTML report: {html_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
