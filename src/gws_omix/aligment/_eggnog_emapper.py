#!/usr/bin/env python3
"""
_exécuter eggNOG-mapper (mode DIAMOND)_

Usage :
    python _eggnog_emapper.py <fasta> [--itype proteins|transcripts]
                              [--cpus N] [--output-dir DIR]
"""
import os
import sys
import subprocess
from pathlib import Path
from typing import Final

###############################################################################
# CONSTANTES / DÉFAUTS
###############################################################################
DEFAULT_CPUS:   Final[int]  = 25
DEFAULT_ITYPE:  Final[str]  = "proteins"
DEFAULT_DATA:   Final[Path] = Path(
    os.getenv(
        "EGGNOG_DATA_DIR",
        Path(sys.prefix) / "share" / "eggnog-mapper" / "data",
    )
).expanduser()

###############################################################################
# HELPER 1 : s’assurer que la base DIAMOND est là
###############################################################################
def ensure_database(data_dir: Path) -> Path:
    data_dir.mkdir(parents=True, exist_ok=True)
    db_file = data_dir / "eggnog_proteins.dmnd"

    if db_file.exists():
        print(f"[INFO] DIAMOND DB trouvé : {db_file}")
        return db_file

    print(f"[INFO] DIAMOND DB absent → téléchargement…")
    subprocess.run(
        ["download_eggnog_data.py", "-D", "-f", "-y", "--data_dir", str(data_dir)],
        check=True,
    )
    print(f"[INFO] Téléchargement terminé : {db_file}")
    return db_file

###############################################################################
# HELPER 2 : lancer emapper.py
###############################################################################
def run_emapper(
    fasta: Path,
    prefix: str,
    data_dir: Path,
    output_dir: Path,
    itype: str,
    cpus: int,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "emapper.py",
        "-i", str(fasta),
        "-o", prefix,
        "--output_dir", str(output_dir),
        "--data_dir", str(data_dir),
        "--itype", itype,
        "-m", "diamond",
        "--cpu", str(cpus),
        "--override",
    ]
    print(f"[CMD] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    ann = output_dir / f"{prefix}.emapper.annotations"
    if not ann.exists():
        raise FileNotFoundError(f"Résultat attendu manquant : {ann}")
    return ann

###############################################################################
# MAIN
###############################################################################
def main() -> None:
    if len(sys.argv) < 2:
        sys.exit("Usage : _eggnog_emapper.py <fasta> [--itype …] [--cpus …] [--output-dir …]")

    fasta       = Path(sys.argv[1]).resolve()
    itype       = DEFAULT_ITYPE
    cpus        = DEFAULT_CPUS
    output_dir  = Path(f"{fasta.stem}_output").resolve()

    # parsing ultra-léger des options
    if "--itype" in sys.argv:
        itype = sys.argv[sys.argv.index("--itype")+1]
    if "--cpus" in sys.argv:
        cpus = int(sys.argv[sys.argv.index("--cpus")+1])
    if "--output-dir" in sys.argv:
        output_dir = Path(sys.argv[sys.argv.index("--output-dir")+1]).resolve()

    data_dir = DEFAULT_DATA
    ensure_database(data_dir)
    ann = run_emapper(fasta, fasta.stem, data_dir, output_dir, itype, cpus)
    print(f"[✓] Fin – annotations : {ann}")


if __name__ == "__main__":  # pragma: no cover
    main()
