import subprocess
import sys
from pathlib import Path


class EggNOGPipeline:
    def __init__(self,
                 input_fasta: str,
                 output_dir: str,
                 cpus: int = 25,
                 itype: str = 'proteins'):

        self.input_fasta = Path(input_fasta)
        self.output_dir = Path(output_dir)
        self.db_dir = Path(sys.prefix) / 'share' / 'eggnog-mapper' / 'data'
        self.cpus = cpus
        self.itype = itype
        self.prefix = self.input_fasta.stem

        self.dmnd_url = 'http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz'
        self.dmnd_gz_path = self.db_dir / 'eggnog_proteins.dmnd.gz'
        self.dmnd_path = self.db_dir / 'eggnog_proteins.dmnd'

        self.required_db_files = [
            'eggnog.db',
            'eggnog.taxa.db',
            'eggnog.taxa.db.traverse.pkl'
        ]

    def run_cmd(self, cmd):
        print(f"[CMD] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    def ensure_dir(self, path: Path):
        if not path.exists():
            print(f"[INFO] Creating directory: {path}")
            path.mkdir(parents=True, exist_ok=True)

    def ensure_core_db(self):
        missing = [f for f in self.required_db_files if not (self.db_dir / f).exists()]
        if missing:
            print(f"[INFO] Missing core DB files: {missing}")
            self.run_cmd([
                'download_eggnog_data.py',
                '-D', '-f', '-y',
                '--data_dir', str(self.db_dir)
            ])
        else:
            print(f"[INFO] Core DB files already present in {self.db_dir}")

    def ensure_dmnd(self):
        if self.dmnd_path.exists():
            print(f"[INFO] DIAMOND DB already present: {self.dmnd_path}")
            return
        if not self.dmnd_gz_path.exists():
            print("[INFO] Downloading eggnog_proteins.dmnd.gz...")
            self.run_cmd(['wget', self.dmnd_url, '-O', str(self.dmnd_gz_path)])
        print("[INFO] Extracting DIAMOND DB...")
        self.run_cmd(['gunzip', '-f', str(self.dmnd_gz_path)])

    def run_emapper(self) -> Path:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            'emapper.py',
            '-i', str(self.input_fasta),
            '-o', self.prefix,
            '--output_dir', str(self.output_dir),
            '--data_dir', str(self.db_dir),
            '--itype', self.itype,
            '-m', 'diamond',
            '--cpu', str(self.cpus),
            '--override'
        ]
        self.run_cmd(cmd)
        return self.output_dir / f"{self.prefix}.emapper.annotations"

    def clean_output(self, annotation_file: Path):
        print(f"[INFO] Cleaning annotation file: {annotation_file}")
        lines = annotation_file.read_text().splitlines()
        header = next((line for line in lines if line.startswith('#') and 'query' in line), None)
        data = [line for line in lines if not line.startswith('#')]

        if header is None:
            raise ValueError("Missing header line with '#query'.")

        out_file = self.output_dir / 'annotation.tsv'
        out_file.write_text('\n'.join([header] + data) + '\n')
        print(f"[INFO] Saved cleaned TSV: {out_file}")

        annotation_file.unlink()
        print(f"[INFO] Deleted original annotation file: {annotation_file}")

    def run_pipeline(self):
        self.ensure_dir(self.db_dir)
        self.ensure_dir(self.output_dir)

        self.ensure_core_db()
        self.ensure_dmnd()

        annotation = self.run_emapper()
        self.clean_output(annotation)

        print(f"[INFO] Pipeline complete. Final output: {self.output_dir / 'annotation.tsv'}")


# --- Run the pipeline ---
pipeline = EggNOGPipeline(
    input_fasta=sys.argv[1],
    output_dir=sys.argv[2],
    cpus=sys.argv[3],
    itype=sys.argv[4]  # CDS,proteins,genome,metagenome
)
pipeline.run_pipeline()
