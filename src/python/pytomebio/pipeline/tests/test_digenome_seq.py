"""Tests for the digenome-seq pipeline"""

from typing import Any
from typing import Dict

from py._path.local import LocalPath as TmpDir
from pytomebio.pipeline.tests.util import run_snakemake
from pathlib import Path


def _touch_path(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w"):
        pass
    return path


def test_digenome_seq(tmpdir: TmpDir) -> None:
    """Basic unit test that runs the snakefile in dry-run mode to ensure it
    parses correctly.
    """

    # Set up reference data
    ref_fasta: Path = _touch_path(Path(tmpdir) / "ref" / "ref.fasta")
    for ext in ["amb", "ann", "bwt", "pac", "sa"]:  # BWA files
        _touch_path(Path(f"{ref_fasta}.{ext}"))

    # Set up FASTQs
    samples = ["foo", "bar", "two"]
    fq_dir: Path = Path(tmpdir) / "fastqs"
    for sample in samples:
        _touch_path(fq_dir / f"{sample}_R1_001.fastq.gz")
        _touch_path(fq_dir / f"{sample}_R2_001.fastq.gz")

    # Set up the digenomitas JAR
    digenome_jar: Path = _touch_path(Path(tmpdir) / "digenome.jar")
    fgsv_jar: Path = _touch_path(Path(tmpdir) / "fgsv.jar")

    rules: Dict[str, int] = {
        "align": 3,
        "all": 1,
        "digenomitas_identify_cut_sites": 3,
        "fastq_to_bam": 3,
        "fastqc": 6,
        "fgsv_aggregatesvpileup": 3,
        "fgsv_svpileup": 3,
        "mark_duplicates": 3,
        "multiqc": 1,
        "picard_collect_multiple_metrics": 3,
        "picard_collect_wgs_metrics": 3,
    }

    config: Dict[str, Any] = {
        "digenome_jar": digenome_jar,
        "fgsv_jar": fgsv_jar,
        "settings": [
            {
                "name": "first",
                "fq_dir": fq_dir,
                "guide": "GATTACA",
                "ref_fasta": ref_fasta,
                "enzyme": "TheEnzyme",
                "pam_three_prime": "NGG",
                "samples": ["foo"],
            },
            {
                "name": "second",
                "fq_dir": fq_dir,
                "ref_fasta": ref_fasta,
                "clipped_start_sequences": ["GATTACA", "ACATTAG"],
                "samples": ["bar", "two"],
            },
        ],
    }

    run_snakemake(
        pipeline="digenome-seq",
        workdir=tmpdir,
        rules=rules,
        config=config,
    )
