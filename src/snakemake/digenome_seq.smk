################################################################################
# Pipeline for Digenome Seq Analysis
################################################################################

from pathlib import Path
from typing import List
from typing import Optional

from pytomebio.pipeline import snakemake_utils
import os
import copy
from attrs import define


################################################################################
# Utility methods and variables
################################################################################

@define
class Sample:
    name: str
    group: str
    fq_dir: Path
    ref_fasta: Path
    guide: Optional[str] = None
    enzyme: Optional[str] = None
    pam_three_prime: Optional[str] = None
    overhang: int = 0 
    max_offset: int = 2
    min_forward_reads: int = 4
    min_reverse_reads: int = 4
    max_insert_size: int = 1200
    clipped_start_sequences: Optional[List[str]] = None

    @property
    def fq1(self) -> Path:
        return self.fq_dir / f"{self.name}_R1_001.fastq.gz"
    
    @property
    def fq2(self) -> Path:
        return self.fq_dir / f"{self.name}_R2_001.fastq.gz"


def str_from_config(key: str) -> str:
    if key not in config:
        raise Exception(f"Required parameter '{key}' not found in config")
    return config[key]


def path_from_config(key: str) -> Path:
    return Path(str_from_config(key))


def maybe_param(name: str, value: Optional[str]) -> str:
    if value is None:
        return ""
    else:
        return f"--{name} {value}"


# Required parameters
digenome_jar: Path = path_from_config("digenome_jar")
fgsv_jar: Path = path_from_config("fgsv_jar")

# Read in the samples from the config
samples = []
for group in config['settings']:
    sample_names = group['samples']
    args_dict = copy.deepcopy(group)
    del args_dict['samples']
    args_dict['group'] = args_dict['name']
    # convert keys that need to be type(Path)
    for key in ['fq_dir', 'ref_fasta']:
        args_dict[key] = Path(args_dict[key])
    for name in sample_names:
        # set sample name 
        args_dict['name'] = name
        sample = Sample(**args_dict)
        samples.append(sample)

sample_dict = {sample.name: sample for sample in samples}
sample_names = [sample.name for sample in samples]
assert len(sample_dict) == len(sample_names), (
    "Sample names are not unique across ALL samples"
)


################################################################################
# Terminal files
################################################################################

# File extensions to generate
extensions = [
    "alignment_summary_metrics.txt",
    "cut_sites.txt",
    "mark_duplicates.txt",
    "svpileup.aggregated.txt",
    "wgs_metrics.txt",
]


all_terminal_files: List[Path] = [
    Path(f"{sample.group}/{sample.name}/{sample.name}.{ext}")
    for sample in samples
    for ext in extensions
]

all_terminal_files.append(Path("sequencing_quality_report.html"))

################################################################################
# Snakemake rules
################################################################################

onerror:
    """Block of code that gets called if the snakemake pipeline exits with an error."""
    snakemake_utils.on_error(snakefile=Path(__file__), config=config, log=Path(log))


rule all:
    input:
        all_terminal_files,


rule fastq_to_bam:
    """Converts the input FASTQs to an unmapped BAM
    
    Runs:
    - fgbio FastqToBam
    """
    input:
        fq1 = lambda wildcards: sample_dict[wildcards.sample].fq1,
        fq2 = lambda wildcards: sample_dict[wildcards.sample].fq2,
    output:
        bam = "{group}/{sample}/{sample}.unmapped.bam"
    params:
        rs1 = "+T",
        rs2 = "+T"
    log: "logs/{group}/{sample}.fastq_to_bam.log"
    resources:
        mem_gb = 2
    shell:
        """
        fgbio -Xmx{resources.mem_gb}g \
          --async-io=true \
          FastqToBam \
          --input {input.fq1} {input.fq2} \
          --read-structures {params.rs1} {params.rs2} \
          --output {output.bam} \
          --sample {wildcards.sample} \
          --library {wildcards.sample} \
        &> {log}
        """


rule align:
    """Aligns the reads in the unmapped BAM to the reference genome and coordinate sorts

    Important: runs with -Y so hard-clipping is not performed on supplementary alignments,
    but instead soft-clipping.  This is needed for when digenomitas IdentifyCutSites is
    run given a set of allowable sequences that can be soft-clipped.

    Runs:
    - samtools fastq
    - bwa mem
    - fgbio ZipperBams
    - samtools sort
    """
    input:
        bam = "{group}/{sample}/{sample}.unmapped.bam",
        ref_fasta = lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
        bwa_files = lambda wildcards: [Path(f"{sample_dict[wildcards.sample].ref_fasta}.{ext}") for ext in ["amb", "ann", "bwt", "pac", "sa"]]
    output:
        bam = "{group}/{sample}/{sample}.mapped.bam",
        csi = "{group}/{sample}/{sample}.mapped.bam.csi"
    log: "logs/{group}/{sample}.align.log"
    threads: 18
    resources:
        mem_gb = 32,  # 8 for bwa, 4 for fgbio, and 16 for samtools sort (8 threads times 2G per thread)
        fgbio_mem_gb = 4,
        sort_mem_gb = 2,
        sort_threads = 8
    shell:
        """
        (samtools fastq -N {input.bam} \
          | bwa mem -Y -K 320000000 -t {threads} -p {input.ref_fasta} - \
          | fgbio -Xmx{resources.fgbio_mem_gb}g --compression 0 ZipperBams --unmapped {input.bam} --ref {input.ref_fasta} \
          | samtools sort -m {resources.sort_mem_gb}G -@ {resources.sort_threads} -o {output.bam} --write-index
        ) &> {log}
        """


# TODO: use --BARCODE_TAG RX if UMIs are present
rule mark_duplicates:
    """Marks PCR duplicates
    
    Runs:
    - picard MarkDuplicates
    """
    input:
        bam = "{group}/{sample}/{sample}.mapped.bam"
    output:
        bam = "{group}/{sample}/{sample}.deduped.bam",
        txt = "{group}/{sample}/{sample}.mark_duplicates.txt",
        tmp_dir = temp(directory("{group}/{sample}/tmp_mark_duplicates_dir"))
    log: "logs/{group}/{sample}.mark_duplicates.log"
    resources:
        mem_gb=12,
        jvm_gb=8
    shell:
        """
        picard \
          -Dsamjdk.use_async_io_read_samtools=true \
          -Duse_async_io_write_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          MarkDuplicates \
          --INPUT {input.bam} \
          --OUTPUT {output.bam} \
          --METRICS_FILE {output.txt} \
          --CREATE_INDEX \
          --TMP_DIR {output.tmp_dir} \
        &> {log}
        """


rule picard_collect_multiple_metrics:
    """Collect basic sequencing quality metrics.
    
    Runs:
    - picard CollectMultipleMetrics
    """
    input:
        bam="{group}/{sample}/{sample}.deduped.bam",
        ref_fasta = lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
    output:
        "{group}/{sample}/{sample}.alignment_summary_metrics.txt"
    params:
        prefix="{group}/{sample}/{sample}"
    log:
        "logs/{group}/{sample}.picard_collect_multiple_metrics.log"
    resources:
        mem_gb=12,
        jvm_gb=8
    shell:
        """
        picard \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          CollectMultipleMetrics \
          --INPUT {input.bam} \
          --OUTPUT {params.prefix} \
          --REFERENCE_SEQUENCE {input.ref_fasta} \
          --PROGRAM null \
          --PROGRAM CollectAlignmentSummaryMetrics \
          --PROGRAM QualityScoreDistribution \
          --PROGRAM MeanQualityByCycle \
          --PROGRAM CollectBaseDistributionByCycle \
          --PROGRAM CollectGcBiasMetrics \
          --PROGRAM CollectSequencingArtifactMetrics \
          --PROGRAM CollectInsertSizeMetrics \
          --PROGRAM CollectQualityYieldMetrics \
          --FILE_EXTENSION .txt \
         &> {log}
        """


rule picard_collect_wgs_metrics:
    """Collects WGS coverage metrics
    
    Runs:
    - picard CollectWgsMetrics
    """
    input:
        bam = "{group}/{sample}/{sample}.deduped.bam",
        ref_fasta = lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
    output:
        txt = "{group}/{sample}/{sample}.wgs_metrics.txt"
    log: "logs/{group}/{sample}.picard_collect_wgs_metrics.log"
    resources:
        mem_gb=12,
        jvm_gb=8
    shell:
        """
        picard \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          CollectWgsMetrics \
          -I {input.bam} \
          -O {output.txt} \
          -R {input.ref_fasta} \
        &> {log}
        """


rule fastqc:
    """Collect pre-alignment sequencing quality control
    
    Note: run as a single thread since FASTQC only parallelizes
    across files, so adding threading doesn't help here.
    
    Runs:
    - fastqc
    """
    input:
        fq = lambda wildcards: sample_dict[wildcards.sample].fq_dir / "{sample}_R{read_num}_001.fastq.gz"
    output:
        html = "{group}/{sample}/{sample}.{read_num}_fastqc.html",
        zip = "{group}/{sample}/{sample}.{read_num}_fastqc.zip"
    params:
        outdir = "{group}/{sample}"
    log:
        "logs/{group}/{sample}.fastqc.R{read_num}.log"
    shell:
        """
        (mkdir -p {params.outdir} &&
          (cat {input.fq} \
            | zcat -fc \
            | fastqc \
              --outdir {params.outdir} \
              stdin:{wildcards.sample}.{wildcards.read_num} \
          )
        ) &> {log}
        """


# TODO: Customize the MultiQC config?
rule multiqc:
    """Aggregates metrics into a HTML (MultiQC) report
    
    Runs:
    - MultiQC
    """
    input:
        fastqc_html = [f"{sample.group}/{sample.name}/{sample.name}.{r}_fastqc.html" for sample in samples for r in [1, 2]],
        fastqc_zip = [f"{sample.group}/{sample.name}/{sample.name}.{r}_fastqc.zip" for sample in samples for r in [1, 2]],
        dupe = [f"{sample.group}/{sample.name}/{sample.name}.mark_duplicates.txt" for sample in samples],
        asm = [f"{sample.group}/{sample.name}/{sample.name}.alignment_summary_metrics.txt" for sample in samples],
        wgs = [f"{sample.group}/{sample.name}/{sample.name}.wgs_metrics.txt" for sample in samples],
    output:
        multiqc_report = "sequencing_quality_report.html",
    params:
        directory=f"{os.getcwd()}"
    log:
        f"logs/multiqc.log"
    shell:
        "multiqc {params.directory} --force --no-ansi --filename {output.multiqc_report} &> {log}"


# TODO:
# - tune other parameters?
rule digenomitas_identify_cut_sites:
    """Identify putative cut sites
    
    Runs:
    - digenomitas IdentifyCutSites
    """
    input:
        bam = "{group}/{sample}/{sample}.deduped.bam",
        ref_fasta = lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
    output:
        txt = "{group}/{sample}/{sample}.cut_sites.txt",
    params:
        jar = digenome_jar,
        guide = lambda wildcards: maybe_param('guide', sample_dict[wildcards.sample].guide),
        enzyme = lambda wildcards: maybe_param('enzyme', sample_dict[wildcards.sample].enzyme),
        pam_three_prime = lambda wildcards: maybe_param('pam-three-prime', sample_dict[wildcards.sample].pam_three_prime),
        overhang = lambda wildcards: sample_dict[wildcards.sample].overhang,
        max_offset = lambda wildcards: sample_dict[wildcards.sample].max_offset,
        min_forward_reads = lambda wildcards: sample_dict[wildcards.sample].min_forward_reads,
        min_reverse_reads = lambda wildcards: sample_dict[wildcards.sample].min_reverse_reads,
        max_insert_size = lambda wildcards: sample_dict[wildcards.sample].max_insert_size,
        clipped_start_sequences = lambda wildcards: (
            maybe_param("clipped-start-sequences",
                None if sample_dict[wildcards.sample].clipped_start_sequences is None else (
                    " ".join(sample_dict[wildcards.sample].clipped_start_sequences)
                )
            )
        ),
        output_all_with_clipped_support = lambda wildcards: str(sample_dict[wildcards.sample].clipped_start_sequences is not None).lower()
    log: "logs/{group}/{sample}.digenomitas_identify_cut_sites.log"
    resources:
        mem_gb=5,
        jvm_gb=4
    shell:
        """
        java \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          -jar {params.jar} \
          IdentifyCutSites \
          --input={input.bam} \
          --ref={input.ref_fasta} \
          --output={output.txt} \
          {params.guide} \
          {params.enzyme} \
          {params.pam_three_prime} \
          --overhang={params.overhang} \
          --max-offset={params.max_offset} \
          --min-forward-reads={params.min_forward_reads} \
          --min-reverse-reads={params.min_reverse_reads} \
          --max-insert-size={params.max_insert_size} \
          output_all_with_clipped_support
          --output-all-with-clipped-support={params.output_all_with_clipped_support} \
          {params.clipped_start_sequences} \
        &> {log}
        """


rule fgsv_svpileup:
    """Collates a pileup of putative structural variant supporting reads.

    Runs:
    - fgsv SvPileup
    """
    input:
        bam = "{group}/{sample}/{sample}.deduped.bam"
    output:
        txt = "{group}/{sample}/{sample}.svpileup.txt",
        bam = "{group}/{sample}/{sample}.svpileup.bam"
    params:
        jar = fgsv_jar,
        prefix = "{group}/{sample}/{sample}.svpileup",
    log: "logs/{group}/{sample}.fgsv_svpileup.log"
    resources:
        mem_gb=5,
        jvm_gb=4
    shell:
        """
        java \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          -jar {params.jar} \
          SvPileup \
          --input={input.bam} \
          --output={params.prefix} \
        &> {log}
        """


rule fgsv_aggregatesvpileup:
    """Merges nearby pileups of reads supporting putative breakpoints.

    Runs:
    - fgsv AggregateSvPileup
    """
    input:
        txt = "{group}/{sample}/{sample}.svpileup.txt"
    output:
        txt = "{group}/{sample}/{sample}.svpileup.aggregated.txt",
    params:
        jar = fgsv_jar
    log: "logs/{group}/{sample}.fgsv_aggregatesvpileup.log"
    resources:
        mem_gb=5,
        jvm_gb=4
    shell:
        """
        java \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          -jar {params.jar} \
          AggregateSvPileup \
          --input={input.txt} \
          --output={output.txt} \
        &> {log}
        """
