# tbDigIn
Digenome data analysis for both Cas9 and Integrase

<!---toc start-->
   * [Local Setup](#local-setup)
      * [Running Tests](#running-tests)
   * [Run the Pipeline](#run-the-pipeline)
      * [Reference Preparation](#reference-preparation)
      * [Execution](#execution)

<!---toc end-->

## Local Setup

- [Install conda][conda-link]

- Install mamba

```console
conda install -c conda-forge mamba
```

- Get a local copy of the `tbDigIn` repo

```console
git clone git@github.com:tomebio/tbDigIn.git
cd tbDigIn
```

- Create the `tbDigIn` conda environment


```console
mamba env create -f environment.yml
```

- Activate the `tbDigIn` conda environment

```console
mamba activate tbDigIn
```


- Install `pytomebio` (developer mode)


```console
python setup.py develop
```

- Build [`fgsv`][fgsv-link]

```console
git clone git@github.com:fulcrumgenomics/fgsv.git
cd fgsv
./mill _.deployLocal
java -jar jars/fgsv.jar -help
cd ..
 ```

The executable JAR for `fgsv` should be in `fgsv/jars/fgsv.jar`.

### Running Tests

To run tests, execute:

```console
bash ci/precommit.sh
```

This will run:

1. Unit test for Python code and Snakemake plumbing (with `pytest`)
2. Linting of Python code (with `flake8`)
3. Code style checking of Python code (with `black`)
4. Type checking of Python code (with `mypy`)
5. Code style checking of Shell code (with `shellcheck`)

## Run the Pipeline

### Reference Preparation

The reference FASTA must be indexed with `bwa` and `samtools faidx`:

```console
bwa index ref.fasta
samtools faidx ref.fasta
```

### Execution

Execute the following command to run the pipeline

```console
bash src/scripts/run_snakemake.sh \
    -s src/snakemake/digenome_seq.smk \
    -c /path/to/config.yml \
    -o /path/to/output
```

An example `config.yml` is shown below:

```yaml
digenome_jar: /path/to/digenome.jar
fgsv_jar: /path/to/fgsv.jar
settings:
  - name: crispr
    fq_dir: /path/to/directory/containing/fastqs
    ref_fasta: /path/to/ref/ref.fasta
    guide: GGGGCCACTAGGGACAGGAT
    enzyme: AAVS1
    pam_three_prime: NGG
    overhang: 0
    max_offset: 2
    samples:
      - AAVS1-RNP-rep1
      - AAVS1-RNP-rep2
  - name: bxb1
    fq_dir: /path/to/directory/containing/fastqs
    ref_fasta: /path/to/ref/ref.fasta
    overhang: 2
    max_offset: 5
    min_forward_reads: 0
    min_reverse_reads: 0
    max_insert_size: 5000
    clipped_start_sequences:
      - GCCGCTAGCGGTGGTTTGTCTGGTCAACCACCGCG
      - CCCGGGATCCCCGGATGATCCTGACGACGGAG
      - CACCACGCGTGGCCGGCTTGTCGACGACGGCG
      - GACCGGTAGCTGGGTTTGTACCGTACACCACTGAG
    samples:
      - HEK12C-Bxb1-attP-rep1
      - HEK12C-mock-rep1
```

The FASTQs for a given group are assumed to be in the given FASTQ directory, and named `<sample-name>_R1_001.fastq.gz` and `<sample-name>_R2_001.fastq.gz`.

The reference for the Integrase samples should be the GRCh38 assembly with the PL312 whole plasmid appended.
The reference for the CRISPR samples can be either the unmodified GRCh38 assembly or the reference used for the Integrase,
but for ease of analysis should be the modified GRCh38 assembly.

The config is organized at two levels: run and group level.
The run level contains configuration that applies to all groups and samples, for example the paths to the `digenomitas` and `fgsv` JARs.
The group level contains configuration that applies to related samples, for example the guide for specific CRISRP samples, or tool-specific parameters recommended for integrase samples.

| Config Key                | Description                                                                                        | Level  | Required | Default |
|---------------------------|----------------------------------------------------------------------------------------------------|--------|----------|---------|
| `digenome_jar`            | The path to the digenomitas JAR                                                                    | Run    | Yes      | NA      |
| `fgsv_jar`                | The path to the fgsv JAR                                                                           | Run    | Yes      | NA      |
| `fq_dir`                  | The directory containing FASTQs with suffixes `_R<#>_001.fastq.gz`                                 | Group  | Yes      | NA      |
| `ref_fasta`               | The path to the reference FASTA, with accompanying BWA index files and FASTA index                 | Group  | Yes      | NA      |
| `name`                    | The name of the sample (FASTQs must be `<name>_R<1 or 2>_001.fastq.gz`                             | Group  | Yes      | NA      |
| `guide`                   | The guide (nucleotide) sequence                                                                    | Group  | No       | None    |
| `enzyme`                  | The name of the enzyme                                                                             | Group  | No       | None    |
| `pam_three_prime`         | The `--pam-three-prime` to `digenom.jar IdentifyCutSites`                                          | Group  | No       | None    |
| `overhang`                | The `--overhang` to `digenom.jar IdentifyCutSites`                                                 | Group  | No       | 0       |
| `max_offset`              | The `--max-offset` to `digenom.jar IdentifyCutSites`                                               | Group  | No       | 2       |
| `min_forward_reads`       | The `--min-forward-reads` to `digenom.jar IdentifyCutSites`                                        | Group  | No       | 4       |
| `min_reverse_reads`       | The `--min-reverse-reads` to `digenom.jar IdentifyCutSites`                                        | Group  | No       | 4       |
| `max_insert_size`         | The `--max-insert-size` to `digenom.jar IdentifyCutSites`                                          | Group  | No       | 1       |
| `clipped_start_sequences` | Allow reads where the 5' clipped bases (in sequencing order) start with one of the given sequences | Group  | No       | None    |

For CRISPR samples, the `guide`, `enzyme`, and `pam_three_prime` should be specified, whereas for Integrase samples these should be left unspecified (omitted).

The `overhang` and `max_offset` should be `0` and `2` respectively for the CRISPR samples, and `2` and `5` respectively for the Integrase samples.

Additionally, the following options should be specified for the Integrase samples:

1. `clipped_start_sequences` should be specified as:

- `GCCGCTAGCGGTGGTTTGTCTGGTCAACCACCGCG` (leading sequence of attR until the overhang)
- `CCCGGGATCCCCGGATGATCCTGACGACGGAG` (reverse complement of attR until the overhang)
- `CACCACGCGTGGCCGGCTTGTCGACGACGGCG` (leading sequence of attL until the overhang)
- `GACCGGTAGCTGGGTTTGTACCGTACACCACTGAG` (reverse complement of attL until the overhang)

2. `min_forward_reads` and `min_reverse_reads` set to zero, so that we do not need support from both strands (a minimum total depth will still be enforced).
3. `max_insert_size` set to `5000`.  This does not filter reads that map beyond the default (`1200`) due to longer than expected mappings in the plasmid sequence.

[conda-link]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
[fgsv-link]: https://github.com/fulcrumgenomics/fgsv