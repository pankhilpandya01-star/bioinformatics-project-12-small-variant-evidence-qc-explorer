# Methods

## Objective and scope

The workflow identifies SNPs and short insertions or deletions in a haploid
*Escherichia coli* sample, then separates candidates that meet fixed read-
evidence thresholds from candidates that do not. An accepted candidate is
supported under this configuration; it is not automatically a verified
biological mutation.

The authentic analysis reuses the paired-end preprocessing and mapping design
from Project 11. Variant calling begins from its coordinate-sorted paired BAM.

## Software environment

The reproducible environment pins:

| Software | Version | Role |
| --- | --- | --- |
| Python | 3.12 | CLI, validation, accounting, reports, and tests |
| Cutadapt | 5.2 | Upstream paired trimming reproduced from Project 11 |
| Bowtie2 | 2.5.5 | Upstream end-to-end paired mapping |
| Samtools | 1.24 | BAM validation, duplicate marking, indexing, and depth |
| BCFtools | 1.24 | Pileup, haploid calling, normalization, indexing, and statistics |
| Matplotlib | 3.11.1 | Four-panel dashboard |
| SRA Toolkit | 3.4.1 | Public-read acquisition |

The default is one thread to make commands deterministic and directly
comparable. Tool versions are checked before a workflow run begins.

## Input validation

The workflow requires a regular reference FASTA, coordinate-sorted BAM, BAM
index, one sample identifier, and a new output path. The FASTA parser rejects
missing sequences, duplicate names, whitespace inside sequences, and unsupported
symbols. `samtools quickcheck` must accept the BAM.

Reference names and lengths must exactly match the BAM `@SQ` records. Every
`@RG` record must include `ID` and `SM`, and all read groups must name the single
requested sample. Sample identifiers must start with a letter or number and may
otherwise contain letters, numbers, dots, underscores, or hyphens.

## Duplicate marking

Repeated fragments remain in the evidence file but receive the SAM duplicate
flag. Samtools is run in the required order:

```text
samtools collate input.bam
samtools fixmate -m name_collated.bam fixmate.bam
samtools sort fixmate.bam
samtools markdup -c -s position_sorted.bam marked.bam
```

The final BAM is reheadered without program-command records that could contain
machine-specific paths, indexed, and checked again. Input and output alignment
record counts must be identical. Duplicate counts are recorded separately.

## Reliable depth and callable positions

Depth is measured at every reference position:

```text
samtools depth -aa -q 20 -Q 20 variant_ready.sorted.bam
```

For `samtools depth`, `-q 20` is the minimum base quality and `-Q 20` is the
minimum mapping quality. Duplicate-flagged records are excluded by Samtools'
default depth policy. Zero-depth positions remain in the denominator because
`-aa` emits the complete reference.

A position is callable when depth is at least 10 and no greater than:

```text
max(50, ceiling(3 × genome-wide median depth))
```

The dataset gate requires mean reliable depth of at least 15 and at least 90%
of the reference at depth 10 or greater. The complete public run is the first
tested cohort that passes both requirements.

## Candidate calling and normalization

BCFtools first computes a quality-filtered pileup with allele, strand, and depth
annotations. The caller uses haploid ploidy because one chromosome sequence is
expected at each bacterial position.

```text
bcftools mpileup --min-MQ 20 --min-BQ 20 --max-depth 1000000 \
  --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP
bcftools call --multiallelic-caller --variants-only --keep-alts --ploidy 1
bcftools norm --fasta-ref reference.fasta --multiallelics -any --rm-dup exact
```

Normalization checks reference alleles, left-aligns indels, splits
multiallelic sites into biallelic rows, and removes exact duplicate records.
Each remaining record receives a stable `chromosome:position:REF:ALT` key.

## Evidence rules

Each normalized candidate is evaluated independently against six rules:

| VCF filter | Rule |
| --- | --- |
| `LowQual` | QUAL is below 30 |
| `LowDepth` | FORMAT/DP is below 10 |
| `HighDepth` | FORMAT/DP exceeds the chromosome-wide high-depth boundary |
| `LowAltDepth` | Alternate FORMAT/AD is below 5 |
| `LowAltFraction` | Alternate AD divided by DP is below 0.80 |
| `OneStrandOnly` | Alternate ADF or ADR is below 1 |

All applicable reasons are retained. A candidate is accepted only when none of
the six warnings applies. Every normalized key must occur exactly once in the
accepted or rejected partition.

## Reporting and publication

The workflow writes a duplicate-marked BAM, depth outputs, raw and normalized
candidate files, evaluated/accepted/rejected VCFs, evidence and summary CSVs,
BCFtools statistics, sanitized logs, checksums, portable commands, and a
four-panel dashboard. Compressed variant files receive CSI indexes.

All work occurs in temporary sibling directories. The final directory appears
only after BAM checks, VCF readability and indexing, candidate partitioning,
artifact accounting, and source-checksum verification succeed. Any failure
removes the temporary tree.
