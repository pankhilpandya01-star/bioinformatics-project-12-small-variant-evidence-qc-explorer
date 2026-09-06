# Data Provenance

## Public sources

| Material | Source | Selected content |
| --- | --- | --- |
| Paired reads | NCBI SRA experiment [`SRX10301019`](https://www.ncbi.nlm.nih.gov/sra/SRX10301019), run `SRR13921545` | Complete run: 5,094,925 pairs |
| Reference assembly | NCBI RefSeq [`GCF_000005845.2`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/) / ASM584v2 | Chromosome `NC_000913.3`, 4,641,652 bases |
| Upstream workflow | [Project 11 — Paired-End Mapping & Fragment QC Explorer](https://github.com/pankhilpandya01-star/bioinformatics-project-11-paired-end-fragment-qc-explorer) | Q20 trimmed, retained paired BAM design |

The SRA record describes paired genomic DNA sequencing from *E. coli* K-12
MG1655 on an Illumina NovaSeq 6000. Only public experiment identifiers and
scientific metadata are retained here.

## Acquisition

Reads were retrieved on 2026-09-05 using SRA Toolkit 3.4.1:

```text
fastq-dump --split-files --skip-technical SRR13921545
```

Each mate contains 5,094,925 records and 509,492,500 bases. The first 5,000 R1
and R2 records reproduce Project 11 exactly. The project tested 5,000, 500,000,
and 1,000,000 pairs before following the predefined fallback to the complete
run; only the complete run passed both coverage requirements.

## Source checksums

| File | SHA-256 |
| --- | --- |
| Full raw R1 | `2ef045ec96d1576fb33a58f5c53cadaf3e76ee7c3f86732b77bc9052b219bdb3` |
| Full raw R2 | `72439b12192327391cfca0f1c775c8b03c1ac8966e76d9c26017cb8c5d0435b2` |
| MG1655 reference FASTA | `53bb6a51b6e92139ced1e38f74b7938781027c52200922ff03718c2237d23bb4` |

The machine-readable copy is [`data/source_metadata.csv`](../data/source_metadata.csv).

## Upstream preprocessing and mapping

The complete run was processed with Project 11's fixed settings:

- Q20 3′ trimming on both mates;
- the documented R1 and R2 Illumina adapters;
- minimum length 75 with Cutadapt `pair-filter=any`;
- complete retained pairs only for the main paired alignment;
- Bowtie2 `--end-to-end --very-sensitive --seed 0 --reorder`;
- FR orientation and a 0–1,000 nt fragment range; and
- one thread and sample/read-group metadata for `SRR13921545`.

The retained cohort contained 5,082,209 pairs and produced 10,164,418 primary
alignment records. The pre-duplicate BAM SHA-256 is
`5a878cbf329b9dec83a326091d38c043ef199d338b7749244b6d2eb05e934077`.

After duplicate marking, the variant-ready BAM SHA-256 is
`4a8174d9cc609191050bac8176e5f67a71ce037e25e33b49bbb7bec458d8951b`.
Its BAI SHA-256 is
`a5f0d2c4528d88a6d206762cf2687749341f5614b7345d7b612b12e8fb6b7b15`.

## Integrity and publication boundary

Checksums were computed before and after each workflow stage and remained
unchanged. BAM outputs passed `samtools quickcheck`; VCF/BCF outputs were indexed
and read back through BCFtools.

Raw and trimmed FASTQs, BAM/BAI files, raw BCF, native indexes, and the complete
position-depth stream remain local. They are reproducible from the public
accessions but are too large for this portfolio repository. Compact normalized
VCFs, evidence tables, summary reports, callable intervals, depth windows,
checksums, and the dashboard are retained.

External sequence data and reference files remain subject to their source
repositories' terms and are not relicensed by this project.
