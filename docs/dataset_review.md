# Dataset Review

## Decision

The complete public run `SRR13921545`, containing 5,094,925 paired templates,
is selected for Project 12. It is the first tested cohort that passes both fixed
coverage requirements.

This choice does not depend on finding any variants. A zero-variant authentic
result remains valid if the genome was examined with sufficient read depth.

## Why Project 11's subset is too small

Project 11 used 5,000 pairs to teach mapping. With Project 12's minimum mapping
and base qualities of 20, that cohort provides 0.19x mean depth, covers 14.53%
of the chromosome at least once, and covers no position with ten reads. Most of
the genome therefore cannot be evaluated for small differences.

The gate requires:

- at least 15x mean quality-filtered depth; and
- at least 90% of the 4,641,652 chromosome positions at depth 10 or greater.

## Observed scaling evidence

| Input pairs | Mean depth | Genome covered | Genome at depth 10+ | Result |
| ---: | ---: | ---: | ---: | --- |
| 5,000 | 0.19x | 14.53% | 0.00% | Fail |
| 500,000 | 18.84x | 90.59% | 84.82% | Fail |
| 1,000,000 | 37.69x | 90.75% | 89.93% | Fail |
| 5,094,925 | 191.30x | 91.07% | 90.54% | Pass |

The one-million-pair cohort missed the second threshold by 0.07 percentage
points. The threshold was not weakened; the plan's complete-run fallback was
used instead.

The plateau near 91% covered positions suggests that additional reads cannot
resolve every reference region. Possible causes include genuine sequence
differences, inaccessible or repetitive regions, and limitations of the chosen
end-to-end alignment configuration. The dataset gate measures suitability; it
does not claim which explanation is correct.

## Provenance and integrity

The reads were retrieved from NCBI SRA experiment
[`SRX10301019`](https://www.ncbi.nlm.nih.gov/sra/SRX10301019), run
`SRR13921545`, using SRA Toolkit 3.4.1:

```text
fastq-dump --split-files --skip-technical SRR13921545
```

The full run contains 5,094,925 records and 509,492,500 bases in each mate.
The first 5,000 records reproduce Project 11 exactly:

| Prefix | Project 11 SHA-256 |
| --- | --- |
| R1 | `c3d708ba4b1c7eb4bb95dbae6f7189f291250d129514f4d26b50aad272fecc15` |
| R2 | `a6fb66fd23ecf08c1bf77fda0d1cc075e6a4c1c644c5c0b279e753a85fbf5eb6` |

The reference is the unchanged Project 11 MG1655 RefSeq FASTA,
`GCF_000005845.2`, chromosome `NC_000913.3`. Input and reference checksums were
recomputed after processing and remained unchanged. The generated sorted BAM
was indexed and passed `samtools quickcheck`.

## Preprocessing and mapping evidence

The three Project 11 preprocessing stages retained 5,082,209 complete pairs
(99.75%) and preserved 12,716 rejected pairs separately. Q20 trimming removed
629,103 bases. Adapter trimming matched 39,050 R1 reads and 38,794 R2 reads,
removing 1,576,795 bases.

The retained cohort produced 10,164,418 primary alignment records. Of these,
9,400,402 mapped (92.48%), and 4,591,358 templates were proper pairs (90.34%).

Depth was measured with the same thresholds intended for variant calling:

```text
samtools depth -aa -q 20 -Q 20 trimmed_paired.sorted.bam
```

Here, `-q 20` is the minimum base quality and `-Q 20` is the minimum mapping
quality for `samtools depth`.

## Publication boundary

The raw FASTQs, trimmed FASTQs, BAM, BAI, and position-level depth stream remain
local because they are large reproducible artifacts. This repository retains
portable commands, tool versions, source and generated checksums, compact CSV
summaries, and the exact selection decision.

Duplicate marking is deliberately left for Milestone 2. The same coverage gate
will be checked again after duplicates are marked and excluded from evidence.
