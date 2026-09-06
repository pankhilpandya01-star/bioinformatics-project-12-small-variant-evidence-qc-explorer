# Authentic Small-Variant Analysis

## Result

The full `SRR13921545` dataset produced 33,203 distinct normalized candidates.
Under the fixed Project 12 evidence thresholds, 30,973 candidates (93.28%) were
accepted and 2,230 (6.72%) were rejected.

An accepted candidate means that its read evidence passed these settings. It is
not automatically a verified biological mutation.

## Plain-language interpretation

Most candidate differences are strongly supported substitutions of one DNA
letter for another. The accepted calls contain:

| Variant type | Accepted | Rejected | Total |
| --- | ---: | ---: | ---: |
| SNP | 30,786 | 1,776 | 32,562 |
| Insertion | 94 | 408 | 502 |
| Deletion | 93 | 46 | 139 |
| Complex | 0 | 0 | 0 |

The typical accepted candidate has 158 reliable reads, a 100% alternate-allele
fraction, and variant quality 225.42. Accepted depths range from 10 to 354. The
accepted SNP density is about 0.73% of callable reference positions.

The [NCBI SRA record](https://www.ncbi.nlm.nih.gov/sra/SRX10301019) describes
this experiment as an *E. coli* K-12 MG1655 DNA-seq replicate, and the reference
is also MG1655 chromosome `NC_000913.3`. Finding more than thirty thousand
strongly supported differences between nominally matching sample and reference
is therefore an important QC observation. It warrants independent provenance,
reference, and protocol validation before these candidates are described as
true strain mutations. Project 12 reports the evidence without resolving that
discrepancy.

## Callable genome

After duplicate exclusion and Q20/MAPQ20 filtering:

- mean reliable depth is 165.59x;
- median depth is 171x;
- 90.50% of the chromosome has at least 10 reliable reads;
- 90.47% lies between 10x and the 513x high-depth boundary; and
- 9.50% is below 10x, so those positions are not treated as callable.

Only 1,640 positions exceed the high-depth boundary. No normalized candidate
was rejected for `HighDepth`.

## Why candidates were rejected

Filter counts overlap because one candidate may have more than one warning.

| Evidence warning | Candidates |
| --- | ---: |
| Missing support on one strand | 1,859 |
| Alternate depth below 5 | 1,081 |
| Total reliable depth below 10 | 1,000 |
| Alternate fraction below 80% | 793 |
| Variant quality below 30 | 341 |
| Depth above 513 | 0 |

The rejected VCF preserves every one of these candidates. For example, a weak
candidate can simultaneously have low depth, low alternate fraction, and
one-strand-only support; the filter column retains all applicable reasons.

## Substitution pattern

The 30,786 accepted SNPs include 22,408 transitions (`A>G`, `G>A`, `C>T`, or
`T>C`) and 8,378 transversions, a transition/transversion ratio of about 2.67.
All twelve directional substitutions occur, so the result is not dominated by
one single letter-conversion category. This is a descriptive pattern, not proof
of the variants' origin.

## Integrity and accounting

- The source reference and BAM checksums remained unchanged.
- The input BAM passed `samtools quickcheck` again.
- Normalization processed 32,636 raw records, split 545 multiallelic records,
  realigned 335 records, and removed one exact duplicate record.
- Normalized and evaluated VCFs each contain 33,203 records.
- Accepted and rejected VCFs are disjoint and contain 30,973 and 2,230 records.
- All BCF/VCF files are compressed where appropriate, CSI-indexed, and readable.
- VCF, evidence CSV, summary CSV, dashboard, and manifest totals agree.
- Text artifacts and VCF headers contain no native user paths or credentials.

Large BAM, raw BCF, and position-level depth files remain local. Compact VCFs,
indexes, CSV reports, sanitized logs, callable intervals, 10 kb depth windows,
and the dashboard are retained in the project.

## Scope boundary

No gene annotation, effect prediction, resistance claim, clinical
interpretation, structural-variant analysis, or phylogeny is performed. Those
questions require separate methods and additional validation.
