# Controlled Variant Review

## Decision

Milestone 3 passes. The controlled 2 kb chromosome contains four intended DNA
changes and one deliberately weak error. BCFtools 1.24 found all five
candidates. The four intended changes passed every evidence threshold, while
the weak error remained visible in the rejected output with three explicit
warnings.

## What this example tests

- A SNP changes one reference letter to another.
- An insertion adds letters relative to the reference.
- A deletion removes letters relative to the reference.
- Haploid calling expects one chromosome version at each bacterial position.
- Normalization gives insertions and deletions one consistent representation.
- Filtering labels weak evidence without deleting the candidate.

Each intended variant has ten high-quality alternate-supporting reads: five
forward and five reverse. The weak error has six alternate reads out of ten,
all on the forward strand. This construction tests alternate depth, alternate
fraction, strand balance, and variant quality separately from biological data.

## Expected and observed results

| Candidate | Kind | Reliable depth | Alternate evidence | Result |
| --- | --- | ---: | --- | --- |
| `control:300:T:G` | SNP | 10 | 10 reads, 5 forward + 5 reverse | Accepted |
| `control:700:T:A` | SNP | 10 | 10 reads, 5 forward + 5 reverse | Accepted |
| `control:1000:T:TGG` | Insertion | 10 | 10 reads, 5 forward + 5 reverse | Accepted |
| `control:1300:TAC:T` | Deletion | 10 | 10 reads, 5 forward + 5 reverse | Accepted |
| `control:1600:T:C` | Weak error | 10 | 6 reads, forward only | Rejected |

The weak error receives `LowQual`, `LowAltFraction`, and `OneStrandOnly`. It is
not silently discarded: the same normalized record appears in the evaluated
and rejected VCF files and in `variant_evidence.csv`.

## Workflow

The controlled BAM is processed with the same public thresholds planned for the
authentic sample:

```text
bcftools mpileup --min-MQ 20 --min-BQ 20 \
  --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP
bcftools call --multiallelic-caller --variants-only --keep-alts --ploidy 1
bcftools norm --fasta-ref reference.fasta --multiallelics -any --rm-dup exact
```

The caller emits raw candidates first. Normalization then splits multiallelic
records, standardizes allele representation, and removes exact duplicate rows
before Python applies the six named evidence checks. This order ensures that
every row has one stable variant key before it enters the accepted or rejected
partition.

The implementation follows the official [BCFtools 1.24
manual](https://www.htslib.org/doc/1.24/bcftools.html) for pileup generation,
haploid calling, normalization, compressed VCF output, indexing, and statistics.

## Artifacts and accounting

The controlled result contains:

- raw BCF candidates and a CSI index;
- normalized, evaluated, accepted, and rejected compressed VCFs with CSI indexes;
- one evidence CSV row per normalized candidate;
- variant-type and filter-reason summary CSVs;
- BCFtools statistics, depth summaries, portable commands, and checksums.

The manifest reports five normalized candidates, four accepted candidates, and
one rejected candidate. Tests compare those records with the truth table and
confirm that every normalized key appears exactly once in one partition.

## Interpretation boundary

The controlled result proves that the software recognizes these deliberately
constructed examples under the selected thresholds. It does not establish the
authentic sample's variants. That biological analysis begins in Milestone 4,
and any accepted authentic candidate will still mean “supported under these
settings,” not proven truth.
