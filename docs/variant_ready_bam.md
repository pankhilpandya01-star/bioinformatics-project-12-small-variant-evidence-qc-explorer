# Variant-Ready BAM Review

## Decision

Milestone 2 passes. The full authentic dataset remains suitable for Project 12
after duplicate reads are marked and excluded from independent depth evidence.

This step does not delete repeated records. It sets the SAM duplicate flag so
the evidence remains auditable while downstream tools can ignore repeated
fragments when measuring support.

## Method

Samtools 1.24 requires mate data before duplicate marking. The workflow follows
the documented order:

```text
samtools collate input.bam
samtools fixmate -m name_collated.bam fixmate.bam
samtools sort fixmate.bam
samtools markdup position_sorted.bam marked.bam
```

See the official [Samtools duplicate-marking
workflow](https://www.htslib.org/algorithms/duplicate.html) and
[`samtools markdup` 1.24 documentation](https://www.htslib.org/doc/1.24/samtools-markdup.html).
Optical-duplicate detection was not enabled because the SRA identifiers do not
retain the instrument coordinate fields needed for that classification.

The final header retains the reference and read-group records but omits command
provenance lines that contain native machine paths. The manifest records the
same commands in portable form.

## Duplicate evidence

| Measure | Result |
| --- | ---: |
| Input alignment records | 10,164,418 |
| Output alignment records | 10,164,418 |
| Duplicate records | 1,276,720 |
| Duplicate record fraction | 12.56% |
| Preserved nonduplicate records | 8,887,698 |
| Markdup estimated library size | 16,321,154 |

The exact equality of input and output record counts confirms that marking did
not discard evidence. The duplicate percentage describes alignment records,
not unique biological molecules in the original sample.

## Post-duplicate depth evidence

Depth uses minimum mapping and base qualities of 20. Samtools' default duplicate
filter keeps duplicate-flagged records out of this measurement.

| Measure | Before marking | After duplicate exclusion | Change |
| --- | ---: | ---: | ---: |
| Mean depth | 191.30x | 165.59x | -25.71x |
| Reference positions at depth 10+ | 90.54% | 90.50% | -0.04 percentage points |

The average depth fell by 13.44%, as expected when repeated records stopped
counting as independent evidence. Genome breadth barely changed, showing that
the duplicates mainly added coverage where reads were already present.

The chromosome-wide median depth is 171x. The unusually-high-depth boundary is
therefore 513x, three times the median and above the fixed floor of 50x. A total
of 4,199,245 positions (90.47% of the chromosome) fall between 10x and 513x and
are currently callable. Only 1,640 positions have at least 10x depth but exceed
the high-depth boundary. High depth is a caution flag for ambiguous or repeated
evidence; it is not proof that a position is incorrect.

## Integrity and atomicity

- The final BAM passed `samtools quickcheck` and has a BAI index.
- All 10,164,418 input records appear in the output.
- The input BAM and reference SHA-256 values remained unchanged.
- The final BAM header contains only `@HD`, `@SQ`, and `@RG` records.
- Compact text artifacts passed a publication-safety scan.
- Controlled tests force a subprocess failure and confirm that no partial output
  directory remains.

Large BAM, BAI, and position-level depth files remain local. The project retains
their checksums plus compact duplicate, depth, and command evidence.

## Interpretation boundary

Passing this milestone means the dataset has enough reliable breadth for the
planned candidate-calling step. It does not prove that every callable position
contains the correct reference letter, that every duplicate flag is a true PCR
duplicate, or that a future variant call is biologically correct.
