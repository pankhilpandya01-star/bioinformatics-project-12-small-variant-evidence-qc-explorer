# Quality Review

## Result

Milestone 5 passes. All 27 tests complete successfully in the pinned Python
3.12 environment with Samtools 1.24 and BCFtools 1.24.

The review tests failure paths and biological edge cases as well as the normal
workflow. A failed run must leave neither a final output directory nor a hidden
staging directory.

## Test matrix

| Area | Evidence checked |
| --- | --- |
| Reference and BAM validation | Malformed FASTA, missing BAM index, malformed BAM bytes, coordinate order, exact reference-contig agreement, and read groups |
| Sample validation | Empty or unsafe sample identifiers, a wrong sample, and multiple BAM samples |
| Parameter validation | Invalid ploidy, integer boundaries, Boolean values, zero or negative settings, and non-finite floating-point values |
| Tool validation | Missing executables and exact Samtools 1.24 and BCFtools 1.24 version requirements |
| Duplicate handling | Duplicate records are marked but preserved; their flags exclude them from default depth evidence |
| Read-quality evidence | Samtools is exercised directly to confirm Q20 base-quality and MAPQ20 mapping-quality exclusion |
| Depth logic | Fixed reference denominators, missing or reordered positions, low depth, the 50x minimum high-depth boundary, and excessive depth |
| Evidence filters | Low quality, low depth, high depth, low alternate depth, low alternate fraction, and one-strand-only support are assigned independently |
| Variant forms | SNPs, insertions, deletions, multiallelic splitting, left alignment, and exact duplicate removal |
| VCF validation | Missing headers, multiple samples, missing FORMAT evidence, duplicate variant keys, and exact accepted/rejected partitioning |
| Empty result | A valid header-only BAM produces zero candidates, four readable indexed VCF outputs, summaries, and a dashboard |
| Failure recovery | Forced preparation, calling, and complete-pipeline subprocess failures remove every partial staging tree |
| Published evidence | Controlled truth, authentic VCFs, CSV summaries, dashboard counts, duplicate totals, depth totals, and source checksums reconcile |

## Native-tool checks

The native-tool tests use deliberately constructed records so each exclusion is
observable. At one interval, four reads overlap the same bases: one passes, one
has MAPQ 19, one has base quality 19, and one carries the duplicate flag. With
MAPQ20 and Q20 thresholds, depth is exactly one. With both thresholds lowered to
zero, depth is three because the duplicate remains excluded by Samtools' default
depth policy.

A separate VCF contains a repeat-context deletion twice and one multiallelic
SNP. BCFtools left-aligns the deletion from position 4 to position 1, removes its
exact duplicate, and splits the multiallelic SNP into two biallelic records. The
three resulting keys are distinct.

## Atomicity and immutability

Preparation, calling, and the combined CLI each publish through a temporary
sibling directory. Tests inject failures into all three paths and confirm that
the requested output never appears and that staging directories are removed.
Successful tests also verify that source BAM and reference checksums do not
change.

## Remaining scope boundary

These checks establish deterministic software behavior under the documented
settings. They do not independently resolve the unexpectedly large number of
strong candidates in the authentic MG1655-labelled sample. That remains a
provenance and protocol question, not a software-test failure.
