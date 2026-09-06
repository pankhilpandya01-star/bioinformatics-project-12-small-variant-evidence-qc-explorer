# Limitations

## Evidence is not biological proof

`PASS` means that a normalized candidate meets the configured read-evidence
rules. It does not prove that the allele is a true mutation, establish its
origin, or show that it changes phenotype.

The authentic sample is labelled MG1655 and the reference is also MG1655, yet
the analysis finds more than 30,000 strongly supported candidates. That scale
is unexpected for a nominally matching sample and reference. Independent sample
provenance, reference choice, library protocol, and alignment review are needed
before the calls are treated as strain differences.

## Reference and alignment dependence

Results depend on one linear MG1655 reference and Project 11's Bowtie2
end-to-end mapping. Reference bias can reduce evidence for divergent alleles.
Repetitive, missing, duplicated, or structurally different regions may map
poorly or ambiguously. About 9.5% of the chromosome remains below the minimum
reliable depth even with the complete run.

No alternative mapper, local-alignment configuration, de novo assembly, or
pangenome representation is compared.

## Duplicate assumptions

Samtools marks likely duplicate fragments from alignment coordinates. Duplicate
records are preserved for auditing but excluded from independent depth evidence.
Coordinate duplication is not proof of PCR duplication, especially in a small
bacterial genome at high coverage. Optical-duplicate classification is not used
because the public identifiers do not retain the needed instrument coordinates.

## Fixed hard filters

The thresholds are intentionally transparent and educational, not trained or
calibrated for every bacterial sequencing experiment. One read on each strand
is only a minimum strand-support check, not a formal strand-bias model. The
three-times-median high-depth rule is a caution boundary rather than a universal
standard.

Variant QUAL and read counts are produced by one caller. No second caller,
replicate concordance, truth set, or orthogonal laboratory validation is used
for the authentic data.

## Variant and sample scope

The workflow supports one sample and ploidy one or two, with ploidy one used for
this project. It reports SNPs and short insertions/deletions only. It does not
perform:

- structural-variant or copy-number analysis;
- gene annotation or effect prediction;
- antimicrobial-resistance or clinical interpretation;
- recalibration, local realignment, or consensus-genome construction;
- phylogeny, joint calling, or population analysis; or
- rescued-singleton calling from Project 11.

Multiallelic records are split for independent evidence accounting. This is
useful for stable keys but does not model complex haplotypes.

## Reproduction boundary

The public repository contains compact VCF, CSV, report, checksum, and dashboard
artifacts. Raw FASTQ, BAM, raw BCF, and per-position depth files remain local due
to size. Reproduction therefore requires retrieving the public reads and
reference, rerunning Project 11 preprocessing and mapping, and then running this
workflow.
