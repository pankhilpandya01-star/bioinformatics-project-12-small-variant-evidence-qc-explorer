# Publication Checklist

## Completed locally

- [x] The controlled 2 kb example calls two SNPs, one insertion, one deletion,
  and rejects the weak constructed error as expected.
- [x] All 27 tests pass in Python 3.12 with Samtools 1.24 and BCFtools 1.24.
- [x] Source FASTA and BAM checksums remain unchanged.
- [x] Duplicate marking preserves every alignment record.
- [x] BAM outputs pass `samtools quickcheck` and have BAI indexes.
- [x] Normalized, evaluated, accepted, and rejected VCFs are compressed,
  indexed, readable, reference-consistent, and free of duplicate keys.
- [x] Every normalized candidate appears exactly once in accepted or rejected
  output.
- [x] VCF, CSV, dashboard, BCFtools, Samtools, and manifest totals agree.
- [x] Empty-callset behavior, native quality filtering, normalization,
  multiallelic splitting, and exact deduplication are tested.
- [x] Forced failures in preparation, calling, and the complete pipeline leave
  no partial output or staging directory.
- [x] Large FASTQ, BAM, BCF, and position-depth artifacts are excluded by
  `.gitignore`.
- [x] Publishable filenames, text, and all eight compressed VCF headers were
  scanned for private paths, credentials, personal metadata, and unwanted
  authorship markers.
- [x] README, methods, limitations, provenance, references, license, result
  reviews, portfolio links, and GitHub Actions workflow are present.

## Requires publication authorization

- [ ] Create the GitHub repository with the intended public name.
- [ ] Review the exact Git file list and commit locally.
- [ ] Push the repository.
- [ ] Confirm the first GitHub Actions run passes on Python 3.12.
- [ ] Verify the README dashboard, badges, and internal links on GitHub.
- [ ] Add or update the portfolio index and prepare the Project 12 post.

No repository creation, commit, push, or external publication is part of the
completed local review.
