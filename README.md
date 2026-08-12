# HelixGuard

A local command-line tool that strips disease-risk markers out of raw 23andMe files before you share them with anyone.

## The problem

If you want a third-party app to build you a diet plan, it asks you to upload your entire 23andMe export. That file has around 640,000 SNPs in it. The app might need fifty of them. The rest sits in someone else's database, including markers for BRCA1/2, Alzheimer's risk, Huntington's, and hereditary clotting disorders. You cannot change your genome after a breach the way you change a password.

This is not hypothetical. openSNP, one of the largest open genomic data platforms, shut down permanently in April 2025 and deleted all user data. The founder cited the 23andMe bankruptcy and a reassessment of what open-access genetic data actually exposes people to.

HelixGuard produces a file you can hand over freely. Carrier status for the target conditions is gone, and it is still a valid 23andMe-format file that downstream tools parse without complaint.

## How it works

1. Parse the raw file (tab-separated, `#` comment header)
2. Look up each rsID against a ClinVar-derived blocklist covering six genes
3. Replace flagged genotypes with the homozygous population-reference genotype for that specific variant
4. Find SNPs in linkage disequilibrium with the flagged ones, since a correlated neighbor can leak what you just redacted
5. Apply Laplace differential privacy noise to those neighbors
6. Write the sanitized file in the same four-column format

Everything runs on your machine. Nothing is uploaded. The only network call in the whole project was a one-time LD reference lookup against a public NIH API, and it never saw a user genotype.

### Target genes

| Gene | Condition | Chr |
| :--- | :--- | :--- |
| BRCA1 | Hereditary breast and ovarian cancer | 17 |
| BRCA2 | Hereditary breast and ovarian cancer | 13 |
| APOE | Alzheimer's disease | 19 |
| HTT | Huntington's disease | 4 |
| F5 | Factor V Leiden | 1 |
| LDLR | Familial hypercholesterolemia | 19 |

## Install

**Python 3.11 is required.** Not 3.12, not 3.14. `diffprivlib` imports private symbols out of `sklearn.tree._tree` that newer scikit-learn releases removed, so anything past 3.11 breaks the differential privacy layer.

```bash
git clone https://github.com/R1sh1-11/HelixGuard.git
cd HelixGuard

python3.11 -m venv venv311
source venv311/bin/activate        # Windows: .\venv311\Scripts\Activate.ps1

pip install -r requirements.txt
```

Verify the DP layer actually imports before you go further:

```bash
python -c "from diffprivlib.mechanisms import Laplace; print('ok')"
```

If that fails you are on the wrong Python. `scikit-learn==1.4.2` is pinned for the same reason, so do not bump it without rerunning that check.

## Usage

```bash
python -m helixguard input.txt output.txt
```

Options:

```bash
python -m helixguard input.txt output.txt --epsilon 0.5   # tighter privacy budget
python -m helixguard input.txt output.txt --no-dp         # blocklist only, no DP layer
```

`--no-dp` is what generates the naive baseline used in the evaluation.

### One gotcha

Anything under `src/` has to be run as a module or the project root will not be on `sys.path`:

```bash
python -m src.redteam --mode compare --ground-truth data/genome_X.txt --save-csv    # works
python src/redteam.py ...                                                          # ModuleNotFoundError
```

Scripts under `scripts/` run normally.

### Building the blocklist from scratch

The generated blocklist database is around 95MB and is gitignored, so you build it locally. Grab `variant_summary.txt.gz` from the ClinVar FTP first.

```bash
python src/build_blocklist.py      # ClinVar dump -> data/blocklist.csv
python src/build_database.py       # csv -> data/blocklist.db
python scripts/import_ld_neighbors.py   # restores the LD table from the tracked csv
```

That last step matters. The LD table came from a 45-minute LDlink API run and lives in `data/ld_neighbors.csv` so it survives database rebuilds. Do not refetch it unless you have to.

## Validation

```bash
python scripts/check_state.py            # health check, 4 verdicts
python scripts/carrier_check.py          # the headline metric
python scripts/batch_validate.py         # all 10 genomes, 5 metrics each
python scripts/epsilon_sweep.py --trials 20
python src/redteam_ld.py --sanitized results/output.txt --ground-truth data/genome_X.txt --label "HelixGuard" --out results/redteam_ld_helixguard.csv
python src/plot_tradeoff.py
python -m pytest tests/ -v               # 17 tests
```

## Results

Evaluated on 10 raw 23andMe exports from the Personal Genome Project (Harvard).

**Carrier concealment: 51 carriers across the cohort, 0 exposed, 100%.**

Every risk genotype present in the ten genomes was replaced with the population reference. Concretely, from one participant:

```
rs429358  ref TT  truth CT  released TT   # APOE4, Alzheimer's risk
rs6025    ref TT  truth CC  released TT   # Factor V Leiden
rs144848  ref AA  truth AC  released AA   # BRCA2
```

Runtime is 9 to 16 seconds per genome. All 10 sanitized files pass third-party parse validation.

### Reading the two metrics

There is a second number in the results, residual disclosure, sitting around 97 to 98 percent. That looks like a catastrophic failure and is not. It counts every targeted SNP where the released value equals the truth, and the overwhelming majority of those are non-carriers whose real genotype already was the reference. Leaving them alone is the correct behavior. Publishing a non-carrier's true non-carrier genotype leaks nothing. Carrier concealment is the number that measures privacy.

### Differential privacy

The DP layer works and moves in the right direction:

| ε | Residual disclosure | Genotypes changed |
| :--- | ---: | ---: |
| no DP | 99.13% | 0 |
| 1.0 | 97.65% | 8.5 |
| 0.5 | 97.47% | 9.55 |
| 0.1 | 97.25% | 10.8 |

20 trials per epsilon. Monotonic, but the effect is small and error bars between adjacent epsilon values overlap. The binding constraint is LD coverage, not the privacy budget.

### The LD attack found nothing

`redteam_ld.py` tries to predict carrier status from the sanitized genotypes of LD neighbors, weighted by r². Against a trivial baseline that always guesses non-carrier:

**LD uplift over baseline: +0.00 points.**

Only 6 of 550 blocklisted targets in a given genome have any LD signal at all, and none of the true carriers was among them. With correct per-variant reference alleles carriers are rare, roughly 1% of blocklisted sites, so guessing non-carrier is already 99% accurate and LD adds nothing on top.

This is a real negative result, not a broken attack. It bounds the threat model with data instead of assuming it.

## Limitations

Worth reading before you trust this with anything real.

- **LD coverage is narrow.** 1,175 pairs but only 15 distinct target SNPs. The LDlink fetch silently skipped most of the rsIDs it queried. The entire LD and DP arm operates on a small slice of the protected surface.
- **Reference substitution does the actual privacy work.** DP over LD neighbors is a secondary layer and at this scale it is marginal.
- **Your chip version decides how much protection you get.** A v4 file in this cohort yielded around 110 blocklist hits. The v5-era files yielded 525 to 550. That is a five-fold gap from chip coverage alone, nothing to do with the person's genetics. HelixGuard can only protect what the array measured.
- **Replacement is coarse.** A site with reference `TT` and truth `GG` gets released as `TT`, a full homozygous flip rather than a subtle edit. Maximum concealment, but a determined observer could notice the transformation.
- **Allele frequencies come from the 10-genome cohort**, not an external population panel. Both `dp_noise.py` and `redteam_ld.py` derive observed alleles this way, and accuracy improves with cohort size.
- **Around 3.5% of rows per file are unrecognized genotypes** that pass through untouched. Almost certainly single-allele calls on X, Y, and mitochondrial DNA. None of the six target genes sit on those chromosomes.
- **Two possible duplicate pairs in the cohort.** participant6/participant7 and participant4/participant8 have identical SNP counts and identical metrics. Verify independence before treating this as ten fully distinct participants.
- **Chip version is unknown for 7 of the 10 genomes** because those exports use a non-standard header that the vetting script cannot read.
- **The naive baseline replaces rather than blanks.** The original design called for blanking blocklisted SNPs. This implementation substitutes reference genotypes instead, which is a stronger baseline than the spec asked for.

## Repo layout

```
helixguard/
├── data/               # genomes, ClinVar dump, LD csv (most of it gitignored)
├── docs/               # utility validation log, tradeoff chart
├── results/            # metrics csvs (tracked), sanitized outputs (not)
├── scripts/            # batch runs, health checks, genome vetting
├── src/                # pipeline, DP, red team, blocklist building
├── tests/              # 17 edge case tests
├── helixguard.py       # CLI entry point
└── requirements.txt
```

`data/*.txt` is in `.gitignore`, so the genome files are only tracked because they were force-added. Anything new you drop in there needs `git add -f` or it will silently vanish.

## Data sources

- **Genomes:** Personal Genome Project (Harvard), `https://pgp.med.harvard.edu/data`. Citations in `data/Sources.md`. The original project plan called for openSNP, which no longer exists.
- **Pathogenic variants:** ClinVar `variant_summary.txt` from NIH.
- **LD reference:** NIH LDlink `ldproxy` endpoint, CEU population, r² ≥ 0.8, backed by 1000 Genomes. Requires a free API token from `https://ldlink.nih.gov/?tab=apiaccess`, set as `LDLINK_TOKEN` in your shell. Never hardcode it.

## Credits

Built for altREU 2026 as a 10-week undergraduate research project.

- Rishi: pipeline, sanitization, differential privacy, red team evaluation, tooling
- Srikritha Kosuri: parsing, ClinVar blocklist construction, LD neighbor logic, utility validation

## License

MIT
