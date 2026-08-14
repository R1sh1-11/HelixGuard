# HelixGuard

A local command-line tool that strips disease-risk markers out of raw 23andMe files before you share them with anyone.

---

## The problem

If you want a third-party app to build you a diet plan, it asks you to upload your entire 23andMe export. That file has around 640,000 SNPs in it. The app might need fifty of them. The rest sits in someone else's database, including markers for BRCA1/2, Alzheimer's risk, Huntington's, and hereditary clotting disorders. You cannot change your genome after a breach the way you change a password.

This is not hypothetical. openSNP, one of the largest open genomic data platforms, shut down permanently in April 2025 and deleted all user data. The founder cited the 23andMe bankruptcy and a reassessment of what open-access genetic data actually exposes people to.

HelixGuard produces a file you can hand over freely. Carrier status for the target conditions is gone, and it is still a valid 23andMe-format file that downstream tools parse without complaint.

---

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

---

# Quickstart for new users

If you just want to clean your own 23andMe file, this section is all you need.

## 1. Check your Python version

**Python 3.11 is required.** Not 3.12, not 3.13, not 3.14.

`diffprivlib` imports private symbols out of `sklearn.tree._tree` that newer scikit-learn releases removed. On any other Python the differential privacy layer silently fails to import.

```bash
python --version        # must print 3.11.x
```

If it does not, install Python 3.11 before continuing. On Windows, `py -3.11 --version` will tell you whether you already have it.

## 2. Clone and set up

```bash
git clone https://github.com/R1sh1-11/HelixGuard.git
cd HelixGuard
```

**macOS / Linux:**
```bash
python3.11 -m venv venv311
source venv311/bin/activate
pip install -r requirements.txt
```

**Windows PowerShell:**
```powershell
py -3.11 -m venv venv311
.\venv311\Scripts\Activate.ps1
pip install -r requirements.txt
```

Your prompt should now start with `(venv311)`. If it does not, the environment is not active and nothing below will work.

## 3. Confirm the privacy layer loaded

```bash
python -c "from diffprivlib.mechanisms import Laplace; print('ok')"
```

If this prints anything other than `ok`, you are on the wrong Python. Go back to step 1. `scikit-learn==1.4.2` is pinned for the same reason, so do not upgrade it without rerunning this check.

## 4. Build the variant database

The blocklist database is about 95MB, so it is not committed to the repo. You build it once, locally.

Download `variant_summary.txt.gz` from the [ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) and put it in `data/`. It is roughly 440MB compressed.

```bash
python src/build_blocklist.py           # ClinVar dump -> data/blocklist.csv
python src/build_database.py            # csv -> data/blocklist.db
python scripts/import_ld_neighbors.py   # restores the LD table from the tracked csv
```

**Do not skip that last step.** The LD table came from a 45-minute run against the NIH LDlink API. It lives in `data/ld_neighbors.csv` so it survives database rebuilds. Refetching it is slow and requires your own API token.

Confirm everything landed:

```bash
python scripts/check_state.py
```

You want four passing verdicts. If the LD table reports empty, rerun the import step.

## 5. Get your raw data file

From 23andMe: **Settings → Privacy & Data → Access your data → Download raw data.** They email you a link. The file is a `.txt` inside a zip.

The file should look like this, with a `#` comment header followed by four tab-separated columns:

```
# rsid  chromosome  position  genotype
rs4477212   1   82154   AA
rs3094315   1   752566  AG
```

## 6. Sanitize it

```bash
python -m helixguard your_genome.txt cleaned_genome.txt
```

You will see output like:

```
[Parser] Successfully loaded 643535 SNPs.
[Validation] 9600 rows with missing genotype ('--' etc) -- kept as-is, skipped in pipeline.
[Sanitize] 4 SNPs are both blocklisted AND LD neighbors -- blocklist replacement applied.
Replaced 550 flagged SNPs total.
Tagged 25 LD neighbors for DP noise.
[DP] epsilon=1.0 (diffprivlib): perturbed 22, changed 8, skipped 3
Saved sanitized genome to cleaned_genome.txt
```

`cleaned_genome.txt` is the file you share. Keep your original somewhere private.

Runtime is about 11 seconds for a typical file.

## 7. Options

```bash
# tighter privacy budget: more noise on correlated neighbors
python -m helixguard input.txt output.txt --epsilon 0.5

# blocklist substitution only, no differential privacy layer
python -m helixguard input.txt output.txt --no-dp
```

**`--epsilon`** defaults to 1.0. Lower means more noise and less leakage through correlated positions, at the cost of slightly less faithful data. Values of 1.0, 0.5, and 0.1 are what we tested.

**`--no-dp`** turns off the noise layer entirely. This is mainly for generating a comparison baseline, not something a normal user needs.

---

## Understanding your output

**The row count does not change.** HelixGuard replaces genotypes in place, it never deletes rows. Your cleaned file has the same number of lines as the original. That is deliberate, since a file with missing rows is both obviously tampered with and more likely to be rejected by whatever tool you feed it to.

**Very few genotypes actually change.** On a typical v5 file, around 550 positions match the blocklist, but roughly 545 of those already carry the population-reference genotype, because most people are not carriers at most sites. So the real edit count is small. On our reference genome it was 13 total: 5 blocklist replacements plus 8 differential privacy perturbations.

**That is the tool working, not failing.** A non-carrier's true genotype at a risk position reveals nothing sensitive. HelixGuard changes what needs changing and leaves the rest alone.

To see exactly what changed:

```bash
python -c "
import pandas as pd
from src.parser import parse_genome
b = parse_genome('your_genome.txt')
a = pd.read_csv('cleaned_genome.txt', sep='\t', low_memory=False)
m = b.merge(a, on='rsid', suffixes=('_b','_a')).dropna(subset=['genotype_b'])
d = m[(m.genotype_b != m.genotype_a) & (m.genotype_a != '--')]
print(d[['rsid','genotype_b','genotype_a']].to_string())
"
```

---

## For developers

### Running from `src/`

Anything in `src/` must be run as a module or the project root will not be on `sys.path`:

```bash
python -m src.redteam --mode compare --ground-truth data/genome_X.txt --save-csv    # works
python src/redteam.py ...                                                          # ModuleNotFoundError
```

Scripts under `scripts/` run normally with `python scripts/name.py`.

**On Windows PowerShell, do not use `\` for line continuation.** It is a bash convention and PowerShell will throw parser errors. Use a backtick or keep the command on one line.

### Validation and evaluation

```bash
python scripts/check_state.py            # DB and cohort health, 4 verdicts
python scripts/carrier_check.py          # carrier concealment, the headline metric
python scripts/batch_validate.py         # all 8 genomes, 5 metrics each
python scripts/vet_genome.py <file>      # profile a single genome before adding it
python scripts/epsilon_sweep.py --trials 20
python src/plot_tradeoff.py              # regenerates docs/privacy_utility_tradeoff.png
python -m pytest tests/ -v               # 17 tests
```

Red team, two separate things:

```bash
# residual disclosure: a baseline count, NOT an attack
python -m src.redteam --mode compare --ground-truth data/genome_X.txt --save-csv

# the actual LD inference attack
python src/redteam_ld.py --sanitized results/output.txt --ground-truth data/genome_X.txt --label "HelixGuard" --out results/redteam_ld.csv
```

### Repo layout

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

---

## Results

Evaluated on 8 raw 23andMe exports from the Personal Genome Project (Harvard).

### Carrier concealment: 40 carriers, 0 exposed, 100%

Every risk genotype present in the ten genomes was replaced with the population reference.

| Genome | Carriers | Exposed | Concealment |
| :--- | ---: | ---: | ---: |
| James | 5 | 0 | 100% |
| Joshua | 7 | 0 | 100% |
| Marika | 2 | 0 | 100% |
| participant4 | 6 | 0 | 100% |
| participant5 | 5 | 0 | 100% |
| participant6 | 5 | 0 | 100% |
| participant7 | 6 | 0 | 100% |
| participant8 | 4 | 0 | 100% |
| **Cohort** | **40** | **0** | **100%** |

Two files originally used in early testing were found to be exact byte-for-byte duplicates of other files in the set and were removed; the participants originally numbered 9 and 10 were renumbered to 7 and 8 to fill the gap. See `data/Sources.md` for the full note.

Concretely, from one participant:

```
rs429358    truth CT    released TT    # APOE4, Alzheimer's risk
rs6025      truth CC    released TT    # Factor V Leiden
rs144848    truth AC    released AA    # BRCA2
rs405509    truth GG    released TT    # APOE promoter
```

All 8 sanitized files pass third-party parse validation. Mean runtime 11 seconds.

### Reading the two metrics

There is a second number, residual disclosure, sitting at 98.61% with DP and 99.13% without. That looks like catastrophic failure and is not.

It counts every targeted SNP where the released value equals the truth. The overwhelming majority of those are non-carriers whose real genotype already was the reference, and leaving them alone is correct behavior. Publishing a non-carrier's true non-carrier genotype leaks nothing.

Carrier concealment is the number that measures privacy. Residual disclosure was designed before the per-variant reference fix and now measures the wrong thing at these sites.

### Differential privacy

| ε | Residual disclosure | Genotypes changed |
| :--- | ---: | ---: |
| no DP | 99.13% | 0 |
| 1.0 | 97.65% | 8.5 |
| 0.5 | 97.47% | 9.55 |
| 0.1 | 97.25% | 10.8 |

20 trials per epsilon. Monotonic in the right direction, but the effect is small and error bars between adjacent epsilon values overlap. The binding constraint is LD coverage, not the privacy budget.

### The LD inference attack

`redteam_ld.py` predicts carrier status from the sanitized genotypes of LD neighbors, weighted by r², against a baseline that always guesses non-carrier.

| Metric | Value |
| :--- | ---: |
| Targets evaluated | 550 |
| True carriers | 5 |
| Targets with LD signal | 6 |
| Baseline accuracy | 99.09% |
| LD attacker accuracy | 99.09% |
| **Uplift over baseline** | **+0.00 pts** |
| Carrier recall | 20% (1 of 5) |
| Carrier precision | 50% |

**Read this carefully.** Accuracy uplift is genuinely zero: with only 6 of 550 targets carrying any LD signal and carriers making up roughly 1% of blocklisted sites, guessing non-carrier is already 99% accurate and the attacker cannot beat it on accuracy. But the attacker did correctly flag one of the five real carriers. So the honest claim is "LD-based re-identification gained nothing measurable at this coverage level," not "the attack found nothing."

This is a real negative result that bounds the threat model with data instead of assuming it.

---

## Limitations

Read these before trusting HelixGuard with anything that matters.

**LD coverage is narrow.** The LD table has 1,175 pairs but only 15 distinct target SNPs, out of 8,942 blocklist entries. The LDlink fetch silently skipped most of the rsIDs it queried, and we caught it late. The entire LD and differential privacy arm operates on a small slice of the protected surface.

**Reference substitution does the actual privacy work.** Differential privacy over LD neighbors is a secondary layer, and at this scale it is marginal. If you are building something similar, measure which component is carrying the weight before optimizing the sophisticated one.

**Your chip version decides how much protection you get.** A v4 file in this cohort yielded 112 blocklist hits. The v5-era files yielded 525 to 550. That is a five-fold gap from array coverage alone, nothing to do with the person's genetics. HelixGuard can only protect what the chip measured, and users on older tests get proportionally less protection through no fault of the tool.

**Replacement is coarse.** A site with reference `TT` and truth `GG` is released as `TT`, a full homozygous flip rather than a subtle edit. This maximizes concealment but is a detectable, blunt transformation. Someone comparing your file against population allele frequencies could in principle notice.

**Only six genes are covered.** BRCA1, BRCA2, APOE, HTT, F5, LDLR. Every other disease marker in your file passes through untouched. This is not a general-purpose genomic privacy tool.

**Allele frequencies come from the 8-genome cohort**, not an external population panel. Both `dp_noise.py` and `redteam_ld.py` derive observed alleles this way. Accuracy improves with cohort size.

**Around 3.5% of rows per file are unrecognized genotypes** that pass through unchanged. Almost certainly single-allele calls on X, Y, and mitochondrial DNA that the validator does not cover. None of the six target genes sit on those chromosomes, so results should be unaffected, but this was never formally investigated.

**Chip version is unknown for 5 of the 8 genomes** (participant4 through participant8), because those exports use a non-standard header that `vet_genome.py` cannot parse.

**The naive baseline replaces rather than blanks.** The original design called for blanking blocklisted SNPs. This implementation substitutes reference genotypes instead, which is a stronger baseline than specified. Noted so nobody assumes blanking when reading the comparison.

**No formal security review.** This is a ten-week undergraduate research project, not audited software. It reduces a specific, measured category of exposure. It is not a guarantee of genomic privacy, and it should not be the only thing standing between you and a decision about who gets your DNA.

---

## Data sources

- **Genomes:** Personal Genome Project (Harvard), `https://pgp.med.harvard.edu/data`. Per-file citations in `data/Sources.md`. The original project plan called for openSNP, which no longer exists.
- **Pathogenic variants:** ClinVar `variant_summary.txt` from NIH. Filtered to the six target genes, then to variants whose clinical significance contains "pathogenic," "risk factor," or "association," then to single nucleotide variants with a valid reference allele. 8,942 rows after filtering.
- **LD reference:** NIH LDlink `ldproxy` endpoint, CEU population, r² ≥ 0.8, backed by 1000 Genomes phase 3. Requires a free API token from `https://ldlink.nih.gov/?tab=apiaccess`, set as the `LDLINK_TOKEN` environment variable. Never hardcode it.

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'src'`** — you ran a file in `src/` directly. Use `python -m src.modulename` instead.

**`ImportError` mentioning `sklearn.tree._tree`** — wrong Python version. You need 3.11 and `scikit-learn==1.4.2`.

**`FileNotFoundError` on `blocklist.db`** — you skipped the database build. See Quickstart step 4.

**`check_state.py` reports the LD table empty** — run `python scripts/import_ld_neighbors.py`.

**Every genotype in the output looks unchanged** — this is usually correct. See "Understanding your output" above. If you want to confirm the pipeline ran, check that the console printed a nonzero "Replaced N flagged SNPs" line.

**PowerShell parser errors with `--flags`** — you copied a multi-line bash command. Put it on one line.

---

## Credits

Built for altREU as a ten-week undergraduate research project.

- **Rishi:** sanitization pipeline, differential privacy layer, red team evaluation, CLI, tooling
- **Srikritha Kosuri:** file parsing, ClinVar blocklist construction, LD neighbor logic, utility validation, genome sourcing
- **Generative AI:** Claude & Gemini

## License

MIT
