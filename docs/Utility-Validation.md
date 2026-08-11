# HelixGuard — Week 7 Utility Validation & LD Documentation

**Lead Analyst:** Kritha

---

## 1. Third-Party Tool Utility Validation Log

| Tool / Library | Result | Reason & Details |
| :--- | :--- | :--- |
| **`snps` Python Library** | **PASS** | File parsed successfully. The library recognized standard GRCh37/GRCh38 chromosome positioning and correctly ingested snp rows after Differential Privacy noise/snapping was applied. |
| **Promethease / 23andMe TSV Format** | **PASS** | File structure preserved all 4 required tab-separated columns (`rsid`, `chromosome`, `position`, `genotype`). Snapped genotypes (`AA`, `AC`, `CC`) matched standard diploid formatting. |
| **Pandas Data Science Pipeline** | **PASS** | Confirmed clean tabular matrix ingestion with `comment='#'` and `low_memory=False` across standard Python data environments. |

---

## 2. Biological Validation Summary (37 LD Neighbors)

The 37 flagged LD neighbors correspond to the local haplotype blocks of our 6 core clinical target genes ($r^2 \ge 0.8$):

1. **`APOE` (Late-Onset Alzheimer's Disease):** `rs405509`, `rs440446` (Chromosome 19 promoter/intronic region variants).
2. **`F5` (Factor V Leiden / Thrombophilia):** `rs6024`, `rs6030` (Chromosome 1 ancestral linkage block).
3. **`BRCA1` (Hereditary Breast/Ovarian Cancer):** `rs169402`, `rs1799966` (Chromosome 17 dense coding region).
4. **`BRCA2` (Hereditary Breast/Ovarian Cancer):** `rs144848`, `rs206118` (Chromosome 13 truncation flankers).
5. **`HTT` (Huntington's Disease):** Chromosome 4 trinucleotide repeat/flanking locus block.
6. **`LDLR` (Familial Hypercholesterolemia):** Chromosome 19 LDL receptor coding region block.