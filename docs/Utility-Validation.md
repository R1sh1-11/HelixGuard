# HelixGuard — Week 7 Utility Validation & LD Documentation

**Lead Analyst:** Kritha

---

## 1. Third-Party Tool Utility Validation Log

| Tool / Library | Result | Reason & Details |
| :--- | :--- | :--- |
| **`snps` Python Library** | **PASS** | File parsed successfully. The library recognized standard GRCh37/GRCh38 chromosome positioning and correctly ingested snp rows after Differential Privacy noise/snapping was applied. |
| **Promethease / 23andMe TSV Format** | **PASS** | File structure preserved all 4 required tab-separated columns (`rsid`, `chromosome`, `position`, `genotype`). Snapped genotypes (`AA`, `AC`, `CC`) matched standard diploid formatting. |

---

## 2. Biological Validation Summary (37 LD Neighbors)

The 37 flagged LD neighbors correspond to the local haplotype blocks of our 6 core clinical target genes ($r^2 \ge 0.8$):

1. **`APOE` (Alzheimer's Disease):** `rs405509`, `rs440446` (Chromosome 19 promoter/intronic region variants).
2. **`F5` (Factor V Leiden Clotting):** `rs6024`, `rs6030` (Chromosome 1 ancestral linkage block).
3. **`BRCA1` (Hereditary Breast/Ovarian Cancer):** `rs169402`, `rs1799966` (Chromosome 17 dense coding region).
4. **`BRCA2` (Hereditary Breast/Ovarian Cancer):** `rs144848`, `rs206118` (Chromosome 13 truncation flankers).
5. **`LCT` (Lactase Persistence):** `rs41380347`, `rs41456145` (Chromosome 2 selective sweep haplo-block).
6. **`CFTR` (Cystic Fibrosis):** `rs121908745`, `rs21018` (Chromosome 7 transmembrane regulator region).
