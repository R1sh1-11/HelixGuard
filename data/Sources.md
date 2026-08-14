# Project Sources & Datasets

Primary data sources for the HelixGuard genomic privacy tool.

## Participant Genomic Data

Sourced from the Personal Genome Project (PGP). Used for testing and validation of the sanitization pipeline. All figures generated via `scripts/vet_genome.py` and `scripts/carrier_check.py`, verified against a byte-level file comparison.

| Participant | Chip Version | SNPs | Blocklist Hits | Carriers | Collection |
| :--- | :--- | ---: | ---: | ---: | :--- |
| James | v5 | 643,535 | 550 | 5 | [link](https://745d71146d59a622dc9f936edf97db77-99.collections.ac2it.arvadosapi.com/_/) |
| Joshua | v5 | 631,454 | 535 | 7 | [link](https://46de496072e6785dbe2aed5aa92fa119-101.collections.ac2it.arvadosapi.com/_/) |
| Marika | v4 | 601,783 | 112 | 2 | [link](https://5048e1a37aefb72d900d62d091fab991-103.collections.ac2it.arvadosapi.com/_/) |
| Participant 4 | v5 | 638,531 | 525 | 6 | [PGP Collection](https://my.pgp-hms.org/) |
| Participant 5 | v4 | 601,888 | 110 | 5 | [PGP Collection](https://my.pgp-hms.org/) |
| Participant 6 | v5 | 643,535 | 550 | 5 | [PGP Collection](https://my.pgp-hms.org/) |
| Participant 7 | v5 | 638,468 | 503 | 6 | [PGP Collection](https://my.pgp-hms.org/) |
| Participant 8 | v5 | 638,570 | 526 | 4 | [PGP Collection](https://my.pgp-hms.org/) |

**Cohort size: 8 genomes.** Two files originally used in early testing turned out to be exact byte-for-byte duplicates of other files in the set and were removed. The two participants originally numbered 9 and 10 were renumbered to 7 and 8 to fill the resulting gap; the data itself is unchanged, only the file names differ from earlier drafts of this document (missing actual links to sources after Marika).

*Source: [Personal Genome Project](https://pgp.med.harvard.edu/data)*
