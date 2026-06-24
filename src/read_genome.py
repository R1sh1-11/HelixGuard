def read_genome(filepath):
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
            if not rsid.startswith("rs"):
                continue
            print(rsid, genotype)

read_genome() #Insert Filename
