from snps import SNPs

def parse_genome(filepath):
    s = SNPs(filepath)
    df = s.snps  # returns a dataframe with rsID as index, columns: chrom, pos, genotype
    df = df[df.index.str.startswith("rs")]  # drop i-prefixed rows
    return df

if __name__ == "__main__":
    df = parse_genome("data/genome_James_Jones_v5_Full_20230726173828.txt")
    print(df.head(20))