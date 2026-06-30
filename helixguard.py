import sys
from src.sanitize import sanitize

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    df = sanitize(input_file)
    df.to_csv(output_file, sep="\t")
    print(f"Saved sanitized genome to {output_file}")

if __name__ == "__main__":
    main()