import os
import gzip
import argparse
from Bio import SeqIO

def extract_first_accession(file_path):
    try:
        with gzip.open(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                return record.id
    except Exception as e:
        print(f"❌ Error reading {file_path}: {e}")
    return None
def has_files(folder_path):
    if not os.path.isdir(folder_path):
        return False
    
    return any(
        os.path.isfile(os.path.join(folder_path, f))
        for f in os.listdir(folder_path)
    )


def main():
    parser = argparse.ArgumentParser(description="Extract first accession ID from each .gbff.gz file in a folder.")
    parser.add_argument("-i", help="Path to the folder containing .gbff.gz files")
    parser.add_argument("-o", help="Path to outputfile")
    args = parser.parse_args()

    folder = os.path.abspath(args.i)
    if not has_files(folder):
        print(f"❌ Folder '{folder}' is empty or has already been processed through gembases")
    else:
        print(f"Folder '{folder}' contains files")
        with open(args.o, 'w') as f:
            for file_name in os.listdir(folder):
                if file_name.endswith(".gbff.gz"):
                    file_path = os.path.join(folder, file_name)
                    acc = extract_first_accession(file_path)
                    if acc:
                        f.write(f"{acc}\t{file_path}\n")
        print(f"✅ Successfully stored the accessions for files in Folder '{folder}' into '{args.o}'")

if __name__ == "__main__":
    main()
