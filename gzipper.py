import gzip
import shutil
import os
import glob
import argparse

def gzip_gbff_files(directory,suffix):
    gbff_files = glob.glob(os.path.join(directory, f"*.{suffix}"))
    
    for file_path in gbff_files:
        gz_path = file_path + ".gz"
        with open(file_path, 'rb') as f_in, gzip.open(gz_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(file_path)  # Remove original file
        print(f"Compressed and removed: {file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        GZIPS GBFF files because they take too much space
        """
    )
    #Parses the fasta headers and provides a GEMBASES identifier, according to AABB.XXX.YYYY
    parser.add_argument("-d", help="Path to the input directory (multi-FASTA allowed, GBFF allowed, same assembly/species expected)")
    parser.add_argument("-s", help="Suffix files to compress")
    args = parser.parse_args()
    gzip_gbff_files(args.d, args.s)  # or gzip_gbff_files("/your/target/directory")