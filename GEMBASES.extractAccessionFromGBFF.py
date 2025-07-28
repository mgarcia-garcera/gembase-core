import PangenomeFunctions as PF
import argparse
import gzip

parser = argparse.ArgumentParser(
        description="""
        Concatenates Fasta files for clustering
        """
    )
#Parses the fasta headers and provides a GEMBASES identifier, according to AABB.XXX.YYYY
parser.add_argument("-i", "--input", required=True, help="Input NxN matrix file")
parser.add_argument("-id", required=True, help="Input NxN matrix file")
args = parser.parse_args()

accessions = PF.extract_accessions_from_gbff(args.input)
#for accession in accessions:
#    print(f"{args.id}\t{accession}")
print(f"{args.id}\t{accessions[0]}")

