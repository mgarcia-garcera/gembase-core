import argparse
import os
from pathlib import Path
import CoregenomeFunctions as CF
import shutil
import glob
import sys
def main():
    parser = argparse.ArgumentParser(description="Construct a Core-genome from a taxonomic name and a gembases location")
    parser.add_argument("-t", "--taxlevel", required=True, help="Species identifier (e.g., AABB.XXX)")
    parser.add_argument("-g", "--gembases_folder", required=True, help="Path to folder with FASTA files")
    parser.add_argument("-w", "--working_folder", default='./TMP/', help="")
    parser.add_argument("-s", "--software_folder", default="/usr/local/", help="")
    parser.add_argument("-o","--output_folder", default="./OUT/", help="")
    parser.add_argument("-i", "--identity", type=float, default=80, help="")
    parser.add_argument("--synteny", default="NO", help="")
    parser.add_argument("--sum_lim", default=0, help="")
    parser.add_argument("--radius", default=0, help="")
    parser.add_argument("--dist_min", default=0, help="")

    args = parser.parse_args()

    LSTINFO = os.path.join(args.gembases_folder, "LSTINFO")
    PRTFOLDER = os.path.join(args.gembases_folder, "Proteins")
    if not os.path.exists(args.working_folder):
        print(f"creating working directory {args.working_folder}")
        #missing
    
    if not os.path.exists(args.output_folder):
         print(f"creating output directory {args.output_folder}")
         #missing

    #Extract the files that will be used to build the core-genome
    print(f"STEP1: Find all genomes matching the taxonomic level: {args.taxlevel}")
    taxon = args.taxlevel
    Listfiles = CF.find_genomes_by_taxon(args.taxlevel, Path(PRTFOLDER), Path(LSTINFO))

    if not Listfiles:
        raise ValueError("No genome files found matching the specified taxonomic level.")
    ListGenomes = CF.get_basenames(Listfiles)
    #Extract the reference to be used
    reference = CF.find_genome_with_most_genes(Listfiles)
    reference_file = str(f"{PRTFOLDER}/{reference}.prt")
    print(f"STEP2: I will use {reference} as my Reference to build the core genome. The genome is in {reference_file}")
    #define variables for synteny
    if args.synteny.lower() == "yes":
            if not args.sum_lim:
                sum_lim=4
            else:
                sum_lim=args.sum_lim
            if not args.radius:
                radius=5
            else:
                radius=args.radius
            if not args.dist_min:
                dist_min=5
            else:
                dist_min = args.dist_min
    else:
            sum_lim, radius, dist_min=0,0,0
    
    mm = f"{args.software_folder}/lib/BLOSUM60"
    ops = f"{args.software_folder}/bin/opscan"

    if not shutil.which(ops):
         raise ValueError(f"Unable to find {ops}. Not such file or directory")
    
    #build the identity maps
    #output=CF.run_opscan_and_parse(Listfiles, reference, mm, PRTFOLDER, args.output_folder, args.working_folder,sum_lim, radius, dist_min, args.synteny, args.identity)
    output=CF.run_mmseqs_and_parse(Listfiles, reference, PRTFOLDER, args.output_folder, args.working_folder, sum_lim, radius, dist_min, args.synteny, args.identity)
    if output == 0:
         print(f"Orthology pairwise is done")
    else:
         raise ValueError(f"Could not run opscan to build the pairwise similarity files") 
    #build the core-genome file
    output=CF.orthoFinder(args.output_folder, taxon, args.identity,args.synteny)
    if isinstance(output, int):
        print(f"Core genome file for {taxon} writen in {args.output_folder}/CoreGenome-{taxon}.{args.identity}.{args.synteny}.lst with {output} core-genes")
        for file_path in glob.glob(os.path.join(args.output_folder, "*.shrt")):
             os.remove(file_path)
        for file_path in glob.glob(os.path.join(args.output_folder, "*.synt")):
             os.remove(file_path)
    else:
        raise ValueError(f"Unable to build the core-genome file.") 
    

if __name__ == "__main__":
     main()

    