import os
from pathlib import Path
import argparse
import PangenomeFunctions as PF
import sys


def main():
    parser = argparse.ArgumentParser(
        description="""
        Concatenates Fasta files for clustering
        """
    )
    #Parses the fasta headers and provides a GEMBASES identifier, according to AABB.XXX.YYYY
    parser.add_argument("-i", required=True, help="Species identifier")
    parser.add_argument("-f", "--folder", required=True, help="Path to input folder (where args.id should be found")
    parser.add_argument("-o", "--output", help="Output directory for renamed FASTA files")
    parser.add_argument("-r", "--ref", default="gembases.ref", help="species names for AABB.XXX identifiers")
    parser.add_argument("--id", default=0.9, help="Identity value for Vsearch or MMseqs2")
    parser.add_argument("--cov", default=0.9,help="Coverage value for Vsearch or MMseqs2")
    parser.add_argument("--mode", default="mmseqs2", help="Method for clustering")
    
    args = parser.parse_args()
    
    
    species_id=args.i
    species_name = PF.get_value_by_key(args.ref, species_id)
    print(f"I will work with {species_name}")

    input_folder=args.folder
    
    #check existance of input folder
    try:
        os.path.isdir(input_folder)
        print("✅ Input Folder exists")
    except Exception as e:
        print(f"❌ ERROR: {e}")
    pangenome_size = PF.count_files_with_speciesId(input_folder, species_id)
    print(f"{pangenome_size} genomes are found for {species_name}")
    if not args.output:
        cwd = os.getcwd()
        output_file=os.path.join(cwd,f"{species_id}.concat.fasta")
    else:
        try:
            os.path.isdir(args.output)
            print("✅ Output Folder exists")
        except Exception as e:
            print(f"❌ ERROR: {e}")
        output_folder = args.output
        concatenate_file = os.path.join(output_folder, f"{species_id}.concat.fasta")
    
    PF.concatenate_species_fastas(species_id, input_folder, concatenate_file)
    print(f"✅ SUCCESS 1: All genes from {species_name} have been concatenated into {concatenate_file}")
    if args.mode.lower() is "vsearch":
        PF.run_vsearch_clustering(output_file,args.id,str(50))
    elif args.mode.lower() is "mmseqs2":
        PF.run_mmseq2_clustering(species_id,output_file,output_folder,args.id,args.cov,force=T)
    else: 
        print(f"❌ ERROR: Unknown mode: {args.mode}")


if __name__ == "__main__":
    main()

        
