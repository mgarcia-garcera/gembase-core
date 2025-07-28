import os
import argparse
import GembaseFunctionsGBFF as GF
import subprocess
import sys
import shutil
from collections import defaultdict
from Bio import SeqIO
import warnings


warnings.filterwarnings(
    "ignore",
    module="numpy"
)

def main():
    parser = argparse.ArgumentParser(
        description="""
        Builds a GEMBASES entry for a defined Genome. 
        It detects the type of replicon, runs bakta and builds the different gembases files according to the following directory strucure:
        - Replicons
        - Proteins
        - Genes
        - RNA
        - LSTINFO
        """
    )
    #Parses the fasta headers and provides a GEMBASES identifier, according to AABB.XXX.YYYY
    parser.add_argument("-i", help="Path to the input file (multi-FASTA allowed, GBFF allowed, same assembly/species expected)")
    parser.add_argument("-r", "--reference", default="species_reference.tsv", help="Path to reference file (species code)")
    parser.add_argument("-hi", "--history", default="species_history.tsv", help="Path to history file (species+strain identifiers)")
    parser.add_argument("-o", "--output", default="renamed_fastas", help="Output directory for renamed FASTA files")
   #parser.add_argument("--db", required=False, help="Path to Bakta database (if not provided, uses $BAKTA_DB)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("-f","--force", action="store_true", help="Force overwrite of output directory if it exists")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    parser.add_argument("-e", "--email", required=True, help="Email address for Entrez API")
    parser.add_argument("-l","--log", required=False, help="LogFile")
    
    args = parser.parse_args()

    subdirectories = ['Replicons','Genes','Proteins','LSTINFO','RNA']

    tmpdir = os.path.join(args.output, "TMP")
    os.makedirs(tmpdir, exist_ok=True)
    
    #create log file if not as argument
    if not args.log:
        logfile = GF.define_log_file(args.i, args.output)
        if not os.path.exists(logfile):
            open(logfile, 'a').close()
    else:
        logfile = args.log
        if not os.path.exists(logfile):
            open(logfile, 'a').close()

    if GF.check_gembases_output_noid(logfile):
        print("✅ Pipeline already completed Previously!")
        sys.exit()
    else:
        print("❌ Pipeline did not complete successfully or log missing.")
    loghandle = open(logfile, 'a')
    sys.stdout = loghandle
    sys.stderr = loghandle
    print(f"my input is {args.i}")
    mode=GF.detect_file_format(args.i)
    print(f"This is my mode: {mode}")
    identifier = GF.GBFFSpeciesIdentifier(args.reference, args.history)
    id = identifier.assign_identifier_from_gbff(args.i)

   # if GF.check_gembases_output(logfile, id):
   #     print("✅ Pipeline already completed Previously!")
   #     sys.exit()
   # else:
   #     print("❌ Pipeline did not complete successfully or log missing.")

    print(f"This is my identifier {id}\n")
    GF.ensure_directories(args.output, subdirectories)
    ## 1. Replicon
    suffix=".fna"
    out_file = f"{id}{suffix}"
    outdir_repl = os.path.join(args.output, "Replicons")
    replicon_file = os.path.join(outdir_repl,out_file)
    if os.path.exists(replicon_file):
        print(f"Replicon file:{replicon_file} already done")
    elif os.path.exists(f"{replicon_file}.gz"):
        print(f"Replicon file {replicon_file} already done AND compressed")
    else:
        try:
            GF.export_renamed_fasta_from_gbff(args.i, id, replicon_file)
        except Exception as e:
            print(f"❌ ERROR: {e}")
            sys.exit(1)

    ## 2.mapping

    try: 
        mapping = GF.build_old_to_new_replicon_mapping(replicon_file)
    except Exception as e:
        print(f"❌ ERROR: {e}")
        sys.exit(1)

    ## 3. TSV
    #input definition for GF.reformat_tsv
    suffix=".tsv"
    tmpsfx=".tmp"
    TSVfile=f"{id}{suffix}"
    #output definition for GF.reformat_tsv
    outdir_lstinfo = os.path.join(args.output, "LSTINFO")
    out_tsv = os.path.join(outdir_lstinfo, TSVfile)
    if os.path.exists(out_tsv):
        print(f"TSV file {out_tsv} exists")
    else:
        print(f"Starting TSV conversion for {id} ")
        try:
            GF.extract_features_from_gbff(args.i,id,out_tsv,mapping)
        except Exception as e:
            print(f"❌ ERROR: {e}")
            sys.exit(1)
    
    ## 4.GBFF 
    suffix=".gbff"
    GBFFoutfile=f"{id}{suffix}"
    outfile=os.path.join(outdir_lstinfo, GBFFoutfile)
    if os.path.exists(outfile):
        print(f"GBFF file {outfile} exists")
    elif os.path.exists(f"{outfile}.gz"):
        print(f"GBFF file {outfile} exists AND is compressed")
    else:
        GF.rename_features_in_gbff(gbff_input=args.i,gbff_output=outfile,tsv_mapping_file=out_tsv,renamed_fasta=replicon_file)
    
    ## 5. fastas
    check=os.path.join(args.output, f"Proteins/{id}.prt")
    if os.path.exists(check):
        print(f"Extraction of DNA and Protein features for {id} is already done")
    else:
        GF.extract_feature_fastas(outfile, id, args.output)

    ## 6. Features for .inf
    suffix=".inf"
    INFoutputfile=f"{id}{suffix}"
    INFfilehandle=os.path.join(outdir_lstinfo,INFoutputfile)
    if os.path.exists(INFfilehandle):
        print(f"Genome metadata for {id} is already processed")
    else:
        INFfeatures = GF.summarize_gbff_assembly(outfile)
        INFfeatures["taxline"]=GF.get_taxonomy_lineage(INFfeatures["taxID"] ,args.email)

        #extract the common identifier 
        descriptor=GF.extract_common_description_prefix_auto(args.i)
        print(f"this is my descriptor: {descriptor}")
        with open(INFfilehandle, 'w') as f:
        #print header
            f.write(f"GEBMASES_Id: {id}\n")
            f.write(f"NAME: {descriptor}\n")
            for key, value in INFfeatures.items():
                f.write(f"{key}: {value}\n")

    print(f"The pipeline was successful") 
    
if __name__ == "__main__":
    main()