import os
import argparse
import GembaseFunctions as GF
import subprocess
import sys
import shutil
from collections import defaultdict
from Bio import SeqIO
import warnings

warnings.filterwarnings(
    "ignore",
    message="The value of the smallest subnormal for <class 'numpy.float32'> type is zero.",
    module="numpy._core.getlimits"
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
    parser.add_argument("-i", help="Path to the input FASTA file (multi-FASTA allowed, same assembly/species expected)")
    parser.add_argument("-r", "--reference", default="species_reference.tsv", help="Path to reference file (species code)")
    parser.add_argument("-hi", "--history", default="species_history.tsv", help="Path to history file (species+strain identifiers)")
    parser.add_argument("-o", "--output", default="renamed_fastas", help="Output directory for renamed FASTA files")
    parser.add_argument("--db", required=False, help="Path to Bakta database (if not provided, uses $BAKTA_DB)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("-f","--force", action="store_true", help="Force overwrite of output directory if it exists")
    parser.add_argument("-s", "--subdirectories", nargs="*", default=['Replicons','Genes','Proteins','LSTINFO','RNA'], help="List of subdirectories to create. Default is 5 standard directories.")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    parser.add_argument("-e", "--email", required=True, help="Email address for Entrez API")
    
    args = parser.parse_args()
    
    email=args.email
    
    tmpdir = os.path.join(args.output, "TMP")
    os.makedirs(tmpdir, exist_ok=True)
    
    identifier = GF.SpeciesIdentifier(
        reference_file=args.reference,
        history_file=args.history,
        output_dir=tmpdir
    )
    
    id,tmpfasta=identifier.assign_identifier_to_fasta(args.i)
    #With the re-named file, it runs bakta and reannotates.
    print(f"This is my identifier {id} and this is my output fasta {tmpfasta}\n")
    bakta_dir=os.path.join(tmpdir, id)
    bakta_check = GF.check_bakta_output_exists(bakta_dir, id)
    
    if bakta_check == 0:
        print("⚙️ Running Bakta...")
        GF.run_bakta(fasta_file=tmpfasta,output_dir=bakta_dir,prefix=id,threads=args.threads,db_path=args.db,force=args.force)
    else: 
        print(f"✔️ All Bakta output files exist for '{id}'")
    
    #once finished the bakta, we start building the different folders
    GF.ensure_directories(args.output, args.subdirectories)
    ### START OF OUTPUT
    ## 1. Replicon
    old=tmpfasta
    suffix=".fna"
    out_file = f"{id}{suffix}"
    new=os.path.join(bakta_dir, out_file)
    outdir_repl = os.path.join(args.output, "Replicons")

    print(f"this is my directory: {outdir_repl}; and this is my file: {out_file}\n")
    replicon_file = os.path.join(outdir_repl,out_file)
    try:
        GF.reheader_fasta(old, new, replicon_file)
    except Exception as e:
        print(f"❌ ERROR: {e}")
        sys.exit(1)
    ## 2. mapping from fasta
    try: 
        mapping = GF.build_key_mapping(replicon_file)
    except Exception as e:
        sys.exit(1)

    ## 3. TSV
    #input definition for GF.reformat_tsv
    suffix=".tsv"
    tmpsfx=".tmp"
    TSVfile=f"{id}{suffix}"
    TMPfile=f"{id}{tmpsfx}"
    oldTSV=os.path.join(bakta_dir,TSVfile) 
    #output definition for GF.reformat_tsv
    outdir_lstinfo = os.path.join(args.output, "LSTINFO")
    out_tmp = os.path.join(outdir_lstinfo,TMPfile)
    out_tsv = os.path.join(outdir_lstinfo, TSVfile)
    
    try:
        
        GF.reformat_tsv(tsv_file=oldTSV, mapping=mapping, output_file=out_tmp)
        GF.switch_columns_tsv(out_tmp,out_tsv,1,6)
    except Exception as e:
        print(f"❌ ERROR: {e}")
        
    ## 4. GBFF
    suffix=".gbff"
    GFFinput = os.path.join(bakta_dir,f"{id}{suffix}")
    GFFoutput = os.path.join(outdir_lstinfo,f"{id}{suffix}")
    try:
        
        locus_tag_map = GF.build_locus_tag_map(out_tsv)
        GF.replace_gbff_fields(GFFinput, mapping, locus_tag_map, GFFoutput)
    except Exception as e:
        print(f"❌ ERROR: {e}")

    ## 5. FFN & FAA
    ### define inputs
    suffixFFN=".ffn"
    suffixFAA=".faa"
    FFNinput=os.path.join(bakta_dir,f"{id}{suffixFFN}")
    FAAinput=os.path.join(bakta_dir,f"{id}{suffixFAA}")

    
    features, file_prefix = GF.parse_tsv(out_tsv, debug=args.debug)

    print("\n========== HEADER CONSISTENCY CHECK ==========")
    ffn_headers = GF.parse_fasta_headers(FFNinput)
    GF.check_headers(features, ffn_headers, os.path.basename(FFNinput), debug=args.debug)
    faa_headers = GF.parse_fasta_headers(FAAinput)
    GF.check_headers(features, faa_headers, os.path.basename(FAAinput), debug=args.debug)
    print("==============================================\n")

    nuc_dir = os.path.join(args.output, "Genes")
    rna_dir = os.path.join(args.output, "RNA")
    protein_dir = os.path.join(args.output, "Proteins")

    suffix_map = {
        'cds': '.nuc',
        'rrna': '.rrna',
        'trna': '.trna',
        'tmrna': '.tmrna',
        'ncrna': '.ncrna',
        'ncrna-region': '.ncrna'
    }

    out_dirs = {'nuc': nuc_dir, 'rna': rna_dir}
    count_by_type = defaultdict(int)

    print(f"[INFO] Processing FFN file...")
    GF.process_ffn(FFNinput, features, file_prefix, suffix_map, out_dirs, count_by_type, debug=args.debug)

    print(f"[INFO] Processing FAA file...")
    GF.process_faa(FAAinput, features, file_prefix, protein_dir, count_by_type, debug=args.debug)

    print(f"[INFO] Completed. Files written to:")
    print(f"       {nuc_dir}  (nucleotide)")
    print(f"       {rna_dir}  (RNA)")
    print(f"       {protein_dir}  (protein)")

    ## 6. INF file
    TXTsuffix = '.txt'
    TXTinfile = os.path.join(bakta_dir,f"{id}{TXTsuffix}")
    header = next(SeqIO.parse(replicon_file, "fasta")).description
    strain_id, genus, species, strain_info, replicons = GF.parse_fasta_header(header)
    try:
        taxid = GF.get_taxid_from_entrez(genus, species,email)
        taxline = GF.get_taxonomy_lineage(taxid,email)
        INFsuffix='.inf'
        INFfile = os.path.join(outdir_lstinfo, f"{id}{INFsuffix}")
        with open(INFfile, 'w') as out_f:
                GF.process_txt_file(TXTinfile, strain_id, genus, species, strain_info, replicons, taxid, taxline, out_f)
        print(f"✔️ Gembases entry for '{id}' performed successfully")
        
    except Exception as e:
        print(f"❌ ERROR: {e}")
        sys.exit(1)
        
    ## 7. Cleanup

    shutil.rmtree(bakta_dir)
    os.remove(tmpfasta)

if __name__ == "__main__":
    main()