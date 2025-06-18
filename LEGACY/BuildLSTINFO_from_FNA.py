import argparse
from Bio import Entrez
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Process fasta header and annotation txt file, append taxID using Entrez.")
    parser.add_argument("-f", "--fasta", required=True, help="Input fasta file path")
    parser.add_argument("-t", "--txt", required=True, help="Input annotation txt file path")
    parser.add_argument("-o", "--output", help="Output file path (optional, prints to stdout if omitted)")
    parser.add_argument("-e", "--email", required=True, help="Email address for Entrez API")
    return parser.parse_args()

def parse_fasta_header(header):
    header = header.lstrip('>')
    parts = header.split()
    
    strain_id = parts[0]
    genus = parts[1]
    species = parts[2]
    strain = parts[3] if len(parts) > 3 else ""
    replicons = re.findall(r'\|([^|]+)\|', header)
    
    return strain_id, genus, species, strain, replicons

def get_taxid_from_entrez(genus, species):
    query = f"{genus} {species}"
    handle = Entrez.esearch(db="taxonomy", term=query)
    record = Entrez.read(handle)
    handle.close()
    id_list = record.get("IdList", [])
    if len(id_list) == 0:
        print(f"Warning: No taxID found for {query}", file=sys.stderr)
        return "NA"
    return id_list[0]

def process_txt_file(txt_file_path, strain_id, genus, species, strain_info, replicons, taxid, taxline, output_handle):
    skip_lines = {
        "Bakta:",
        "Software: v1.11.0",
        "Database: v6.0, light",
        "DOI: 10.1099/mgen.0.000685",
        "URL: github.com/oschwengers/bakta"
    }
    
    with open(txt_file_path) as f:
        lines = f.readlines()
    
    for line in lines:
        if line.strip() in skip_lines:
            continue
        if line.startswith("Sequence(s):"):
            replicon_str = " ".join([f"|{r}|" for r in replicons])
            new_line = f"Sequence(s): {strain_id} - {genus} {species} {strain_info} {replicon_str}\n"
            output_handle.write(new_line)
            
            output_handle.write(f"taxID:{taxid}\n")
            output_handle.write(f"taxLine:{taxline}\n")
        else:
            output_handle.write(line)

def get_taxonomy_lineage(taxid):
    """
    Given a taxID, fetch taxonomy info from NCBI and return a string:
    "TAXline: Kingdom; Phylum; Class; Order; Family; Genus"
    If some ranks are missing, use 'NA' for those.
    """
    
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Error fetching taxonomy for taxID {taxid}: {e}")
        return None
    
    if not records or len(records) == 0:
        print(f"No taxonomy record found for taxID {taxid}")
        return None
    
    record = records[0]
    
    # Ranks to extract in order
    desired_ranks = ["kingdom", "phylum", "class", "order", "family", "genus"]
    
    # Initialize dict for rank: name
    rank_dict = {rank: "NA" for rank in desired_ranks}
    
    # The lineage is a list of dicts with 'Rank' and 'ScientificName'
    lineage = record.get("LineageEx", [])
    
    for taxon in lineage:
        rank = taxon.get("Rank", "").lower()
        if rank in desired_ranks:
            rank_dict[rank] = taxon.get("ScientificName", "NA")
    
    # Sometimes the current taxon itself is at one of the desired ranks
    current_rank = record.get("Rank", "").lower()
    current_name = record.get("ScientificName", "NA")
    if current_rank in desired_ranks and rank_dict[current_rank] == "NA":
        rank_dict[current_rank] = current_name
    
    # Format output string
    lineage_str = "; ".join(rank_dict[rank].capitalize() if rank_dict[rank] != "NA" else "NA" for rank in desired_ranks)
    return lineage_str

def main():
    args = parse_args()
    Entrez.email = args.email
    
    with open(args.fasta) as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()
                break
                
    strain_id, genus, species, strain_info, replicons = parse_fasta_header(header)
    taxid = get_taxid_from_entrez(genus, species)
    taxline = get_taxonomy_lineage(taxid)
    if args.output:
        with open(args.output, 'w') as out_f:
            process_txt_file(args.txt, strain_id, genus, species, strain_info, replicons, taxid, taxline, out_f)
    else:
        process_txt_file(args.txt, strain_id, genus, species, strain_info, replicons, taxid, taxline, sys.stdout)

if __name__ == "__main__":
    main()
