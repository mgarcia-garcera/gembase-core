import os
import sys
import gzip
import argparse
import GembaseFunctionsGBFF as gembases
import re



def extract_species_info(fasta_header):
        """
        Extracts genus, species, and remaining description from a FASTA header.
        Uses regex to find the first occurrence of a species name.
        """
        ignored_words = ['Candidatus', 'MAG', 'Uncultured', 'uncultured', 'environmental', 'bacterium']
        species_regex = re.compile(
            r'(?:(?:' + '|'.join(map(re.escape, ignored_words)) + r')\s+)*'  # optional prefixes
            r'(?P<genus>[A-Z][a-z]+)\s+'               # Genus (e.g., Abditibacterium)
            r'(?P<species>sp\.|[a-z]+)\s+'             # Species (e.g., utsteinense or sp.)
            r'(?:strain\s+|isolate\s+)?'               # optional keywords before strain ID
            r'(?P<strain_id>[A-Z0-9\-_/]+(?:\s+[A-Z0-9\-_/]+)*)',  # Strain ID
            re.IGNORECASE
        )
        # Remove '>' and sequence ID first
    

        match = species_regex.search(fasta_header)
        if match:
            genus = (match.group("genus") or "").strip()
            species = (match.group("species") or "").strip()
            strain = (match.group("strain_id") or 'no strain').strip() 
            return genus, species, strain

        raise ValueError(f"Cannot detect species name in header: {fasta_header}")

def parse_inf_file(inf_path):
    gembases_id = None
    name_line = None

    with open(inf_path, 'r') as f:
        for line in f:
            if line.startswith("GEBMASES_Id:"):
                gembases_id = line.strip().split(":", 1)[1].strip()
            elif line.startswith("NAME:"):
                name_line = line.strip().split(":", 1)[1].strip()
            if gembases_id and name_line:
                break

    if not gembases_id or not name_line:
        raise ValueError(f"Missing required fields in {inf_path}")

    return gembases_id, name_line

def main():
    parser = argparse.ArgumentParser(description="Extract GembasesID and normalized species name from .inf using GembaseFunctions")
    parser.add_argument("-i", required=True, help=".inf file path")
    args = parser.parse_args()
    input_file=args.i
        # Parse .inf file
    gembases_id, name_string = parse_inf_file(input_file)
    # Load GembaseFunctions and instantiate SpeciesIdentifier

    print(f"{gembases_id}\t{name_string}")

if __name__ == "__main__":
    main()