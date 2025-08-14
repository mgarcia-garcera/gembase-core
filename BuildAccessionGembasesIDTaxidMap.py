import csv
from collections import defaultdict

def normalize_prefix(genus, species):
    return genus[:2].upper() + species[:2].upper()

def extract_genus_species(full_name):
    """Extract just the 'Genus species' from full scientific name."""
    parts = full_name.strip().split()
    if len(parts) >= 2:
        return parts[0], parts[1]
    else:
        raise ValueError(f"Invalid scientific name: {full_name}")

def load_taxid_map(taxid_file):
    """Load taxid file: taxid<TAB>Genus species"""
    name_to_taxid = {}
    with open(taxid_file, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 2:
                continue
            elif len(row) >2:
                taxid, name, whateverelse = row
            else:
                taxid, name = row
            name_to_taxid[name.strip()] = taxid.strip()
    return name_to_taxid

def generate_gembases_ids(file1, file2, output_file):
    name_to_taxid = load_taxid_map(file2)

    prefix_to_x = {}
    strain_ids = defaultdict(dict)  # {full_prefix: {strain: YYYY}}
    species_done = {}
    accession_output = []

    for line in open(file1):
        if not line.strip():
            continue
        if len(line.strip().split("\t")) == 2:
            accession, full_name = line.strip().split("\t")
        elif len(line.strip().split("\t")) > 2:
            accession, full_name, whateverelse = line.strip().split("\t")
        else:
            raise ValueError(f"This line is stupid: {line}")
        genus, species = extract_genus_species(full_name)
        strain = full_name.strip()  # treat everything as strain string
        SPname=f"{genus} {species}"
        aabb = normalize_prefix(genus, species)

        if aabb not in prefix_to_x:
            prefix_to_x[aabb] = 1
            species_done[SPname] = prefix_to_x[aabb]
            xxx = species_done[SPname]
        elif SPname not in species_done:
            prefix_to_x[aabb] = prefix_to_x[aabb] +1
            species_done[SPname] = prefix_to_x[aabb]
            xxx = species_done[SPname]
        else:
            xxx = species_done[SPname] 
            
        xxx = f"{xxx:03d}"

        gembase_prefix = f"f{aabb}.{xxx}"

        if strain not in strain_ids[gembase_prefix]:
            strain_ids[gembase_prefix][strain] = len(strain_ids[gembase_prefix]) + 1
        yyyy = f"{strain_ids[gembase_prefix][strain]:04d}"

        gembases_id = f"{gembase_prefix}.{yyyy}"

        taxid = name_to_taxid.get(f"{genus} {species}", "NA")

        accession_output.append((accession, full_name, gembases_id, taxid))

    with open(output_file, "w") as out:
        out.write("Accession\tScientific name\tGembasesID\tTaxID\n")
        for acc, name, gid, tax in accession_output:
            out.write(f"{acc}\t{name}\t{gid}\t{tax}\n")

    print(f"Written output to: {output_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate gembases IDs from accession/scientific name files.")
    parser.add_argument("-a", help="File1: accession<TAB>scientific name (with strain info)")
    parser.add_argument("-t", help="File2: taxid<TAB>scientific name (Genus species only)")
    parser.add_argument("-o", "--output", default="gembases_output.tsv", help="Output TSV file")
    args = parser.parse_args()

    generate_gembases_ids(args.a, args.t, args.output)