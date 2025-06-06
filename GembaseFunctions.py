import struct
from Bio.Seq import Seq
import csv
import re
from collections import defaultdict

def br(gi, file_path):
    """
    Reads a 4-byte big-endian integer from a binary file at a position determined by `gi`.

    Parameters:
        gi (int): Index used to calculate byte position (gi * 4).
        file_path (str): Path to the binary file.

    Returns:
        int: The 4-byte integer value read from the file.
    
    Raises:
        ValueError: If fewer than 4 bytes are read (e.g. end of file).
        IOError: If the file cannot be opened or read.
    """
    pos = gi * 4
    with open(file_path, 'rb') as f:
        f.seek(pos)
        buffer = f.read(4)
        if len(buffer) < 4:
            raise ValueError(f"Unexpected end of file or insufficient data at position {pos}.")
        return struct.unpack('>I', buffer)[0]

def extract_genus(species):
    """
    Extracts the genus from a species name.
    """
    return species.split()[0] if species else None

def taxline(species, genera, taxnames, nodes, phylum_info, class_info, order_info, family_info):
    """
    Constructs a taxonomy line from a species name.

    Parameters:
        species (str): Full species name (e.g., "Homo sapiens").
        genera (dict): Map from genus to taxid.
        taxnames (dict): Map from taxid to name.
        nodes (str): Path to the binary file used in br().
        phylum_info, class_info, order_info, family_info (str): Delimited strings (e.g. "|123|456|") containing taxids.

    Returns:
        str: A taxonomy line like "Phylum, Class, Order, Family, Genus, Species".
    """
    genus = extract_genus(species)
    taxid = genera.get(genus, 0)
    taxinf = {}

    while taxid and taxid not in (0, 1):
        taxid = br(taxid, nodes)
        tag = f"|{taxid}|"
        if tag in phylum_info:
            taxinf['phylum'] = taxnames.get(taxid, "")
        elif tag in class_info:
            taxinf['class'] = taxnames.get(taxid, "")
        elif tag in order_info:
            taxinf['order'] = taxnames.get(taxid, "")
        elif tag in family_info:
            taxinf['family'] = taxnames.get(taxid, "")

    return ', '.join([
        taxinf.get('phylum', ''),
        taxinf.get('class', ''),
        taxinf.get('order', ''),
        taxinf.get('family', ''),
        genus or '',
        species or ''
    ])


def seqtrans(dna_seq,code)
    return dna_seq.translate(table=code)

class SpeciesIdentifier:
    def __init__(self, history_file, species_reference_file):
        self.history_file = history_file
        self.species_reference_file = species_reference_file

        self.history = {}  # identifier -> species strain name
        self.name_to_code = {}  # "Genus species" -> AABB.XXX
        self.aabb_to_used_xxx = defaultdict(set)  # AABB -> set of used XXX
        self.code_to_strain_ids = defaultdict(set)  # AABB.XXX -> used YYYY

        self._load_species_reference()
        self._load_history()

    def _load_species_reference(self):
        try:
            with open(self.species_reference_file, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    species, code = row
                    self.name_to_code[species.lower()] = code
                    aabb, xxx = code.split('.')
                    self.aabb_to_used_xxx[aabb].add(int(xxx))
        except FileNotFoundError:
            pass

    def _save_species_reference_entry(self, species_name, aabb_xxx):
        with open(self.species_reference_file, 'a') as f:
            f.write(f"{species_name}\t{aabb_xxx}\n")

    def _load_history(self):
        try:
            with open(self.history_file, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    identifier, species_strain = row
                    self.history[identifier] = species_strain
                    aabb, xxx, yyyy = identifier.split('.')
                    self.code_to_strain_ids[f"{aabb}.{xxx}"].add(int(yyyy))
        except FileNotFoundError:
            pass

    def _save_history_entry(self, identifier, species_strain):
        with open(self.history_file, 'a') as f:
            f.write(f"{identifier}\t{species_strain}\n")

    def _parse_species_from_header(self, header):
        match = re.match(r'^>(\S+)\s+(.+)', header)
        if not match:
            raise ValueError(f"Invalid FASTA header format: {header}")
        return match.group(2).strip()

    def _generate_identifier(self, species_strain):
        tokens = species_strain.split()
        if len(tokens) < 3:
            raise ValueError(f"Species string must include genus, species, and strain: '{species_strain}'")

        genus = tokens[0]
        species = tokens[1]
        strain = ' '.join(tokens[2:])
        species_key = f"{genus} {species}".lower()

        # Check if strain already exists
        for ident, existing in self.history.items():
            if species_strain.lower() == existing.lower():
                raise ValueError(f"Strain already exists in history: {species_strain}")

        aa = genus[:2].upper()
        bb = species[:2].upper()
        aabb = aa + bb

        # Get or assign AABB.XXX
        if species_key in self.name_to_code:
            aabb_xxx = self.name_to_code[species_key]
        else:
            used = self.aabb_to_used_xxx[aabb]
            next_xxx = 1
            while next_xxx in used:
                next_xxx += 1
            xxx = f"{next_xxx:03}"
            aabb_xxx = f"{aabb}.{xxx}"
            used.add(next_xxx)
            self.name_to_code[species_key] = aabb_xxx
            self._save_species_reference_entry(species_key, aabb_xxx)

        # Assign new YYYY for strain
        used_yyyy = self.code_to_strain_ids[aabb_xxx]
        next_yyyy = 1
        while next_yyyy in used_yyyy:
            next_yyyy += 1
        yyyy = f"{next_yyyy:04}"
        used_yyyy.add(next_yyyy)

        identifier = f"{aabb_xxx}.{yyyy}"
        self._save_history_entry(identifier, species_strain)
        return identifier

    def get_identifier_from_fasta_header(self, header):
        species_strain = self._parse_species_from_header(header)
        return self._generate_identifier(species_strain)


