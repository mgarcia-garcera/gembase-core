import struct
from Bio.Seq import Seq
from Bio import Entrez
import csv
import os
import re
import subprocess
import shutil
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


def seqtrans(dna_seq,code):
    return dna_seq.translate(table=code)

class SpeciesIdentifier:
    ignored_words = ['Candidatus', 'MAG', 'Uncultured', 'uncultured', 'environmental', 'bacterium']
    pattern = r'((?:' + '|'.join(map(re.escape, ignored_words)) + r')\s+)*([A-Z][a-z]+)\s+(sp\.|[a-z]+)'
    species_regex = re.compile(pattern, re.IGNORECASE)

    def __init__(self, reference_file, history_file, output_dir):
        self.reference_file = reference_file
        self.history_file = history_file
        self.output_dir = output_dir

        os.makedirs(self.output_dir, exist_ok=True)

        self.reference = self._load_reference()
        self.history = self._load_history()

    def _load_reference(self):
        reference = {}
        if os.path.isfile(self.reference_file):
            with open(self.reference_file, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    if len(row) >= 2:
                        reference[row[1]] = row[0]
        return reference

    def _load_history(self):
        history = {}
        if os.path.isfile(self.history_file):
            with open(self.history_file, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')
                for row in reader:
                    if len(row) >= 2:
                        history[row[1]] = row[0]
        return history

    def _get_aabb(self, genus, species):
        return (genus[:2] + species[:2]).upper()

    def _next_species_code(self, aabb):
        existing_codes = [
            code for code in self.reference.values()
            if code.startswith(aabb)
        ]
        numbers = [int(code.split('.')[1]) for code in existing_codes if '.' in code]
        next_number = max(numbers) + 1 if numbers else 1
        return f"{aabb}.{next_number:03d}"

    def _next_strain_code(self, aabb_xxx):
        existing = [
            full_id for full_id in self.history.values()
            if full_id.startswith(aabb_xxx + ".")
        ]
        numbers = [int(code.split('.')[-1]) for code in existing]
        next_number = max(numbers) + 1 if numbers else 1
        return f"{next_number:04d}"

    def _extract_species_info(self, fasta_header):
        """
        Extracts genus, species, and remaining description from a FASTA header.
        Uses regex to find the first occurrence of a species name.
        """
        # Remove '>' and sequence ID first
        parts = fasta_header[1:].split(maxsplit=1)
        description = parts[1] if len(parts) > 1 else ""

        match = self.species_regex.search(description)
        if match:
            genus = (match.group(2) or "").strip()
            species = (match.group(3) or "").strip()
            strain = description.replace(match.group(0), '').strip() or "no_strain"
            return genus, species, strain

        raise ValueError(f"Cannot detect species name in header: {fasta_header}")

    def assign_identifier_to_fasta(self, fasta_path):
        # Extract species and strain info from the FIRST header
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    genus, species, strain = self._extract_species_info(line)
                    break
            else:
                raise ValueError("No FASTA header found in file.")
    
        species_key = f"{genus} {species}"            # species-level key
        strain_key = f"{species_key} {strain}"        # full strain key
    
        # Step 1: Get or create species-level code (AABB.XXX)
        if species_key in self.reference:
            aabb_xxx = self.reference[species_key]
        else:
            aabb = self._get_aabb(genus, species)
            aabb_xxx = self._next_species_code(aabb)
            self.reference[species_key] = aabb_xxx
            with open(self.reference_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow([aabb_xxx, species_key])
    
        # Step 2: Get or create strain-level code (YYYY)
        if strain_key in self.history:
            identifier = self.history[strain_key]
        else:
            yyyy = self._next_strain_code(aabb_xxx)
            identifier = f"{aabb_xxx}.{yyyy}"
            self.history[strain_key] = identifier
            with open(self.history_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow([identifier, strain_key])
    
        # Generate renamed FASTA
        output_fasta = os.path.join(
            self.output_dir,
            os.path.basename(fasta_path).rsplit('.', 1)[0] + '.renamed.fasta'
        )
    
        with open(fasta_path, 'r') as infile, open(output_fasta, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    parts = line[1:].split(maxsplit=1)
                    original_seq_id = parts[0]
                    rest = parts[1] if len(parts) > 1 else ""
                    new_header = f">{identifier} {rest.strip()} |OLD:{original_seq_id}|"
                    outfile.write(new_header + '\n')
                else:
                    outfile.write(line)
    
        print(f"[INFO] Created: {output_fasta}")
        print(f"[INFO] Identifier for assembly: {identifier} ({strain_key})")
        return identifier

def run_bakta(fasta_file, output_dir, prefix, threads=4, db_path=None, force=False):
    """
    Run Bakta on the provided FASTA file.
    
    :param fasta_file: Path to input FASTA
    :param output_dir: Directory for Bakta output
    :param prefix: Prefix for output files
    :param threads: Number of threads to use
    :param db_path: Path to Bakta database (if None, uses $BAKTA_DB environment variable)
    """
    if not prefix:
        raise ValueError("Prefix is required for Bakta output naming.")
        
    if db_path is None:
        db_path = os.getenv("BAKTA_DB")
        if db_path is None:
            raise ValueError("No Bakta database provided and BAKTA_DB environment variable not set.")
    if os.path.exists(output_dir):
        if force:
            print(f"[INFO] Removing existing output directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            raise FileExistsError(f"Output directory '{output_dir}' already exists. Use --force to overwrite.")
            
    #os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        "bakta",
        "--db", db_path,
        "--output", output_dir,
        "--prefix", prefix,
        "--threads", str(threads),
        fasta_file
    ]

    print(f"[INFO] Running Bakta:\n{' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True)
        print("[INFO] Bakta completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Bakta failed with return code {e.returncode}.")
        raise
