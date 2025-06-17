import os
import Bio import SeqIO




def generate_replicon_file (fna_file, strain_id, assembly_type, output_file):
    """
    Generates a replicon file with reannotated headers
    
    Args:
        fna_file (str): Path to the .fna file from Bakta
        strain_id (str): Strain identifier in the format AABB.XXX.YYYY.
        Assembly_type (str): type of assembly: "Complete", "draft", or "MAG"
        output_file (str): output file to save the new file.

    """

    chromosomes = [ ]
    plasmids = [ ]
    contigs = [ ]

    #load sequences
    records = list(SeqIO.parse(fna_file, "fasta"))

    