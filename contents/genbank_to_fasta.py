#!/usr/bin/env python

"""
Script to convert .genbank file to .fasta format.
"""

import sys

from Bio import SeqIO


def convert(genbank_in_path, fasta_out_path):
    """Converts a genbank file to a Fasta file.
    """
    with open(genbank_in_path) as in_handle:
        with open(fasta_out_path, 'w') as out_handle:
            SeqIO.write(SeqIO.read(in_handle, 'genbank'), out_handle, 'fasta')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'Usage: genbank_to_fasta.py in.genbank out.fasta'
        sys.exit()
    genbank_in_path = sys.argv[1]
    fasta_out_path = sys.argv[2]
    convert(genbank_in_path, fasta_out_path)
