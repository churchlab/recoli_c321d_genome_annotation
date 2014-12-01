"""
Script to add the realignment SNP calls as annotations to the genome.
"""

import csv

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature

INPUT_GENOME = 'rec1_c321d.gbk'

OUTPUT_GENOME = 'rec1_c321d_with_realignment_snps.gbk'

SNP_CSV_DATA_FILE = 'realignment_snps.csv'

SNP_FIELD_NAMES = [
    'position',
    'ref',
    'alt',
]

REALIGNED_SNP_TYPE = 'realignedSnp'

def main():
    genome_record = SeqIO.read(INPUT_GENOME, 'genbank')

    with open(SNP_CSV_DATA_FILE) as input_fh:
        reader = csv.DictReader(input_fh, SNP_FIELD_NAMES)
        reader.next() # Ignore header.
        for row in reader:
            feature_ref = row['ref']
            feature_alt = row['alt'].replace('[', '').replace(']', '')
            feature_start = int(row['position']) - 1 # pythonic
            feature_end = feature_start + len(feature_alt)
            feature_location = FeatureLocation(feature_start, feature_end)
            feature = SeqFeature(location=feature_location,
                    type=REALIGNED_SNP_TYPE,
                    strand=1)
            feature.qualifiers['ref'] = row['ref']
            feature.qualifiers['alt'] = row['alt']
            genome_record.features.append(feature)

    with open(OUTPUT_GENOME, 'w') as output_fh:
        SeqIO.write(genome_record, output_fh, 'genbank')


if __name__ == '__main__':
    main()
