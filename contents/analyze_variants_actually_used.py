"""
Script to figure out from the vcf_to_genbank.py logs which variants
were actually acted on in creating the c321D Genbank.
"""

import csv

POSITIONS_ACTIONED_LOG = 'vcf_to_genbank.log'

ORIGINAL_REC1_C321D_SNPS = 'rec1_c321d_snps.csv'

USED_VARIANTS_OUTPUT = 'rec1_c321d_variants_affecting_genbank.csv'


def main():
    # Get the variant positions that we actually used to create the Genbank.
    positions_actioned = set([])
    with open(POSITIONS_ACTIONED_LOG) as position_fh:
        for line in position_fh:
            positions_actioned.add(int(line.strip()))

    # Filter the original snps file and write out the ones that were actually
    # used to create the Genbank as recorded in the log file parsed above.
    with open(USED_VARIANTS_OUTPUT, 'w') as output_fh:
        with open(ORIGINAL_REC1_C321D_SNPS) as original_snps_fh:
            reader = csv.DictReader(original_snps_fh)
            writer = csv.DictWriter(output_fh, reader.fieldnames)
            writer.writeheader()
            for read_row in reader:
                if int(read_row['POS']) in positions_actioned:
                    writer.writerow(read_row)


if __name__ == '__main__':
    main()

