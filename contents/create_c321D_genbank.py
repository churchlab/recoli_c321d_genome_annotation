"""
Script to create the c321D Genbank.
"""

import os
import sys

from Bio import SeqIO

# Append genome-refactor src to the path.
PWD = os.path.dirname(os.path.realpath(__file__ ))
SRC_DIR = os.path.join(PWD, '../../../src')
sys.path.append(SRC_DIR)

from biopython_util import make_gene_misc_feature
from biopython_util import remove_gene_features
from biopython_util import update_feature_note
from genome_region_transplant import excise_genome_region
from genome_region_transplant import overlay_features_in_genome_at_position
from misc.genbank_to_ncbi_five_column import convert_genbank_to_five_column
from refactor_config import GENOMES_DIR
import vcf_to_genbank

# Toggle re-running the whole vcf_to_genbank flow.
# Set to False to make everything run.
GENERATE_TABLE_ONLY = False

REC1_C321D_ROOT = os.path.join(GENOMES_DIR, 'rec1_c321d')
MG1655_GENBANK = os.path.join(GENOMES_DIR, 'mg1655', 'mg1655.genbank')
RECOLI_ALL_SNPS = os.path.join(REC1_C321D_ROOT, 'recoli_all_snps.vcf')
C321D_SNPS_VCF = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_snps.vcf')
REC1_C321D_SAMPLE_ID = 'recoli_misq_c31_321D'
REC1_C321D_OUTPUT_ROOT = os.path.join(REC1_C321D_ROOT, 'rec1_c321d')
REC1_C321D_OUTPUT = os.path.join(REC1_C321D_ROOT, 'rec1_c321d') + '.genbank'
PICKLE_DEST = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_liftover.pickle')
MANUAL_UPDATES = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_manual_fixes.txt')
FINAL_TARGET = '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/rec1_c321d_with_lambda_prophage_annotation.genbank'


if not GENERATE_TABLE_ONLY:
    # HACK: Manually remove some gene feature annotations.
    # Eventually this be handled more elegantly.
    print "Performing preliminary fixes ..."
    with open(MG1655_GENBANK) as fh:
        genome_record = SeqIO.read(fh, 'genbank')
    remove_gene_features(genome_record, {
        'prfA': {'types': ('gene', 'CDS')},
        'bioB': {'types': ('gene', 'CDS')},
        'dacB': {'types': ('mat_peptide',)},

        # SEQ_FEAT.PeptideFeatureLacksCDS
        'mhpC': {'types': ('mat_peptide',)},

        # SEQ_FEAT.PeptideFeatureLacksCDS
        'livJ': {'types': ('mat_peptide',)},
    })
    print "Done."

    # First filter the edits that we actually want to keep.
    # TODO: Make the run() method take as an optional argument a user-specified
    # argument for filtering which vcf rows are incorporated into the Genbank.
    CSV_WITH_POS_TO_KEEP = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_snps.csv')
    vcf_to_genbank.create_filtered_vcf(RECOLI_ALL_SNPS, C321D_SNPS_VCF,
            CSV_WITH_POS_TO_KEEP)

    LOG_FILE = os.path.join(REC1_C321D_ROOT, 'vcf_to_genbank.log')

    # Now run the genbank creator.
    kwargs = {
        'liftover_pickle_dest': PICKLE_DEST,
        'manual_updates_filepath': MANUAL_UPDATES,
        'output_format': 'genbank',
        'verbose': True,
        'log_file': LOG_FILE
    }
    genome_record = vcf_to_genbank.run(genome_record, REC1_C321D_OUTPUT_ROOT,
            C321D_SNPS_VCF, REC1_C321D_SAMPLE_ID, **kwargs)
    print "Done running vcf_to_genbank."


    ### Other manual fixes following this run. This is kind of hacky but doing this
    ### kind of stuff in a script makes it reproducible.

    print "Perform manual fixes ..."

    # Make these genes into misc_features.
    make_gene_misc_feature(genome_record, ['galM'])

    # Update lacZ feature with note about mutation.
    for feature in genome_record.features:
        if (feature.type == 'CDS' and 'gene' in feature.qualifiers and
                feature.qualifiers['gene'][0] == 'lacZ'):
            update_feature_note(feature,
                    "Internal stop codon introduced intentionally for benchmarking recombination frequency by Lajoie et al., 2013")
        elif (feature.type == 'variation' and
                feature.location.start == 2134884):
            # Make it clear that this variation goes with wcaC gene since it
            # falls in overlapping region.
            feature.qualifiers['gene'] = 'wcaC'

    # HACK: Write the results out to match interface below.
    SeqIO.write(genome_record, REC1_C321D_OUTPUT, 'genbank')

    print "Done with manual fixes."


    ### Add the lambda prophage annotations.
    # NOTE: Removed excision code after commit
    # 8fdf9cbd7bdc9de38364e48c9041981942c1e703

    print "Adding lambda phage ..."

    LAMBDA_START = 805241
    LAMBDA_END = 817465
    LAMBDA_PROPHAGE_PATH = os.path.join(REC1_C321D_ROOT,
            'lambda_prophage_recoli.genbank')
    lambda_prophage_seq_record = SeqIO.read(LAMBDA_PROPHAGE_PATH, 'genbank')

    # Manually remove the misc_variant placeholder for lambda.
    prelim_final_target = SeqIO.read(REC1_C321D_OUTPUT, 'genbank')
    remaining_feature_list = []
    for feature in prelim_final_target.features:
        if feature.location.start == LAMBDA_START and feature.type == 'misc_variant':
            continue
        remaining_feature_list.append(feature)
    prelim_final_target.features = remaining_feature_list
    SeqIO.write(prelim_final_target, FINAL_TARGET, 'genbank')

    # Add the lambda positions.
    overlay_features_in_genome_at_position(lambda_prophage_seq_record,
            FINAL_TARGET,
            LAMBDA_START,
            FINAL_TARGET)

    print "Done adding lambda phage ..."

print "Generating .tbl representation required for tbl2asn ..."

# Create the NCBI 5-column format.
NCBI_FIVE_COL_OUTPUT = '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/rec1_c321d_with_lambda_prophage_annotation.tbl'

LOCUS_PREFIX = 'N840'

convert_genbank_to_five_column(FINAL_TARGET, NCBI_FIVE_COL_OUTPUT,
        'CP006698', LOCUS_PREFIX)

print "Done."

print "Running tbl2asn."
import shutil
import subprocess
CURRENT_TBL_FILE = '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/ncbi_submission_prep/2013_08_26.rec1_c321d.tbl'
shutil.copyfile(NCBI_FIVE_COL_OUTPUT, CURRENT_TBL_FILE)
os.chdir('ncbi_submission_prep')
subprocess.call('./run_tbl2asn.sh', shell=True)
