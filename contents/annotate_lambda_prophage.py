"""
Script to transfer NCBI annotations to our prophage.

Our lambda prophage cassette is different than the canonical one available in
NCBI so we need to be clever about transferring annotations.

Notes:
Start sequence in rEcoli:
    GTTGGCATTATAAAAAAGCATTGCTTATCAATTTGTTGCAACG

Corresponds to start position in NCBI 27739.
"""

import os
import sys

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

# Append genome-refactor src to the path so we can import
# genome-refactor modules.
PWD = os.path.dirname(os.path.realpath(__file__ ))
SRC_DIR = os.path.join(PWD, '../../../src')
sys.path.append(SRC_DIR)

from biopython_util import maybe_get_feature_gene


RAW_RECOLI_LAMBDA_PROPHAGE = 'lambda_prophage_raw.genbank'

NCBI_LAMBDA_PROPHAGE = 'lambda_prophage_ncbi.genbank'

FINAL_RECOLI_LAMBDA_PROPHAGE = 'lambda_prophage_recoli.genbank'

NCBI_POSITION_OF_RECOLI_START_PYTHONIC = 27738
NCBI_POSITION_OF_RECOLI_END_PYTHONIC = 38002 # NOTE: approximate

RELEVANT_NCBI_FEATURE_TYPES = set(['CDS', 'gene'])

BLA_QUALIFIERS = {
    'gene': ['bla'],
    'locus_tag': ['bla_1']
}


TETR_QUALIFIERS = {
    'gene': ['tetR'],
    'locus_tag': ['tetR_1']
}


def create_gene_feature(gene_name, feature_location, feature_qualifiers):
    """Creates a minimal SeqFeature to represent a gene.
    """
    gene_feature = SeqFeature(feature_location, type='gene')
    gene_feature.qualifiers = {'gene': [gene_name]}
    gene_feature.qualifiers = dict(gene_feature.qualifiers.items() +
            feature_qualifiers.items())
    return gene_feature


def main():
    recoli_phage_record = SeqIO.read(RAW_RECOLI_LAMBDA_PROPHAGE, 'genbank')
    ncbi_phage_record = SeqIO.read(NCBI_LAMBDA_PROPHAGE, 'genbank')

    # Filter the NCBI features to the region of lambda phage that we actually
    # have, and only relevant features.
    ncbi_phage_relevant_features = []
    for ncbi_feature in ncbi_phage_record.features:
        if (ncbi_feature.type in RELEVANT_NCBI_FEATURE_TYPES and
                ncbi_feature.location.start >=
                        NCBI_POSITION_OF_RECOLI_START_PYTHONIC and
                ncbi_feature.location.end <=
                        NCBI_POSITION_OF_RECOLI_END_PYTHONIC):

            # Adjust the location and append.
            ncbi_phage_relevant_features.append(ncbi_feature._shift(
                    -1 * NCBI_POSITION_OF_RECOLI_START_PYTHONIC))

    # For the NCBI lambda prophage record, create a map from gene name
    # to features with that gene.
    ncbi_gene_to_feature_map = {}
    for feature in ncbi_phage_relevant_features:
        maybe_gene = maybe_get_feature_gene(feature)
        if maybe_gene:
            current_features = ncbi_gene_to_feature_map.get(maybe_gene, [])
            current_features.append(feature)
            ncbi_gene_to_feature_map[maybe_gene] = current_features

    # Now we find the features in the NCBI Genbank record that correspond to
    # those in our raw record.
    updated_features = []
    for feature in recoli_phage_record.features:
        # Try to find a feature with the same length in the NCBI prophage.
        # Copy the properties from the corresponding CDS and add the
        # corresponding gene feature.
        replaced = False
        if feature.type == 'CDS':
            for ncbi_feature in ncbi_phage_relevant_features:
                if ncbi_feature.type == 'gene':
                    continue
                if (len(ncbi_feature) == len(feature) and
                        ncbi_feature.location.start == feature.location.start):
                    replaced = True

                    # Get the corresponding gene feature, if possible.
                    maybe_gene = maybe_get_feature_gene(ncbi_feature)
                    added_gene_feature = False
                    if maybe_gene:
                        gene_features = ncbi_gene_to_feature_map[maybe_gene]
                        for gene_feature in gene_features:
                            if gene_feature.type == 'gene':
                                updated_features.append(gene_feature)
                                added_gene_feature = True
                                break

                    assert added_gene_feature, (
                            "No corresponding gene feature found. Handling " +
                            "case has not been implemented.")

                    # Replace the raw feature with the NCBI feature.
                    updated_features.append(ncbi_feature)

                    # Done replacing this feature.
                    break

        # If not replaced, keep the feature.
        if not replaced:
            updated_features.append(feature)
    recoli_phage_record.features = updated_features

    # Manually update tetR and blah.
    updated_features = []
    for feature in recoli_phage_record.features:
        maybe_gene = maybe_get_feature_gene(feature)
        if maybe_gene and maybe_gene == 'bla':
            updated_features.append(create_gene_feature(
                'bla', feature.location, BLA_QUALIFIERS))
        elif maybe_gene and maybe_gene == 'tetR':
            updated_features.append(create_gene_feature(
                'tetR', feature.location, TETR_QUALIFIERS))
        updated_features.append(feature)
    recoli_phage_record.features = updated_features

    # Validation.

    ### Make sure that each CDS has a corresponding gene feature with a matching
    ### locus tag.

    # Build a map from locus tag to list of features with that tag.
    locus_tag_to_feature_list = {}
    for feature in recoli_phage_record.features:
        if not feature.type in RELEVANT_NCBI_FEATURE_TYPES:
            continue
        locus_tag = feature.qualifiers['locus_tag'][0]
        if locus_tag in locus_tag_to_feature_list:
            locus_tag_to_feature_list[locus_tag].append(feature)
        else:
            locus_tag_to_feature_list[locus_tag] = [feature]

    # Now make sure there are at least two features per locus tag and at least
    # one of them is a gene.
    for locus_tag, feature_list in locus_tag_to_feature_list.iteritems():
        assert len(feature_list) > 1, (
                "Not enough features for locus_tag %s.\nFeatures: %s\n" %
                        (locus_tag, feature_list))

    # Write the result.
    SeqIO.write(recoli_phage_record, FINAL_RECOLI_LAMBDA_PROPHAGE, 'genbank')


if __name__ == '__main__':
    main()

