"""
Script that converts a Genbank file to the five column format required by NCBI
for submitting a genome to Genbank. Yes, they do not accept .genbank format.

See: http://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html
"""

import re

from Bio import SeqIO

from biopython_util import does_feature_end_with_valid_stop_codon
from biopython_util import does_feature_start_with_valid_start_codon
from biopython_util import maybe_get_feature_gene
from biopython_util import update_feature_note
from refactor_checker import check_feature_for_internal_stops


TYPES_WITH_LOCUS_TAG = set(['rRNA', 'tRNA', 'CDS', 'mat_peptide', 'tmRNA',
        'STS', 'gene', 'ncRNA'])

GO_TERM_KEY_SET = set(['GO_component', 'GO_process', 'GO_function'])

DB_XREF_KEY = 'db_xref'

OLD_LOCUS_TAG_KEY = 'old_locus_tag'


class LocusTagGenerator(object):
    """Generates incrementing locus tags.

    Implements the generator pattern.
    """

    def __init__(self, locus_prefix):
        # Count from 0.
        self.current_tag_num = 0

        # Prefix for each locus tag.
        self.locus_prefix = locus_prefix

        # Structure for storing the loci for particular genes since we usually
        # have a 'gene' and 'CDS' feature for a given gene.
        self.locus_dict = {}

    def next(self, feature):
        locus_id = self._get_locus_id(feature)
        if locus_id and locus_id in self.locus_dict:
            return self.locus_dict[locus_id]
        current_tag  = self.locus_prefix + '_' + str(self.current_tag_num).zfill(4)
        self.current_tag_num += 1
        if locus_id:
            self.locus_dict[locus_id] = current_tag
        return current_tag

    def _get_locus_id(self, feature):
        """Return the locus id for the feature.

        This helps determine whether to create a new locus_tag. Generally,
        every 'CDS' has a corresponding 'gene' feature.
        """
        # If the feature already has a locus (e.g. existing geneome getting
        # recoded), then use that as the locus id.
        if 'locus_tag' in feature.qualifiers:
            return feature.qualifiers['locus_tag'][0]

        feature_gene = maybe_get_feature_gene(feature)
        if feature_gene is None:
            return None
        return (feature_gene + '_' + str(feature.location.start) + '_' +
                str(feature.location.end))


def convert_genbank_to_five_column(input_genbank_path, output_five_col_path,
        override_record_name, locus_prefix, ignore_feature_types=['source']):
    """Performs the conversion.

    The input Genbank must contain a single genome record.
    """
    genome_record = SeqIO.read(input_genbank_path, 'genbank')

    # Reduce the features to only the minimum viable set.
    minimum_viableize_features(genome_record)

    # Generate locus_tags.
    locus_tag_generator = LocusTagGenerator(locus_prefix)
    update_locus_tags(genome_record, locus_tag_generator, ignore_feature_types)

    with open(output_five_col_path, 'w') as output_fh:
        # First write the record name.
        output_fh.write('>Feature ' + override_record_name + '\n')
        for feature in genome_record.features:
            if feature.type in ignore_feature_types:
                continue
            _write_single_feature(feature, output_fh)


def minimum_viableize_features(genome_record):
    """Since Genbank has a lot of strict rules to adhere to, we want to be able
    to initially create a minimum-viable Genbank which, while perhaps not as
    complete as it could be, does pass all of the rules.

    We can then go back and open this up to more data until we make Genbank happy.
    """
    # Leave these genes unchanged. I've manually verified them.
    MG1655_OKAY = set([
        # Contains frame-shift necessary for RF-2.
        'prfB',

        # Internal stop codon intentionally introduced for evaluating
        # recombination frequency
        'lacZ',
    ])

    # Get rid of all qualifiers except the following.
    QUALIFIERS_TO_KEEP = set([
        'gene',
        'locus_tag',
        'codon_start',
        'transl_table',
        'transl_except',
        'product',
        'note',
        'pseudo',
        'ncRNA_class',
        'mobile_element_type',
        'experiment'
    ])

    # Variables for logging.
    num_pseudo_features_added = 0

    updated_features = []
    for feature in genome_record.features:
        pseudo_added = False

        # If manually verified, just keep it.
        gene = maybe_get_feature_gene(feature)
        if gene and gene in MG1655_OKAY:
            # Manual handling for lacZ
            if gene == 'lacZ':
                if feature.type in ['CDS', 'mat_peptide']:
                    # Drop CDS feature.
                    continue
                elif feature.type == 'gene':
                    update_feature_note(feature,
                            "Internal stop codon introduced intentionally for benchmarking recombination frequency by Lajoie et al., 2013; version shown in this record is in 'off' state and contains the stop codon.")
        else:
            # Remove all but qualifiers that we are keeping.
            minimum_qualifiers = {}
            for key in QUALIFIERS_TO_KEEP:
                if key in feature.qualifiers:
                    minimum_qualifiers[key] = feature.qualifiers[key]
            feature.qualifiers = minimum_qualifiers

            # Make any features with internal stops pseudogenes.
            if feature.type == 'CDS':
                if (not 'pseudo' in feature.qualifiers and
                        len(check_feature_for_internal_stops(
                                feature, genome_record)) > 0):
                    feature.qualifiers['pseudo'] = ['']
                    update_feature_note(feature,
                            "Known internal stop codons resulting from off-target mutations introduced during UAG recoding in Lajoie et al., 2013.")
                    pseudo_added = True

                # Make any features with no in-frame stop codon pseudogenes.
                if not does_feature_end_with_valid_stop_codon(feature,
                        genome_record):
                    feature.qualifiers['pseudo'] = ['']
                    update_feature_note(feature,
                            "Known absence of valid in-frame stop codon resulting from off-target mutations introduced during UAG recoding in Lajoie et al., 2013.")
                    pseudo_added = True

                # Check invalid stop codons.
                if not does_feature_start_with_valid_start_codon(feature,
                        genome_record):
                    feature.qualifiers['pseudo'] = ['']
                    update_feature_note(feature,
                            "Known invalid start codon resulting from off-target mutations introduced during UAG recoding in Lajoie et al., 2013.")
                    pseudo_added = True

                # Make sure pseudo genes do not have a product.
                if 'pseudo' in feature.qualifiers:
                    if 'product' in feature.qualifiers:
                        del feature.qualifiers['product']
                else:
                    # Not pseudo, so if missing product, add a placeholder.
                    assert 'product' in feature.qualifiers, str(feature)

        if pseudo_added:
            num_pseudo_features_added += 1 # logging

        updated_features.append(feature)
    genome_record.features = updated_features
    print 'Added %d pseudo features.' % num_pseudo_features_added


def update_locus_tags(genome_record, locus_tag_generator,
        ignore_feature_types=['source']):
    """Re-indexes the locus_tags on the features in the genome_record.
    """
    updated_features = []
    for feature in genome_record.features:
        if feature.type in ignore_feature_types:
            continue
        elif not 'locus_tag' in feature.qualifiers:
            # Skip adding locus_tags where they didn't exist before, for now.
            pass
        elif feature.type in TYPES_WITH_LOCUS_TAG:
            feature.qualifiers['locus_tag'] = [locus_tag_generator.next(feature)]
        else:
            # Remove locus_tag from qualifiers as there should no longer
            # be one for whatever reason.
            del feature.qualifiers['locus_tag']
        updated_features.append(feature)
    genome_record.features = updated_features


def _write_single_feature(feature, output_fh, ignore_feature_keys=['translation']):
    """Writes the rows for a single feature.
    """
    feature_lines = []

    # HACK: Special handling to deal with prfB slip until I can make
    # this more general.
    is_prfB_hack = False
    if feature.type == 'CDS':
        gene = maybe_get_feature_gene(feature)
        if gene and gene == 'prfB':
            is_prfB_hack = True
            feature_lines.append('3038585\t3038511\tCDS')
            feature_lines.append('3038509\t3037487')

    if not is_prfB_hack:
        start = feature.location.start
        end = feature.location.end

        # Adjust feature indeces. The .tbl format is 1-indexed, end-inclusive
        # while BioPython 0-indexed, is upper-bound-exclusive.
        # So just updated the start. End is updated implicitly.
        start += 1

        # Switch start and end.
        # The way to represent a feature on the negative strand in this format is to
        # reverse the order in which the positions are recorded.
        if feature.strand != 1:
            temp = end
            end = start
            start = temp

        # Cast these to string, since below we might add a '^' character to the
        # start position in order to represent a deletion. We cast end for
        # consistency.
        start = str(start)
        end = str(end)

        # HACK: Deletions are represented by adding a '^' after the first base.
        if 'source_deletion_interval' in feature.qualifiers:
            start += '^'
            del feature.qualifiers['source_deletion_interval']

        feature_position_line = '%s\t%s\t%s' % (start, end, feature.type)
        feature_lines.append(feature_position_line)

    for key in feature.qualifiers:
        if key in ignore_feature_keys:
            continue
        for value in feature.qualifiers[key]:
            if key in GO_TERM_KEY_SET:
                # TODO: Support this.
                # feature_lines.extend(_handle_go_term(key, value))
                continue
            elif key == DB_XREF_KEY:
                # TODO: Support as much of this as possible.
                continue
            elif key == OLD_LOCUS_TAG_KEY:
                # TODO: Figure out what TODO.
                continue
            else:
                feature_lines.append(
                        '\t\t\t%s\t%s' % (key, value))
    feature_string = '\n'.join(feature_lines) + '\n'
    output_fh.write(feature_string)


# Regular expressions for parsing GO terms.
GO_TERM_ID_REGEX = re.compile('(?<=GO:)\w+')

GO_TERM_DESC_REGEX = re.compile('(?<=- ).*')

def _handle_go_term(key, value):
    """Returns list of lines since multiple GO terms are possible here.
    """
    go_id_part = GO_TERM_ID_REGEX.search(value).group()
    desc_part = GO_TERM_DESC_REGEX.search(value).group()
    value_part = desc_part + '|' + go_id_part
    return '\t\t\t%s\t%s' % (key, value_part)


if __name__ == '__main__':
    convert_genbank_to_five_column(
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/mg1655/mg1655.genbank',
            '/home/glebk/Projects/churchlab/genome-refactor/data/genomes/rec1_c321d/test_sqn_prep_mg1655/mg1655.tbl',
            'unassigned_id',
            'N840')
