"""
Tests for BioPython whose purpose is to experiment with various techniques
for using it.

For example, for refactoring, we want to make changes to specific features.
This leads to questions like whether changing the position of a particular
feature changes the position of all features.
"""

import os
import unittest

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from biopython_util import add_custom_seq_record_methods
from biopython_util import add_feature_to_seq_record
from biopython_util import calc_interval_list_to_features_overlapped
from biopython_util import delete_single_base
from biopython_util import find_features_ending_at
from biopython_util import find_features_starting_at
from biopython_util import get_codon_feature_location
from biopython_util import get_feature_by_id
from biopython_util import get_genome_record
from biopython_util import get_region_codon_indeces_in_feature
from biopython_util import insert_sequence
from biopython_util import InsertType
from biopython_util import remove_gene_features
from biopython_util import shift_seq_feature_right_with_head_copy
from biopython_util import swap_region_seq
from refactor_context import RefactorContext


PWD = os.path.dirname(os.path.realpath(__file__ ))
TEST_DATA_DIR = os.path.abspath(PWD + '/../test_data')
TEST_GBK_FILE = os.path.join(TEST_DATA_DIR, 'aggr_in.g')

class TestBioPythonUtil(unittest.TestCase):

    def _assert_feature_seq(self, feature, seq_record, feature_id_to_seq_map):
        """Utility method for assert that a sequence is preserved in a
        SeqRecord.
        """
        original_seq = feature_id_to_seq_map[feature.id]
        self.assertEqual(original_seq, str(feature.extract(seq_record).seq))

    def test_get_genome_record(self):
        genome_record = get_genome_record(TEST_GBK_FILE)

    def test_insert_sequence_no_features(self):
        seq = Seq('ATGTTTGGG', generic_dna)
        seq_record = SeqRecord(seq)

        insert_seq = Seq('CCC', generic_dna)
        insert_start_position = 6
        updated_seq_record = insert_sequence(
                seq_record, insert_seq, insert_start_position)
        self.assertEqual('ATGTTTCCCGGG', str(updated_seq_record.seq))


    def test_insert_sequence_update_features(self):
        seq = Seq('ATGTTTGGGAAATTT', generic_dna)
        seq_record = SeqRecord(seq)

        feature_id_to_seq_map = {
            1: 'ATG',
            2: 'GGG',
            3: 'AAATTT'
        }

        def _assert_feature_seq(feature, seq_record):
            original_seq = feature_id_to_seq_map[feature.id]
            self.assertEqual(original_seq, str(feature.extract(seq_record).seq))

        # Make a mix of features, some before, some after the insertion.
        feature_1_loc = FeatureLocation(0, 3)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)
        _assert_feature_seq(feature_1, seq_record)

        feature_2_loc = FeatureLocation(6, 9)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)
        _assert_feature_seq(feature_2, seq_record)

        feature_3_loc = FeatureLocation(9, 15)
        feature_3 = SeqFeature(feature_3_loc, type='CDS', id=3)
        seq_record.features.append(feature_3)
        _assert_feature_seq(feature_3, seq_record)

        # The sequence to insert.
        insert_seq = Seq('CCC', generic_dna)
        insert_start_position = 6

        # Perform the insertion.
        updated_seq_record = insert_sequence(
                seq_record, insert_seq, insert_start_position)

        # Assert conditions that we want to preserve.
        self.assertEqual('ATGTTTCCCGGGAAATTT', str(updated_seq_record.seq))
        self.assertEqual(
                len(seq_record.features), len(updated_seq_record.features))
        for feature in updated_seq_record.features:
            _assert_feature_seq(feature, updated_seq_record)


    def test_insert_sequence_unexpected_inside_feature(self):
        seq = Seq('ATGTTTGGGAAATTT', generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(6, 9)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=2)
        seq_record.features.append(feature_1)

        # The sequence to insert.
        insert_seq = Seq('CCC', generic_dna)
        insert_start_position = 7

        # Try the insertion.
        self.assertRaises(
                ValueError, insert_sequence,
                seq_record, insert_seq, insert_start_position)


    def test_delete_single_base(self):
        feature_1_seq = 'ATGAAAGGGATG'
        blah_1 = 'AGAGATAGA'
        feature_2_seq = 'ATGTTT'
        blah_2 = 'GACCCGCTA'
        feature_3_seq = 'CCC'
        blah_3 = 'TTTTT'

        expected_feature_id_to_seq_map = {
            1: feature_1_seq,
            2: feature_2_seq,
            3: feature_3_seq,
        }

        whole_seq = Seq(
                feature_1_seq + blah_1 + feature_2_seq + blah_2 + feature_3_seq +
                        blah_3,
                generic_dna)
        seq_record = SeqRecord(whole_seq)

        feature_1_loc = FeatureLocation(0, len(feature_1_seq))
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        feature_2_start = len(feature_1_seq) + len(blah_1)
        feature_2_end = feature_2_start + len(feature_2_seq)
        feature_2_loc = FeatureLocation(feature_2_start, feature_2_end)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)

        feature_3_start = feature_2_end + len(blah_2)
        feature_3_end = feature_3_start + len(feature_3_seq)
        feature_3_loc = FeatureLocation(feature_3_start, feature_3_end)
        feature_3 = SeqFeature(feature_3_loc, type='CDS', id=3)
        seq_record.features.append(feature_3)

        # Test the features before making the change.
        for feature in seq_record.features:
            self._assert_feature_seq(
                    feature, seq_record, expected_feature_id_to_seq_map)

        # Make the change.
        delete_single_base(seq_record, feature_2_start + 1, 'T')
        expected_feature_id_to_seq_map[2] = 'AGTTT'

        # Test the features after making the change.
        for feature in seq_record.features:
            self._assert_feature_seq(
                    feature, seq_record, expected_feature_id_to_seq_map)


    def test_shift_seq_feature_right_with_head_copy(self):
        seq = Seq('ATGAAAGGGATGTTTCCC', generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_seq = 'ATGAAAGGGATG'
        feature_2_seq = 'ATGTTT'
        feature_3_seq = 'CCC'
        expected_feature_id_to_seq_map = {
            1: feature_1_seq,
            2: feature_2_seq,
            3: feature_3_seq,
            '2_overlap_head_cp': 'ATG'
        }

        def _assert_feature_seq(feature, seq_record):
            original_seq = expected_feature_id_to_seq_map[feature.id]
            self.assertEqual(original_seq, str(feature.extract(seq_record).seq))

        feature_1_loc = FeatureLocation(0, 12)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)
        _assert_feature_seq(feature_1, seq_record)

        feature_2_loc = FeatureLocation(9, 15)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id=2)
        seq_record.features.append(feature_2)
        _assert_feature_seq(feature_2, seq_record)

        feature_3_loc = FeatureLocation(15, 18)
        feature_3 = SeqFeature(feature_3_loc, type='CDS', id=3)
        seq_record.features.append(feature_3)
        _assert_feature_seq(feature_3, seq_record)

        original_num_features = len(seq_record.features)

        updated_seq_record = shift_seq_feature_right_with_head_copy(
                feature_2, seq_record, 12)

        self.assertEqual('ATGAAAGGGATGATGTTTCCC', str(updated_seq_record.seq))

        # A new feature was created to represent the pulled apart feature.
        self.assertEqual(original_num_features + 1,
                len(updated_seq_record.features))
        for feature in updated_seq_record.features:
            _assert_feature_seq(feature, updated_seq_record)


    def test_find_features_ending_at(self):
        before = 'TACTAGTCGAT'
        feature_1_seq = 'ATGAAAGGGATG'
        after = 'CTATCTAGCTAGCT'
        whole_seq = Seq(before + feature_1_seq + after, generic_dna)
        seq_record = SeqRecord(whole_seq)

        feature_1_loc = FeatureLocation(
                len(before),
                len(before) + len(feature_1_seq),
                strand=1
        )
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)

        # A variety of tests.
        self.assertEquals(set([]), set(find_features_ending_at(0, seq_record)))
        self.assertEquals(
                set([feature_1]), set(find_features_ending_at(
                        feature_1_loc.end, seq_record)))
        self.assertEquals(
                set(), set(find_features_ending_at(
                        feature_1_loc.end, seq_record, ['misc'])))


    def test_find_features_starting_at(self):
        before = 'TACTAGTCGAT'
        feature_1_seq = 'ATGAAAGGGATG'
        after = 'CTATCTAGCTAGCT'
        whole_seq = Seq(before + feature_1_seq + after, generic_dna)
        seq_record = SeqRecord(whole_seq)

        feature_1_loc = FeatureLocation(
                len(before),
                len(before) + len(feature_1_seq),
                strand=1
        )
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)

        self.assertEquals(set([]), set(find_features_starting_at(0, seq_record)))
        self.assertEquals(
                set([feature_1]), set(find_features_starting_at(
                        feature_1_loc.start, seq_record)))
        self.assertEquals(
                set(), set(find_features_starting_at(
                        feature_1_loc.end, seq_record, ['misc'])))


    def test_get_region_codon_indeces_in_feature__positive_strand(self):
        """Test that the region captures the appropriate indeces in the feature.
        """
        feature_1_loc = FeatureLocation(0, 12, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')

        EXPECTED_RESULT = [0, 1]
        for region in [(0, 5), (0, 6), (1, 6)]:
            result = get_region_codon_indeces_in_feature(feature_1, region)
            self.assertEqual(set(EXPECTED_RESULT), set(result))

        EXPECTED_RESULT = [0, 1, 2]
        for region in [(0, 7), (1, 7)]:
            result = get_region_codon_indeces_in_feature(feature_1, region)
            self.assertEqual(set(EXPECTED_RESULT), set(result))

        feature_2_loc = FeatureLocation(1, 12, strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id='1')
        EXPECTED_RESULT = [0, 1]
        for region in [(0, 5), (0, 6), (1, 6), (0, 7), (1, 7)]:
            result = get_region_codon_indeces_in_feature(feature_2, region)
            self.assertEqual(set(EXPECTED_RESULT), set(result))

        EXPECTED_RESULT = [0, 1, 2]
        for region in [(0, 8), (1, 8)]:
            result = get_region_codon_indeces_in_feature(feature_2, region)
            self.assertEqual(set(EXPECTED_RESULT), set(result))


    def test_get_region_codon_indeces_in_feature__negative_strand(self):
        """Test that the region captures the appropriate indeces in the feature.
        """
        feature_1_loc = FeatureLocation(0, 12, strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')

        EXPECTED_RESULT = [2, 3]
        for region in [(0, 5), (0, 6), (1, 6)]:
            result = get_region_codon_indeces_in_feature(feature_1, region)
            self.assertEqual(set(EXPECTED_RESULT), set(result))

        EXPECTED_RESULT = [1, 2, 3]
        for region in [(0, 7), (1, 7)]:
            result = get_region_codon_indeces_in_feature(feature_1, region)
            self.assertEqual(set(EXPECTED_RESULT), set(result))


    def test_calc_interval_list_to_features_overlapped(self):
        INTERVAL_LIST = [
                (5, 18),
                (25, 70),
                (100, 105),
                (125, 145),
                (500, 505)
        ]

        FEATURE_LIST = [
                SeqFeature(FeatureLocation(0, 12)),
                SeqFeature(FeatureLocation(20, 40)),
                SeqFeature(FeatureLocation(45, 67)),
                SeqFeature(FeatureLocation(100, 120)),
                SeqFeature(FeatureLocation(200, 500)),
        ]

        result = calc_interval_list_to_features_overlapped(
                INTERVAL_LIST, FEATURE_LIST)

        self.assertEqual(len(INTERVAL_LIST), len(result))

        self.assertEqual(1, len(result[0]))
        self.assertTrue(FEATURE_LIST[0] in result[0])

        self.assertEqual(2, len(result[1]))
        self.assertTrue(FEATURE_LIST[1] in result[1])
        self.assertTrue(FEATURE_LIST[2] in result[1])

        self.assertEqual(1, len(result[2]))
        self.assertTrue(FEATURE_LIST[3] in result[2])

        self.assertEqual(0, len(result[3]))

        self.assertEqual(0, len(result[4]))


    def test_swap_region_seq(self):
        SEQ = Seq('ATGTTTGGG', generic_dna)
        SEQ_RECORD = SeqRecord(SEQ)
        REFACTOR_CONTEXT = RefactorContext(SEQ_RECORD)

        NEW_SEQ = 'TTAGGA'

        updated_seq_record = swap_region_seq(
                REFACTOR_CONTEXT, 3, 'TTTGGG', NEW_SEQ)

        self.assertEqual('ATGTTAGGA', str(updated_seq_record.seq))


    def test_swap_region_seq__fuzzy(self):
        SEQ = Seq('ATGTTTGGG', generic_dna)
        SEQ_RECORD = SeqRecord(SEQ)
        REFACTOR_CONTEXT = RefactorContext(SEQ_RECORD)

        NEW_SEQ = 'TTAGGA'

        updated_seq_record = swap_region_seq(
                REFACTOR_CONTEXT, 1, 'TTTGGG', NEW_SEQ, fuzzy=True)

        self.assertEqual('ATGTTAGGA', str(updated_seq_record.seq))


    def test_add_feature_to_seq_record__regular(self):
        SEQ = Seq('ATGTTTGGGTAGAGTA', generic_dna)
        seq_record = SeqRecord(SEQ)
        FEATURE_1_ID = '1'
        FEATURE_1_LOC = FeatureLocation(4, 7)
        feature_1 = SeqFeature(FEATURE_1_LOC, type='CDS', id=FEATURE_1_ID)
        add_feature_to_seq_record(seq_record, feature_1)

        lookup_feature_1 = get_feature_by_id(seq_record, FEATURE_1_ID)
        self.assertEqual(feature_1, lookup_feature_1)


    def test_add_feature_to_seq_record__feature_index_map(self):
        SEQ = Seq('ATGTTTGGGTAGAGTA', generic_dna)
        seq_record = SeqRecord(SEQ)
        add_custom_seq_record_methods(seq_record)
        FEATURE_1_ID = '1'
        FEATURE_1_LOC = FeatureLocation(4, 7)
        feature_1 = SeqFeature(FEATURE_1_LOC, type='CDS', id=FEATURE_1_ID)
        add_feature_to_seq_record(seq_record, feature_1)

        lookup_feature_1 = get_feature_by_id(seq_record, FEATURE_1_ID)
        self.assertEqual(feature_1, lookup_feature_1)


    def test_get_codon_feature_location(self):
        SEQ = Seq('CCCCCCATGTTTGGGTAGAGTA', generic_dna)
        seq_record = SeqRecord(SEQ)
        FEATURE_1_ID = '1'
        FEATURE_1_LOC = FeatureLocation(6, 18)
        feature_1 = SeqFeature(FEATURE_1_LOC, type='CDS',
                id=FEATURE_1_ID, strand=1)
        add_feature_to_seq_record(seq_record, feature_1)

        CODON_0_FEATURE_LOCATION = get_codon_feature_location(feature_1, 0)
        self.assertEqual(6, CODON_0_FEATURE_LOCATION.start)
        self.assertEqual(9, CODON_0_FEATURE_LOCATION.end)

        CODON_1_FEATURE_LOCATION = get_codon_feature_location(feature_1, 1)
        self.assertEqual(9, CODON_1_FEATURE_LOCATION.start)
        self.assertEqual(12, CODON_1_FEATURE_LOCATION.end)


    def test_get_codon_feature_location__reverse_strand(self):
        FEATURE_SEQ_RAW = reverse_complement('ATGTTTGGGTAG')
        SEQ = Seq('CCCCCC' + FEATURE_SEQ_RAW + 'AGTA', generic_dna)
        seq_record = SeqRecord(SEQ)
        FEATURE_1_ID = '1'
        FEATURE_1_LOC = FeatureLocation(6, 18)
        feature_1 = SeqFeature(FEATURE_1_LOC, type='CDS',
                id=FEATURE_1_ID, strand=-1)
        add_feature_to_seq_record(seq_record, feature_1)

        # Sanity check for the feature sequence.
        feature_seq = str(feature_1.extract(seq_record.seq))
        self.assertEqual('ATG', feature_seq[0:3])
        self.assertEqual('TTT', feature_seq[3:6])

        CODON_0_FEATURE_LOCATION = get_codon_feature_location(feature_1, 0)
        self.assertEqual(15, CODON_0_FEATURE_LOCATION.start)
        self.assertEqual(18, CODON_0_FEATURE_LOCATION.end)

        CODON_1_FEATURE_LOCATION = get_codon_feature_location(feature_1, 1)
        self.assertEqual(12, CODON_1_FEATURE_LOCATION.start)
        self.assertEqual(15, CODON_1_FEATURE_LOCATION.end)

    def test_remove_gene_features(self):
        from refactor_config import GENOMES_DIR
        MG1655_GENBANK = os.path.join(GENOMES_DIR, 'mg1655', 'mg1655.genbank')
        GENE_DICT = {'prfA': {'types': ('CDS', 'gene')}}
        with open(MG1655_GENBANK) as fh:
            genome_record = SeqIO.read(fh, 'genbank')
        remove_gene_features(genome_record, GENE_DICT)

        for feature in genome_record.features:
            if ('gene' in feature.qualifiers and
                    feature.qualifiers['gene'][0] == 'prfA'):
                self.fail("Feature %s not removed." % feature)


if __name__ == '__main__':
    unittest.main()
