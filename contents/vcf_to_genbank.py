"""
Methods for impressing the changes in a VCF file onto an existing
Genbank file.

NOTES:
    * Look into insertion annotations being 1 off (specifically galk insertion at SIR.30.31i
    * Annotation of deletion d_mg1655_66564_66733  partial deletion of repeat region 66564 66733 seems to be one off.
"""

import copy
import csv
import os
import pickle

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import vcf

from biopython_util import add_feature_to_seq_record
from biopython_util import delete_interval
from biopython_util import insert_sequence_and_update_features
from refactor_config import GENOMES_DIR

###############################################################################
# File locations
###############################################################################

REC1_C321D_ROOT = os.path.join(GENOMES_DIR, 'rec1_c321d')


###############################################################################
# Constants
###############################################################################

VARIANT_ANNOTATION_TYPE = 'variation'

MAX_REPLACE_CHARS = 12


###############################################################################
# Helper objects used by the main procedure
###############################################################################

class RuntimeLiftover(object):
    """Object that aids in dynamically updating a genome record with a list of
    positions that are relative to the original genome record.

    An example of a case where this is useful is when you are creating a
    mapping from a vcf record which shows SNPs and other variants relative
    to a reference genome.

    For example, say you have two insertions:
        * A - position: 100, size: 3
        * B - position: 200, size: 4

    When we introduce the first insertion, the frame of the underlying genome
    has shifted, so that the second insertion should really be added at
    position 200 + 3.
    """

    def __init__(self, original_genome_record):
        """Constructor.
        """
        # The original record. This remains unchanged throughout the
        # mapping process.
        # TODO: Do we even need to be keeping track of this? Or are intervals
        # sufficient?
        self.source_genome_record = original_genome_record

        # Structure that maintains the mapping of intervals.
        # Let's say we guarantee that it's sorted and exclusive, for now.
        # NOTE: Each interval maintains a one-to-one mapping and so is
        # inclusive bounds on both ends.
        self._interval_mapping = self._initialize_interval_mapping()


    def _initialize_interval_mapping(self):
        """Initializes the interval mapping.
        """
        # The initial mapping is a list with a single element, which is
        # a pair of tuples representing the correspondence between the original
        # whole sequence interval and a copy of itself.
        original_interval = (0, len(self.source_genome_record) - 1)
        initial_mapping_pair = (original_interval, copy.copy(original_interval))
        return [initial_mapping_pair]


    @classmethod
    def from_pickled_intervals(cls, original_genome_record, pickle_dest):
        """Factory method that creates a RuntimeLiftover object and sets
        the intervals from a pickle file.

        The original genome record still has to be provided.
        """
        runtime_liftover = cls(original_genome_record)
        with open(pickle_dest) as pickle_fh:
            runtime_liftover._interval_mapping = pickle.load(pickle_fh)
        return runtime_liftover


    def pickle_interval_mapping(self, pickle_dest):
        """Pickle the interval mapping and write to file.

        This is useful for debugging intervals and developing other output
        formats.
        """
        with open(pickle_dest, 'w') as pickle_fh:
            pickle.dump(self._interval_mapping, pickle_fh)


    def write_chain_file(self, chain_file_dest):
        """Writes the current state of _interval_mapping in the UCSC
        liftover chain file format.

        See: http://genome.ucsc.edu/goldenPath/help/chain.html
        """
        with open(chain_file_dest, 'w') as chain_file_fh:
            # Write the heading.
            chain_file_fh.write('chain\n')

            # Each row is of the form 'size dt dq', separated by spaces.
            #   * size: the size of the ungapped alignment
            #   * dt: the difference between the end of this block and the
            #         beginning of the next block (reference sequence)
            #   * dq: the difference between the end of this block and the
            #         beginning of the next block (query sequence)
            # NOTE: The last line of the alignment section contains only one
            #       number: the ungapped alignment size of the last block.
            interval_index = 0
            num_interval_mappings = len(self._interval_mapping)
            for interval_index in range(num_interval_mappings):
                # I am using the names '*reference*' and '*query*' in the sense
                # that the chain file uses them, where query sequence is the one
                # whose coordinates we are generally trying to convert into the
                # frame of the target. Typically the query sequence is the one
                # we mapped the VCF changes on top of.
                (current_reference_interval, current_query_interval) = (
                        self._interval_mapping[interval_index])
                size = bases_in_interval(current_reference_interval)
                next_interval_index = interval_index + 1
                if next_interval_index < num_interval_mappings:
                    (next_reference_interval, next_query_interval) = (
                            self._interval_mapping[next_interval_index])
                    dt = (next_reference_interval[0] -
                            current_reference_interval[1] - 1)
                    dq = (next_query_interval[0] -
                            current_query_interval[1] - 1)
                    chain_file_fh.write('%d %d %d\n' % (size, dt, dq))
                else:
                    # This is the last line. Just write the block size.
                    chain_file_fh.write('%d\n' % (size,))


    def convert_source_position_to_target(self, source_position, or_next=False):
        """Converts a single position in the source genome to the corresponding
        position in the target genome (the one being updated).

        Args:
            source_position: Position in the source genome. (0-indexed).
            or_next: If True, when no direct mapping, return the next position.

        Returns:
            The position in the target genome, or None if mapping failed.
        """
        assert isinstance(source_position, int), "source_position must be int."

        # For now, the algorithm is to first search every interval in the
        # internal interval mapping data structure until we find the one that
        # the source position lies in, and then find the corresponding target
        # position by using the relative offset of the source position within
        # the interval.
        for mapping_index in xrange(len(self._interval_mapping)):
            source_interval, target_interval = self._interval_mapping[
                    mapping_index]
            if source_interval[0] <= source_position <= source_interval[1]:
                interval_index = source_position - source_interval[0]
                return target_interval[0] + interval_index

            if or_next and source_position < source_interval[0]:
                return self.convert_source_position_to_target(
                        source_interval[0], or_next=or_next)
        return None


    def convert_target_position_to_source(self, target_position):
        """Converts a single position in the target genome to the corresponding
        position in the source genome (the one being updated).

        Similar, but more limited than convert_source_position_to_target().

        Args:
            target_position: Position in the target genome. (0-indexed).

        Returns:
            The position in the source genome, or None if mapping failed.
        """
        assert isinstance(target_position, int), "target_position must be int."

        for mapping_index in xrange(len(self._interval_mapping)):
            source_interval, target_interval = self._interval_mapping[
                    mapping_index]
            if target_interval[0] <= target_position <= target_interval[1]:
                interval_index = target_position - target_interval[0]
                return source_interval[0] + interval_index
        return None


    def handle_insertion(self, variant_data):
        """Handles an insertion with the given data spec.
        """
        # Create a new interval mapping and replace the member attribute
        # at the end.
        new_interval_mapping = []

        # Parse the insert data object.
        insert_position = variant_data['position']
        insert_sequence = variant_data['sequence']
        len_insert_sequence = len(insert_sequence)

        # We use a state machine strategy to first find the interval
        # to insert, and then update all downstream target intervals.
        STATE_SEARCHING = 'SEARCHING'
        STATE_UPDATING_TARGET_DOWNSTREAM = 'UPDATING_TARGET_DOWNSTREAM'

        state = STATE_SEARCHING
        for idx, (source_interval, target_interval) in enumerate(
                self._interval_mapping):
            if state == STATE_SEARCHING:
                if source_interval[0] <= insert_position <= source_interval[1]:
                    insert_position_index = insert_position - source_interval[0]

                    # The source simply gets split.
                    new_source_interval_upstream = (source_interval[0],
                            insert_position - 1)
                    new_source_interval_downstream = (insert_position,
                            source_interval[1])

                    # The target gets split, with the downstream interval
                    # shifted by the size of the insertion sequence.
                    new_target_interval_upstream = (target_interval[0],
                            target_interval[0] + insert_position_index - 1)
                    new_target_interval_downstream = (target_interval[0] +
                            insert_position_index + len_insert_sequence,
                            target_interval[1] + len_insert_sequence)

                    # Append the split sequence pairs.
                    new_interval_mapping.append((new_source_interval_upstream,
                            new_target_interval_upstream))
                    new_interval_mapping.append((new_source_interval_downstream,
                            new_target_interval_downstream))

                    # Update the state for remaining iterations.
                    state = STATE_UPDATING_TARGET_DOWNSTREAM
                elif insert_position < source_interval[0]:
                    # The insert_position was deleted. Shift the target
                    # interval downstream by the size of the insertion.
                    new_source_interval = (
                            source_interval[0],
                            source_interval[1])
                    new_target_interval = (
                            target_interval[0] + len(insert_sequence),
                            target_interval[1] + len(insert_sequence))
                    assert (bases_in_interval(new_source_interval) ==
                            bases_in_interval(new_target_interval))
                    new_interval_mapping.append((new_source_interval,
                            new_target_interval))
                    state = STATE_UPDATING_TARGET_DOWNSTREAM
                else:
                    new_interval_mapping.append(
                            (source_interval, target_interval))
            else:
                # Shift all remaining target intervals.
                new_target_interval = (
                        target_interval[0] + len_insert_sequence,
                        target_interval[1] + len_insert_sequence)
                new_interval_mapping.append((source_interval,
                        new_target_interval))

        if state == STATE_SEARCHING:
            raise RuntimeError("Error updating RuntimeLiftover with %s", (
                    str(variant_data,)))

        self._interval_mapping = new_interval_mapping


    def handle_deletion(self, variant_data):
        """Handles a deletion with the given data spec.

        Args:
            variant_data: Dictionary with keys:
                * interval: A two-tuple representing pythonic interval for the
                    deletion, i.e. (inclusive_start, exclusive_end).
                    e.g. (100, 102) is a deletion of the 2 bases at positions
                    100 and 101. Relative to the source genome.
        """
        interval = variant_data['interval']

        for source_position in range(*interval):

            delete_position = self.convert_source_position_to_target(
                    source_position, or_next=True)

            delete_position = interval[0]
            delete_interval_size = 1

            # Create a new interval mapping and replace the member attribute
            # at the end.
            new_interval_mapping = []

            # We use a state machine strategy to first find the interval
            # to delete, and then update all downstream target intervals.
            STATE_SEARCHING = 'SEARCHING'
            STATE_UPDATING_TARGET_DOWNSTREAM = 'UPDATING_TARGET_DOWNSTREAM'

            state = STATE_SEARCHING
            for source_interval, target_interval in self._interval_mapping:
                if state == STATE_SEARCHING:
                    if source_interval[0] <= delete_position <= source_interval[1]:
                        delete_position_index = delete_position - source_interval[0]

                        # The source simply gets split, dropping a base.
                        new_source_interval_upstream = (source_interval[0],
                                delete_position - 1)
                        new_source_interval_downstream = (
                                delete_position + delete_interval_size,
                                source_interval[1])

                        # The target gets split, including the position, but
                        # reducing the size of this and all following intervals.
                        new_target_interval_upstream = (target_interval[0],
                                target_interval[0] + delete_position_index - 1)
                        new_target_interval_downstream = (
                                target_interval[0] + delete_position_index,
                                target_interval[1] - delete_interval_size)

                        # Append the split sequence pairs.
                        new_interval_mapping.append((new_source_interval_upstream,
                                new_target_interval_upstream))
                        new_interval_mapping.append((new_source_interval_downstream,
                                new_target_interval_downstream))

                        # Update the state for remaining iterations.
                        state = STATE_UPDATING_TARGET_DOWNSTREAM
                    elif delete_position < source_interval[0]:
                        # The interval this delete_position would have fallen
                        # into has been deleted, effectively delete the first
                        # position in this current interval.
                        new_source_interval = (source_interval[0] + 1,
                                source_interval[1])
                        new_target_interval = (target_interval[0],
                                target_interval[1] - 1)
                        new_interval_mapping.append((new_source_interval,
                                new_target_interval))
                        state = STATE_UPDATING_TARGET_DOWNSTREAM
                    else:
                        new_interval_mapping.append(
                                (source_interval, target_interval))
                else: # state == STATE_UPDATING_TARGET_DOWNSTREAM
                    # Shift all remaining target intervals.
                    new_target_interval = (
                            target_interval[0] - delete_interval_size,
                            target_interval[1] - delete_interval_size)
                    new_interval_mapping.append(
                            (source_interval, new_target_interval))

            if state == STATE_SEARCHING:
                raise RuntimeError("Error updating RuntimeLiftover for %s" % (
                        str(variant_data,)))

            self._interval_mapping = new_interval_mapping


def bases_in_interval(interval):
    """Returns the number of bases in the liftover interval.
    These are inclusive on both ends, not Pythonic that is.
    """
    return interval[1] - interval[0] + 1


class VCFToGenbankMaker(object):
    """Object that encapsulates the logic for updating a genbank file
    with changes from a vcf file.

    Usage:
        1. Construct an instance according to the constructor signature.
        2. Call run().
    """

    def __init__(self, genome_record, vcf_path, sample_id,
            manual_updates_filepath=None):
        """Constructor.
        """
        # Keep a copy of the original genome record.
        self.original_genome_record = copy.deepcopy(genome_record)

        # The record that is mutated as we progress.
        self.genome_record = genome_record

        # Save the path to the vcf. We'll use the vcf.Reader stream reader
        # object when we actually need to handle it.
        self.vcf_path = vcf_path

        # The specific sample in the vcf.
        self.sample_id = sample_id

        # Location with manual updates. Might be None.
        self.manual_updates_filepath = manual_updates_filepath

        # Object used to track the dynamically changing interval mapping
        # positions in the original genome record with respect to which
        # vcf positions were identified to the most up-to-date genome.
        self.runtime_liftover = RuntimeLiftover(self.original_genome_record)


    def run(self, verbose=False, log_file=None):
        """Performs the actual updating.
        """
        # Manually add annotations for TAG.
        add_TAG_annotations(self.genome_record)

        # Add changes made manually.
        if self.manual_updates_filepath:
            if verbose: print '...Handling manual updates...'
            self.handle_manual_updates()
            if verbose: print '...Done handling manual updates.'

        # Extra data to add to features.
        position_to_data_map = get_vcf_metadata()

        # Keep track of which vcf changes were actually made.
        vcf_positions_updated = []

        if verbose: print 'Handling vcf...'
        with open(self.vcf_path) as vcf_fh:
            vcf_reader = vcf.Reader(vcf_fh)
            sample_index = vcf_reader.samples.index(self.sample_id)
            for idx, record in enumerate(vcf_reader):
                if verbose: print idx, record
                # metadata = position_to_data_map.get(record.POS, None)
                metadata = None
                was_change_made = self.handle_vcf_record(
                        record, sample_index, metadata)
                assert isinstance(was_change_made, bool), (
                        "handle_vcf_record() must return a boolean.")
                if was_change_made:
                    vcf_positions_updated.append(record.POS)

        # Write debug output for which changes were actually made.
        if log_file:
            with open(log_file, 'w') as log_fh:
                for pos in vcf_positions_updated:
                    log_fh.write(str(pos) + '\n')


    def handle_manual_updates(self):
        """Adds manual updates from an external file with the following
        tab-separated fields:
            * Cassette
            * Comment
            * LP
            * RP
            * Sequence
        """
        assert self.manual_updates_filepath, "No manual updates specified."

        with open(self.manual_updates_filepath) as fh:
            line = fh.readline() # Skip the header.
            line = fh.readline()
            while line:
                # Parse the arguments.
                args = line.split('\t')
                if not len(args) >= 3:
                    line = fh.readline()
                    continue
                fix_id = args[0]
                note = args[1]
                left_bp = int(args[2].strip())
                right_bp = int(args[3].strip())
                if len(args) == 5 and len(args[4].strip()):
                    seq = args[4].strip()
                else:
                    seq = None

                # Process the args to get the data necessary to make the update.
                pythonic_start = left_bp - 1
                pythonic_end = right_bp - 1
                ref = str(self.original_genome_record.seq[
                        pythonic_start:pythonic_end])
                if seq:
                    alt = seq
                else:
                    alt = ''

                # Make the update for this record.
                self._update_genome_record_for_variant(
                        pythonic_start, ref, alt, note)

                # Continue to next line.
                line = fh.readline()


    # def add_cassette_modifications(self, cassette_modifications_csv):
    #     """Method that allows adding modifications to a cassette relative
    #     to the starting position of the cassette.

    #     This might be desired where the cassettes come from some canonical
    #     source, but in the current context they have modifications.

    #     Args:
    #         cassette_modifications_csv: List of modifications to make, described
    #             with the following columns:
    #                 * cassette_id: Unique id for the cassette as provided in the
    #                     manual updates file.
    #                 * comment
    #                 * position: 1-based index of the mutation start.
    #                 * ref: What was there previously.
    #                 * alt: The alternate value for the position. Currently in
    #                     vcf output format with square brackets.
    #     """
    #     with open(variant_data_csv) as csv_fh:
    #         csv_reader = csv.DictReader(csv_fh)
    #         for row in csv_reader:
    #             pass


    def update_from_variant_data_csv(self, variant_data_csv):
        """Updates the genome given a list of variants in a csv.

        Args:
            variant_data_csv: Path to .csv file containing the following cols:
                Required:
                    * position - Position in the starting genome.
                    * ref - Reference sequence at that position.
                    * alt - The alternative to replace with.
                Optional:
                    * note - A note to add to the data.

        """
        with open(variant_data_csv) as csv_fh:
            csv_reader = csv.DictReader(csv_fh)
            for row in csv_reader:
                pythonic_start = int(row['position']) - 1

                ref = row['ref']

                # NOTE: Hacky way of parsing that works for the way that
                # our particular data looks.
                alt = row['alt'][1:-1]

                if 'note' in row:
                    note = row['note']
                else:
                    note = None

                # Make the update for this record.
                self._update_genome_record_for_variant(
                        pythonic_start, ref, alt, note=note)


    def handle_vcf_record(self, record, sample_index, metadata=None):
        """Decides what to do with a single VCF call.
        """
        # The specific sample for this record.
        # NOTE: The vcf may have been generated across many samples at the
        # same time but this script only operates on a single sample.
        sample = record.samples[sample_index]

        # If not called, then nothing to do.
        if not sample.called:
            return False

        # Get the reference and alternate for this call.
        # NOTE: We reduce generality from what pyvcf since we are dealing
        # with a single sample.
        phase_char = sample.gt_phase_char()
        alts = sample.gt_bases.split(phase_char)

        # TODO: Figure out proper way to handle homozygous vs heterozygous.
        # assert len(set(alts)) == 1, (
        #         "Error while processing %s.\n"
        #         "We presently only support homozygous calls" %
        #                 (str(record)))

        # HACK: For now, if we see a record that we should handle, we
        # assume that we should take the first alt that is different than
        # the ref.
        ref = record.REF
        for alt_candidate in alts:
            alt = alt_candidate
            if alt == ref:
                continue

        pythonic_position = record.POS - 1

        try:
            return self._update_genome_record_for_variant(
                    pythonic_position, ref, alt, metadata=metadata)
        except AssertionError as e:
            raise AssertionError(
                    "AssertionError while fixing record %s\n%s" % (
                            str(record), str(e)))


    def _update_genome_record_for_variant(self, pythonic_position, ref, alt,
            note=None, metadata=None):
        """Updates self.genome_record with the passed in data.

        Logic extracted into own method for testing.
        """
        ref = ref.upper()
        alt = alt.upper()

        if ref == alt:
            # Nothing to do.
            return False

        # First, check whether the genome already looks like what it would
        # after the variant is fixed. We do this by fake-removing the ref
        # and adding the alt.
        # NOTE: One specific case that causes this issue is when we are
        # manually making the TAG changes before applying the rest of the
        # variants called in the VCF. This may bite us down the road, but this
        # whole script should probably be re-written in the near future.
        first_base_position = (
                self.runtime_liftover.convert_source_position_to_target(
                        pythonic_position, or_next=True))
        fake_seq = str(
                self.genome_record.seq[:first_base_position] +
                alt +
                self.genome_record.seq[first_base_position + len(ref):])
        if (fake_seq == str(self.genome_record.seq)):
            # Nothing to do.
            return False

        # The reason we don't just switch out the sequence above is that
        # we need to get the annotations right.
        # NOTE: Or is the above actually a more elegant way to do what follows?

        # Now determine the kind of mutation.
        if _is_snp(ref, alt):
            return self.handle_snp({
                'position': pythonic_position,
                'ref': ref,
                'alt': alt
            }, metadata=metadata)

        elif _is_deletion(ref, alt):
            deleted_subseq = _get_deletion(ref, alt)
            deleted_subseq_index_start = ref.rindex(deleted_subseq)
            assert len(deleted_subseq) == len(ref) - deleted_subseq_index_start
            deleted_subseq_start = (pythonic_position +
                    deleted_subseq_index_start)
            return self.handle_deletion({
                'interval': (deleted_subseq_start,
                        deleted_subseq_start + len(deleted_subseq)),
                'validation_seq': deleted_subseq
            }, note=note, metadata=metadata)

        elif _is_insertion(ref, alt):
            insertion_seq = _get_insertion(ref, alt)
            alt_insertion_start_index = alt.rindex(insertion_seq)
            assert len(insertion_seq) == len(alt) - alt_insertion_start_index, (
                    "Error handling insertion: ref: %s, alt: %s, position: d" %
                    (ref, alt, pythonic_position))
            insertion_start = pythonic_position + alt_insertion_start_index
            return self.handle_insertion({
                'position': insertion_start,
                'sequence': insertion_seq
            }, note=note, metadata=metadata)

        else:
            # Since we can't exactly tell, just delete ref and insert alt.
            validation_seq = str(self.original_genome_record.seq[
                        pythonic_position:pythonic_position + len(ref)])
            self.handle_deletion({
                'interval': (pythonic_position,
                        pythonic_position + len(ref)),
                'validation_seq': validation_seq
            }, add_annotation=False)

            self.handle_insertion({
                'position': pythonic_position,
                'sequence': alt
            }, add_annotation=False)

            #### Calculate data for the annotation.

            # The source interval is the interval that was removed
            # from the source.
            source_interval = (pythonic_position, pythonic_position + len(ref))

            # The insertion is not mapped in the liftover after the insertion,
            # so grab the starting position.
            target_genome_start = (
                    self.runtime_liftover.convert_source_position_to_target(
                            pythonic_position - 1, or_next=True) + 1)

            feature_id = 'misc_variant_source_%d-%d' % source_interval
            feature_location = FeatureLocation(
                    target_genome_start,
                    target_genome_start + len(alt))
            feature = SeqFeature(
                    type=VARIANT_ANNOTATION_TYPE,
                    location=feature_location,
                    strand=1,
                    id=feature_id
            )
            if len(ref) <= MAX_REPLACE_CHARS:
                feature.qualifiers['replace'] = ref.lower()
            else:
                feature.qualifiers['replace'] = '%d base replacement' % len(ref)
            if note:
                feature.qualifiers['note'] = note
            if metadata:
                for key, value in metadata.iteritems():
                    if value:
                        feature.qualifiers[key] = value
            add_feature_to_seq_record(self.genome_record, feature)
            return True


    def handle_snp(self, variant_data, add_annotation=True, note=None,
            metadata=None):
        """Handle a single nucleotide position change.
        """
        source_snp_position = variant_data['position']
        ref_base = variant_data['ref']
        alt_base = variant_data['alt']
        snp_size = 1

        # First, translate the position to the frame of the updated genome.
        snp_position = (
                self.runtime_liftover.convert_source_position_to_target(
                        source_snp_position))

        if not snp_position:
            # Nothing to do. This exact position has probably been deleted.
            return False

        # Make sure the ref is what is expected. This is a non-thorough
        # but reasonable and bug check.
        assert ref_base == self.genome_record.seq[snp_position], (
                "Error fixing SNP at "
                "source position %d, "
                "target position %d, "
                "Expected: %s, observed: %s" % (
                        source_snp_position, snp_position, ref_base,
                        self.genome_record.seq[snp_position]))

        new_seq = (
                self.genome_record.seq[:snp_position] +
                alt_base +
                self.genome_record.seq[snp_position + 1:])
        self.genome_record.seq = new_seq

        if add_annotation:
            # Add feature marking SNP.
            snp_feature_location = FeatureLocation(
                    snp_position, snp_position + snp_size)
            snp_feature_id = 'snp_source_%d_%s_to_%s' % (
                    source_snp_position, ref_base, alt_base)
            snp_feature = SeqFeature(
                    type=VARIANT_ANNOTATION_TYPE,
                    location=snp_feature_location,
                    strand=1,
                    id=snp_feature_id
            )
            snp_feature.qualifiers['replace'] = ref_base.lower()
            if note:
                snp_feature.qualifiers['note'] = note
            if metadata:
                for key, value in metadata.iteritems():
                    if value:
                        snp_feature.qualifiers[key] = value
            add_feature_to_seq_record(self.genome_record, snp_feature)

        # Change as made.
        return True


    def handle_insertion(self, variant_data, add_annotation=True, note=None,
            metadata=None):
        """Handles an insertion at the position relative to the original
        genome.

        Args:
            variant_data: Dictionary with keys:
                * position: Pythonic position for the insertion relative
                    to the original genome record.
                * sequence: The sequence being inserted. One or more bases.
        """
        source_position = variant_data['position']
        seq = variant_data['sequence']

        # First, translate the position to the frame of the updated genome.
        target_genome_position = (
                self.runtime_liftover.convert_source_position_to_target(
                        source_position, or_next=True))

        # Insert the sequence at the provided position.
        insert_sequence_and_update_features(self.genome_record, seq,
                target_genome_position, extend_feature_ends=True)

        # Update the liftover interval mapping.
        self.runtime_liftover.handle_insertion(variant_data)

        if add_annotation:
            # Add a feature annotating the insertion.
            feature_id = 'insertion_source_%s' % (source_position,)
            feature_location = FeatureLocation(target_genome_position,
                    target_genome_position + len(seq))
            feature = SeqFeature(
                    type=VARIANT_ANNOTATION_TYPE,
                    location=feature_location,
                    strand=1,
                    id=feature_id
            )
            # TODO: This doesn't work with the .tbl format.
            # Figure out how to fix this.
            # feature.qualifiers['replace'] = ''
            if note:
                feature.qualifiers['note'] = note
            if metadata:
                for key, value in metadata.iteritems():
                    if value:
                        feature.qualifiers[key] = value
            add_feature_to_seq_record(self.genome_record, feature)

        # Change as made.
        return True


    def handle_deletion(self, variant_data, add_annotation=True, note=None,
            metadata=None):
        """Handles a deletion.

        After this operation, the genome_record reflects the deletion.

        Args:
            variant_data: Dictionary with keys:
                * interval: A two-tuple representing pythonic interval for the
                    deletion, i.e. (inclusive_start, exclusive_end).
                    e.g. (100, 102) is a deletion of the 2 bases at positions
                    100 and 101.
                * validation_seq: If provided, used to validate that the
                    interval being deleted is this sequence.
        """
        interval = variant_data['interval']

        # Inclusive-bounds interval for the target.
        target_genome_interval = [
                    self.runtime_liftover.convert_source_position_to_target(
                            bound, or_next=True)
                    for bound in (interval[0], interval[1] - 1)]

        assert (bases_in_interval(target_genome_interval) ==
                interval[1] - interval[0])

        target_genome_interval_pythonic = (
                target_genome_interval[0],
                target_genome_interval[1] + 1)

        delete_interval(self.genome_record, target_genome_interval_pythonic,
                validation_seq=variant_data.get('validation_seq', None))

        # Update the liftover mapping.
        self.runtime_liftover.handle_deletion({
            'interval': interval
        })

        if add_annotation:
            # Add a feature annotating the deletion.

            # Calculate the target genome interval for the annotation.
            # Annotate from the position before the deletion to the position
            # after.

            target_genome_interval_after_deletion = [
                    self.runtime_liftover.convert_source_position_to_target(
                            bound, or_next=True)
                    for bound in (interval[0] - 1, interval[1])]
            feature_id = 'deletion_source_%d-%d' % (
                    interval[0], interval[1])
            feature_location = FeatureLocation(
                    target_genome_interval_after_deletion[0],
                    target_genome_interval_after_deletion[1])
            feature = SeqFeature(
                    type=VARIANT_ANNOTATION_TYPE,
                    location=feature_location,
                    strand=1,
                    id=feature_id
            )
            ref = variant_data.get('validation_seq', '')
            if len(ref) <= MAX_REPLACE_CHARS:
                feature.qualifiers['replace'] = ref.lower()
            else:
                feature.qualifiers['replace'] = '%d base deletion' % len(ref)
            feature.qualifiers['source_deletion_interval'] = str(interval)
            if note:
                feature.qualifiers['note'] = note
            if metadata:
                for key, value in metadata.iteritems():
                    if value:
                        feature.qualifiers[key] = value
            add_feature_to_seq_record(self.genome_record, feature)

        # Change as made.
        return True



###############################################################################
# Main procedure entrypoint
###############################################################################

def run(original_genbank_path, output_root, vcf_path, sample_id,
        **kwargs):
    """Creates a modified genbank file starting from the original genbank
    and applying the changes indicated in the vcf file.

    Args:
        original_genbank_path: Path to the original genbank file.
        output_root: Root of filename, without extension. The extension
            will be appended to the name depending on the output type.
        vcf_path: Path to the vcf file.
        sample_id: Id of the targete in the vcf, e.g. recoli_misq_c31_321D.
        kwargs: Optional keyword args. Supported keys:
            * liftover_pickle_dest: The output file to write the pickled
                liftover interval mapping to.
            * output_format: One of 'fasta' or 'genbank'. Defaults to
                'genbank'.
            * variant_data_csv: If included, will use
                update_from_variant_data_csv() rather than standard run.
            * verbose: Toggle for amount of informative print statements during
                processing.

        Returns:
            The final SeqRecord that was also written to output.
    """

    # Strategy:
    # Iterate through the calls in the VCF file and incrementally
    # update the genome record. There are tricky subtleties including:
    #     * The frame of the target genome is constantly changing.

    # Nuances: When adding insertions/deletions, this may shift the overall
    # frame of the genome downstream from that particular position. We need a
    # liftover-like intermediate representation that allows us
    # to keep track of these accumulated shifts. For example, every successive
    # change that we want to make should have its position updated using
    # this method. That way, the annotation can potentially preserve the
    # position of the SNP relative to the original record, but we can
    # introduce the changes into the underlying sequence and update all
    # features appropriately.

    if isinstance(original_genbank_path, SeqRecord):
        genome_record = original_genbank_path
    else:
        # Read in the original genome.
        genome_record = SeqIO.read(original_genbank_path, 'genbank')

    # Get optional manual updates file.
    if 'manual_updates_filepath' in kwargs:
        manual_updates_filepath = kwargs['manual_updates_filepath']
    else:
        manual_updates_filepath = None

    # Create the object that encapsulates most of the calculation.
    vcf_to_genbank_maker = VCFToGenbankMaker(genome_record, vcf_path,
            sample_id, manual_updates_filepath)

    if 'variant_data_csv' in kwargs:
        vcf_to_genbank_maker.update_from_variant_data_csv(
                kwargs['variant_data_csv'])
    else:
        vcf_to_genbank_maker.run(
                verbose=kwargs.get('verbose', False),
                log_file=kwargs.get('log_file', None))

    # Write the final result.
    DEFAULT_OUTPUT_FORMAT = 'genbank'
    output_format = kwargs.get('output_format', DEFAULT_OUTPUT_FORMAT)
    output_path = output_root + '.' + output_format
    SeqIO.write(genome_record, output_path, output_format)

    # Optional: Pickle the liftover interval mappings.
    if 'liftover_pickle_dest' in kwargs:
        vcf_to_genbank_maker.runtime_liftover.pickle_interval_mapping(
                kwargs['liftover_pickle_dest'])

    return genome_record


###############################################################################
# Helpers to evaluate SNP type.
###############################################################################

def _is_snp(ref, alt):
    return len(ref) == 1 and alt in ['A', 'T', 'G', 'C']


def _is_deletion(ref, alt):
    return _get_deletion(ref, alt) is not None


def _get_deletion(ref, alt):
    """Extracts the portion of ref that is deleted relative to alt.

    Returns None if no valid deletion found.
    """
    if len(ref) <= len(alt):
        return None

    if len(alt) == 0:
        return ref

    # Make sure they are both uppercase for matching procedure below.
    ref = ref.upper()
    alt = alt.upper()

    # Step through the two simultaneously until the first mismatch.
    idx = 0
    while idx < len(alt):
        if ref[idx] != alt[idx]:
            break
        idx += 1

    if idx < len(alt):
        # Our definition of deletion requirex the entire alt to be matched.
        return None

    deletion = ref[idx:]
    if not deletion:
        return None

    return deletion


def _is_insertion(ref, alt):
    return _get_insertion(ref, alt) is not None


def _get_insertion(ref, alt):
    """Extracts the portion of alt that inserted relative to alt.
    """
    # Just call _get_deletion with params reversed.
    return _get_deletion(alt, ref)


###############################################################################
# Other utility methods
###############################################################################

def create_filtered_vcf(vcf_path, out_vcf_path, csv_with_pos_to_keep):
    """Filters the passed in vcf down to the variant calls that we actually
    want to add to the updated genbank file.

    Writes the results to out_vcf_path.

    The reason for this method that cleans up the vcf, rather than just
    using the csv directly is that the logic for going from vcf to genbank
    will hopefully be re-usable, so we might as take a first stab at it here.
    """
    # Positions uniquely identify SNPs (manually confirmed). The provided
    # csv file should only have SNPs that we are keeping.
    positions_to_keep = set([])
    with open(csv_with_pos_to_keep) as csv_fh:
        csv_reader = csv.DictReader(csv_fh)
        for row in csv_reader:
            positions_to_keep.add(int(row['POS']))

    # Now create a filtered a vcf with only the above positions.
    with open(vcf_path) as vcf_fh, open(out_vcf_path, 'w') as out_vcf_fh:
        vcf_reader = vcf.Reader(vcf_fh)
        vcf_writer = vcf.Writer(out_vcf_fh, vcf_reader)
        for record in vcf_reader:
            if record.POS in positions_to_keep:
                vcf_writer.write_record(record)


def add_TAG_annotations(genome_record):
    """Temporary method for adding our UAG mutations manually.

    Mutates the passed in genome_record by adding features for the
    amber SNPs.
    """
    TAG_ANNOTATION_TYPE = VARIANT_ANNOTATION_TYPE

    UAG_LOCATIONS_FILE = os.path.join(
            GENOMES_DIR, 'mg1655', 'mg1655_uag_locations.csv')

    # Import the list of UAG locations and make the positions 0-indexed
    # to be consistent with the BioPython convention.
    uag_location_list = []
    with open(UAG_LOCATIONS_FILE) as fh:
        fh.readline() # Drop the first line
        for line in fh.readlines():
            uag_location_list.append(int(line.strip()) - 1)

    uag_location_list = sorted(uag_location_list)
    for current_uag_position in uag_location_list:
        current_base = genome_record.seq[current_uag_position]
        if current_base == 'G':
            alt_base = 'A'
            feature_location = FeatureLocation(
                    current_uag_position - 2,
                    current_uag_position + 1)
            feature_strand = 1
        elif current_base == 'C':
            alt_base = 'T'
            feature_location = FeatureLocation(
                    current_uag_position,
                    current_uag_position + 3)
            feature_strand = -1
        else:
            raise AssertionError("Invalid base at position %d: %s" % (
                    current_uag_position, current_base))

        # Update the sequence.
        new_seq = (
            genome_record.seq[:current_uag_position] +
            alt_base +
            genome_record.seq[current_uag_position + 1:])
        genome_record.seq = new_seq

        # Add a feature annotation.
        feature_id = 'remove_uag_%d' % current_uag_position
        feature = SeqFeature(
                type=TAG_ANNOTATION_TYPE,
                location=feature_location,
                strand=feature_strand,
                id=feature_id
        )
        feature.qualifiers['replace'] = 'tag'
        feature.qualifiers['note'] = 'Reassigning UAG'
        add_feature_to_seq_record(genome_record, feature)


def get_vcf_metadata():
    """Returns a dictionary with keys being SNP positions and values
    being the data from DBG's analysis with SnpEFF etc.
    """
    VARIANT_METADATA = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_snps.csv')

    position_to_data_map = {}
    with open(VARIANT_METADATA) as csv_fh:
        csv_reader = csv.DictReader(csv_fh)
        for row in csv_reader:
            position_to_data_map[int(row['POS'])] = {
                    'EFF_AA': row['EFF_AA'],
                    'EFF_CODON': row['EFF_CODON'],
                    'EFF_FUNC': row['EFF_FUNC'],
                    'EFF_SEV': row['EFF_SEV'],
            }

    return position_to_data_map


###############################################################################
# Testing and scripting
###############################################################################

if __name__ == '__main__':
    MG1655_GENBANK = os.path.join(GENOMES_DIR, 'mg1655', 'mg1655.genbank')
    RECOLI_ALL_SNPS = os.path.join(REC1_C321D_ROOT, 'recoli_all_snps.vcf')
    C321D_SNPS_VCF = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_snps.vcf')
    REC1_C321D_SAMPLE_ID = 'recoli_misq_c31_321D'
    REC1_C321D_OUTPUT_ROOT = os.path.join(REC1_C321D_ROOT, 'rec1_c321d')
    REC1_C321D_PRELIM_OUTPUT_ROOT = os.path.join(REC1_C321D_ROOT,
            'rec1_c321d.preliminary')
    PICKLE_DEST = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_liftover.pickle')
    MANUAL_UPDATES = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_manual_fixes.txt')

    # # First filter the edits that we actually want to keep.
    # # TODO: Make the run() method take as an optional argument a user-specified
    # # argument for filtering which vcf rows are incorporated into the Genbank.
    # CSV_WITH_POS_TO_KEEP = os.path.join(REC1_C321D_ROOT, 'rec1_c321d_snps.csv')
    # create_filtered_vcf(RECOLI_ALL_SNPS, C321D_SNPS_VCF, CSV_WITH_POS_TO_KEEP)

    # Now run the genbank creator.
    kwargs = {
        'liftover_pickle_dest': PICKLE_DEST,
        'manual_updates_filepath': MANUAL_UPDATES,
        'output_format': 'genbank',
        'verbose': True
    }
    run(MG1655_GENBANK, REC1_C321D_OUTPUT_ROOT, C321D_SNPS_VCF,
            REC1_C321D_SAMPLE_ID, **kwargs)


    ### Developing manual updates

    # REC1_C321D_MANUAL_ONLY_OUT = os.path.join(REC1_C321D_ROOT,
    #         'rec1_c321d_manual_only.genbank')
    # genome_record = SeqIO.read(MG1655_GENBANK, 'genbank')
    # vcf_to_genbank_maker = VCFToGenbankMaker(genome_record, None,
    #         None, MANUAL_UPDATES)
    # vcf_to_genbank_maker.handle_manual_updates()
    # SeqIO.write(genome_record, REC1_C321D_MANUAL_ONLY_OUT, 'genbank')


    ### Developing chain file output.

    # runtime_liftover = RuntimeLiftover(genome_record)
    # runtime_liftover.pickle_interval_mapping(PICKLE_DEST)
    # print runtime_liftover._interval_mapping

    # CHAIN_FILE_DEST = os.path.join(GENOMES_DIR, 'mg1655.to.rec1_c321d.chain')
    # runtime_2 = RuntimeLiftover.from_pickled_intervals(
    #         [], PICKLE_DEST)
    # runtime_2.write_chain_file(CHAIN_FILE_DEST)


    ### Generating source genome positions for c321D Genbank SNPs

    # OUT_CSV_FIELD_NAMES = ['position', 'ref', 'alt', 'mg1655_position']
    # REALIGNED_SNPS_IN = os.path.join(REC1_C321D_ROOT, 'realignment_snps.csv')
    # REALIGNED_SNPS_WITH_MG1655_POS = os.path.join(
    #         REC1_C321D_ROOT, 'realignment_snps_with_mg1655_positions.csv')
    # MG1655_GENOME_RECORD = SeqIO.read(MG1655_GENBANK, 'genbank')
    # runtime_liftover_obj = RuntimeLiftover.from_pickled_intervals(
    #         MG1655_GENOME_RECORD, PICKLE_DEST)
    # with open(REALIGNED_SNPS_WITH_MG1655_POS, 'w') as out_csv_fh:
    #     csv_writer = csv.DictWriter(out_csv_fh, OUT_CSV_FIELD_NAMES)
    #     csv_writer.writeheader()
    #     with open(REALIGNED_SNPS_IN) as in_csv_fh:
    #         csv_reader = csv.DictReader(in_csv_fh)
    #         for row in csv_reader:
    #             out_row = copy.copy(row)

    #             # Convert the position. Note the indexing correction.
    #             c321d_position_pythonic = int(row['position']) - 1
    #             mg1655_position_pythonic = (
    #                     runtime_liftover_obj.convert_target_position_to_source(
    #                             c321d_position_pythonic))
    #             if mg1655_position_pythonic is not None:
    #                 # Validation check that the base at that position is
    #                 # correct.
    #                 assert row['ref'][0] == (
    #                         runtime_liftover_obj.source_genome_record.seq[
    #                                 mg1655_position_pythonic])

    #                 # Make the position one-indexed to be consistent
    #                 # with the vcf standard.
    #                 mg1655_position = mg1655_position_pythonic + 1
    #             else:
    #                 mg1655_position = None

    #             out_row['mg1655_position'] = mg1655_position
    #             csv_writer.writerow(out_row)
