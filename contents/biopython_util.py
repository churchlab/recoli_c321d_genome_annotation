"""
Methods that encapsulate genome refactoring techniques applied to BioPython
objects.
"""

import csv
import copy
import re
from uuid import uuid4

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.Seq import translate
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from annotation_feature_types import FIX_OVERLAP_SIZE_4
from annotation_feature_types import FIX_OVERLAP_SYNONYMOUS
from annotation_feature_types import FORBIDDEN_CODON_FEATURE_TYPE


# When performing an insertion, we add a feature to the SeqRecord
# to track this insertion and give it some semantically-meaningful
# type depending on the context of the insertion.
class InsertType(object):
    FIX_OVERLAP_HEAD_COPY = 'overlap_head_cp'
    FIX_OVERLAP_RBS_COPY = 'rbs_cp'
    PARTITION_GEN9_SEG = 'gen9_seg'
    SYNTHESIS_SEG = 'synthesis_seg'


COMMONLY_IGNORED_FEATURE_TYPES = (
        ['source', 'gene', 'misc_feature'] +
        [FORBIDDEN_CODON_FEATURE_TYPE, FIX_OVERLAP_SYNONYMOUS, FIX_OVERLAP_SIZE_4] +
        [getattr(InsertType, name)
                for name in dir(InsertType)
                if not re.match('__', name)]
)


# Additional keys used to save metadata in the SeqFeature objects
# that will then be written to the genbank file.
class FeatureQualifierExtensionKeys(object):
    GENOME_REFACTOR_ID = 'genome_refactor_id'


def get_genome_record(
        source_gbk_file,
        features_to_ignore=['source', 'gene'],
        only_get_features=[],
        use_old_id_strategy=False,
        assert_unique_ids=True,
        verbose_ids=False,
        drop_duplicate_features=False):
    """Reads the source .gbk file and returns a Bio.SeqRecord.SeqRecord
    object generated from the source file.

    Adds ids to each feature, making sure they are unique.

    Args:
        source_gbk_file: The source Genbank file.
        features_to_ignore: Genbank is a noisy format, so we specify some
            types of features that we usually want to ignore.
        only_get_features: If specified, then only get those features. If this
            parameter is passed, then features_to_ignore is ignored.
        use_old_id_strategy: Whether to use the old strategy for calculating
            the SeqFeature.id.
        assert_unique_ids: Whether to make sure we have unique ids. In some
            use-cases we can safely ignore this.

    Returns:
        A SeqRecord object.
    """
    with open(source_gbk_file) as fh:
        record = SeqIO.read(fh, 'genbank')

    # Add ids to the features.
    new_features = []
    ids = set()
    for feature in record.features:
        if only_get_features:
            if not feature.type in only_get_features:
                continue
        # Skip features that we are fine with just dropping.
        elif feature.type in features_to_ignore:
            continue

        # Create an id for each feature.
        if use_old_id_strategy:
            new_id = _make_feature_id_old_strategy(feature)
        else:
            new_id = _make_feature_id(feature, verbose_ids)

        # Make sure ids are unique.
        if assert_unique_ids:
            if drop_duplicate_features and new_id in ids:
                continue # Drop this feature
            # HACK: Not sure why forbiddenCodon features are doubled.
            if re.match(r'forbiddenCodon', new_id) and new_id in ids:
                new_id = new_id + '_' + str(uuid4())[:8]
            else:
                assert new_id not in ids, "Id %s already exists." % new_id
            ids.add(new_id)

        # Copy the feature and update its id.
        new_feature = copy.deepcopy(feature)
        new_feature.id = new_id
        new_feature.qualifiers[
                FeatureQualifierExtensionKeys.GENOME_REFACTOR_ID] = new_id
        new_features.append(new_feature)

    # Replace the features.
    record.features = new_features

    # Wrap in an object with our custom methods (e.g. support for
    # internal dictionary mapping from feature id to position in the
    # underlying SeqRecord.features list.
    # add_custom_seq_record_methods(record)

    # Return the final SeqRecord.
    return record


def _make_feature_id(feature, verbose_ids=False):
    """Creates an id for SeqFeature object. Ideally this should create
    something unique. Clients of this method are responsible for enforcing
    uniqueness, which may require using a method with different logic,
    or extending this method. This is currently sufficient for working with the
    MDS42 genbank file downloaded from NCBI, but may need to be tweaked for
    other genbank files.

    Args:
        feature: A SeqFeature object read from, say, a genbank file.

    Return:
        A string id for the seq_feature.
    """
    if (FeatureQualifierExtensionKeys.GENOME_REFACTOR_ID in
            feature.qualifiers):
        new_id = feature.qualifiers[
                FeatureQualifierExtensionKeys.GENOME_REFACTOR_ID][0]
    elif 'locus_tag' in feature.qualifiers:
        if verbose_ids:
            # NOTE: Locus tags should be unique for a well-formatted genbank file.
            # UPDATE: Not true for mg1655.genbank.
            new_id = (feature.type + '_' +
                    feature.qualifiers.get('locus_tag', ['UNKNOWN_LOCUS'])[0])
        else:
            new_id = feature.qualifiers.get('locus_tag', ['UNKNOWN_LOCUS'])[0]
    else:
        new_id = _make_feature_id_old_strategy(feature)

    return new_id


def _make_feature_id_old_strategy(feature):
    """Exhaustive id method that all but guarantees unique id.

    May be necessary for analysis of data from before these changes.
    """
    new_id = (
            feature.type + '_' +
            feature.qualifiers.get('gene', 'x')[0] + '_' +
            feature.qualifiers.get('locus_tag', 'x')[0] + '_' +
            str(feature.location.start + 1) + '_' +
            str(feature.location.end + 1)
    )
    return new_id


def shift_seq_feature_right_with_head_copy(
        seq_feature, seq_record, new_start_position):
    """Moves a SeqFeature to a new position to the right of where it was
    previously, copying the the number of bases at the
    head of the seq_feature according to the distance moved and inserting
    them in the gap that would be left behind otherwise.

    This is used, e.g. when separating overlapping features.

    Original sequence:
        N N N N A T G
                A T G N N N N N

    New sequence:
        N N N N A T G A T G N N N N N

    NOTE: seq_record is mutated by this method.
    """
    shift_distance = new_start_position - seq_feature.location.start
    assert shift_distance > 0, "Shift must be to the right."

    # NOTE: Manually extract to make sure we get the forward read in case
    # the feature is on the negative strand.
    head_copy_seq = seq_record.seq[
            seq_feature.location.start:seq_feature.location.start +
            shift_distance]

    seq_record = insert_sequence(
            seq_record,
            head_copy_seq,
            new_start_position,
            safe_features=[seq_feature],
            insert_feature_type=InsertType.FIX_OVERLAP_HEAD_COPY,
            insert_feature_id=(str(seq_feature.id) + '_' +
                    InsertType.FIX_OVERLAP_HEAD_COPY),
            insert_feature_strand=seq_feature.strand
    )

    # Replace the shifted feature.
    shifted_seq_feature = get_feature_by_id(seq_record, seq_feature.id)
    shifted_seq_feature = shifted_seq_feature._shift(shift_distance)
    replace_feature(seq_record, seq_feature.id, shifted_seq_feature)

    return seq_record


def insert_sequence(
        seq_record,
        insert_seq,
        insert_start_position,
        safe_features=[],
        insert_feature_type=None,
        insert_feature_id=None,
        insert_feature_strand=1):
    """WARNING: You probably want to user insert_sequence_and_update_features(),
    unless you know what you are doing. This method was created for complicated
    reasons like splitting overlapping genes and so is not really a clean
    method.

    Insert a sequence into seq_record at the given start_position,
    updating all child SeqFeature positions appropriately. Specifically, all
    SeqFeatures to the right of the insertion need to have their positions
    adjusted. If the insertion is in the middle of a feature, the feature
    must be declared in safe_features, otherwise this method will throw
    an error. It is then up to the caller of this method to handle the
    affected features manually.

    NOTE: seq_record is mutated by this method.

    TODO: Be more careful about inserting in the middle of a feature.

    Args:
        seq_record: The mutable SeqRecord object.
        insert_seq: The sequence to insert, as read along the forward strand.
        insert_start_position: The start position of the insertion.
        safe_features: Features that we are consciously inserting over.
            This is usually the feature's head that is being copied.
        insert_feature_type: Semantically-informative type for the insertion used
            to create a SeqFeature that is added to the SeqRecord in order
            to maintain a record of the insertion. If this is omitted,
            no feature will be created.
        insert_feature_id: Unique id for the feature to be added to represent
            the insertion.
        insert_feature_strand: The strand the inserted feature is one.
            NOTE: The insert_seq should still be in the direction of the
            forward (+1) strand.

    Returns:
        The mutated SeqRecord.
    """
    # Basic validation.
    error_msg = ("Either both or neither of insert_feature_type and "
            "insert_feature_id must be declared.")
    assert bool(insert_feature_type) == bool(insert_feature_id), error_msg

    # Check that the insertion does not overlap any unexpected features.
    safe_feature_ids = map(lambda feature: feature.id, safe_features)
    for feature in seq_record.features:
        if (not feature.id in safe_feature_ids and
                feature.location.start < insert_start_position and
                insert_start_position < feature.location.end):
            raise ValueError(
                    "Unexpected insertion starting at %s inside of feature %s" %
                    (str(insert_start_position), str(feature)))

    insert_size = len(insert_seq)
    new_seq = (
            seq_record.seq[:insert_start_position] +
            insert_seq +
            seq_record.seq[insert_start_position:]
    )
    assert len(new_seq) == len(seq_record.seq) + insert_size
    seq_record.seq = new_seq

    # Update the features to the right of the insertion.
    updated_features = []
    for feature in seq_record.features:
        if feature.location.start >= insert_start_position:
            feature = feature._shift(insert_size)
        updated_features.append(feature)
    seq_record.features = updated_features

    # If insert_feature_type is defined, add a new feature to represent the inserted
    # sequence.
    if insert_feature_type and insert_feature_id:
        new_seq_feature = SeqFeature(
                location=FeatureLocation(insert_start_position,
                        insert_start_position + insert_size),
                type=insert_feature_type,
                strand=insert_feature_strand,
                id=insert_feature_id
        )
        new_seq_feature.qualifiers[
                FeatureQualifierExtensionKeys.GENOME_REFACTOR_ID] = (
                        insert_feature_id)
        add_feature_to_seq_record(seq_record, new_seq_feature)

    # Return the mutated SeqRecord.
    return seq_record


def insert_sequence_and_update_features(seq_record, insert_seq,
        insert_start_position, extend_feature_ends=False,
        no_feature_shift_if_inside=False,
        insert_feature_type=None,
        insert_feature_id=None,
        insert_feature_strand=1):
    """Inserts a sequence at the given position and updates feature locations.

    This method differs from insert_sequence() in that this method doesn't
    "know" what it's being used for.

    Features are updated according to the following rules:
        * If a feature lies upstream of an insertion, nothing to do.
        * If the insertion is inside of a feature, increase the size of the
            feature the appropriate amount.
            * If the feature has sub_features, then recursively apply
                these rules.

    Args:
        seq_record: The mutable SeqRecord object.
        insert_seq: The sequence to insert, as read along the forward strand.
        insert_start_position: The start position of the insertion. Pythonic.
        extend_feature_ends: If True, any features whose end position is at
            insert_start_position will have their end extended to capture
            the entire insertion.
        no_feature_shift_if_inside: Don't shift the feature even if the change
            is inside the feature.

    Returns:
        The mutated seq_record.
    """
    # Update the sequence.
    insert_size = len(insert_seq)
    new_seq = (
            seq_record.seq[:insert_start_position] +
            insert_seq +
            seq_record.seq[insert_start_position:]
    )
    assert len(new_seq) == len(seq_record.seq) + insert_size
    seq_record.seq = new_seq

    # Update the features.
    updated_features = []
    for feature in seq_record.features:
        if feature.location.start >= insert_start_position:
            # Feature is downstream of deletion. Move the whole thing.
            feature = feature._shift(insert_size)
        elif no_feature_shift_if_inside:
            pass
        elif ((feature.location.end > insert_start_position and
                    feature.location.start <= insert_start_position) or
                (extend_feature_ends and
                        feature.location.end == insert_start_position)):
            # Insertion is inside the feature.
            _shift_feature_with_event_at_position(
                    feature, insert_start_position, insert_size,
                    extend_feature_ends=extend_feature_ends)
        updated_features.append(feature)
    seq_record.features = updated_features

    # If insert_feature_type is defined, add a new feature to represent the inserted
    # sequence.
    if insert_feature_type and insert_feature_id:
        new_seq_feature = SeqFeature(
                location=FeatureLocation(insert_start_position,
                        insert_start_position + insert_size),
                type=insert_feature_type,
                strand=insert_feature_strand,
                id=insert_feature_id
        )
        new_seq_feature.qualifiers[
                FeatureQualifierExtensionKeys.GENOME_REFACTOR_ID] = (
                        insert_feature_id)
        add_feature_to_seq_record(seq_record, new_seq_feature)

def delete_single_base(seq_record, base_position, base_value=None):
    """Utility method for deleting a single base in a SeqRecord
    and adjusting all the feature positions accordingly.

    Args:
        seq_record: The SeqRecord that is mutated by this method.
        base_position: Position in the genome of the base to delete.
        base_value: Value of the base. Used as a sanity check.

    Returns:
        The instance of SeqRecord that was passed in and mutated.
    """
    # NOTE: Maybe make this more flexible when the case arises.
    SHIFT_AMOUNT = -1

    # Sanity check.
    if base_value is not None:
        assert base_value == seq_record.seq[base_position]

    # Change the underlying sequence.
    new_seq = (
            seq_record.seq[:base_position] +
            seq_record.seq[base_position + 1:]
    )
    seq_record.seq = new_seq

    # Update affected feature positions.
    updated_features = []
    for feature in seq_record.features:
        if feature.location.start >= base_position:
            # Feature is downstream of deletion. Move the whole thing.
            feature = feature._shift(SHIFT_AMOUNT)
        elif (feature.location.end > base_position and
                feature.location.start <= base_position):
            # Deletion is inside the feature.
            _shift_feature_with_event_at_position(
                    feature, base_position, SHIFT_AMOUNT)
        # NOTE: No change to feature upstream of deletion.
        updated_features.append(feature)
    seq_record.features = updated_features

    # Return the mutated seq_record.
    return seq_record


def delete_interval(seq_record, interval, validation_seq=None,
        validation_start_seq=None,
        validation_end_seq=None):
    """Delete all of the bases in the interval and update features
    appropriately. Mutates seq_record.

    When updating features:
        The feature's end position is shifted upstream for each deleted
        base that overlapped the feature.

    Args:
        seq_record: SeqRecord object that will be mutated to reflect changes.
        interval: Pythonic interval to delete (inclusive lower-bound, exclusive
            upper-bound). Assummes a 0-indexed genome.
        validation_seq: If provided, we validate that this is the sequence
            in the interval. For debugging purposes.

    Returns:
        The mutated SeqRecord.
    """
    if validation_seq:
        assert validation_seq == str(
                seq_record.seq[interval[0]:interval[1]]), (
                "\tValidation sequence doesn't match interval.\n"
                "\tInterval: %s\n"
                "\tActual seq (truncated): %s\n"
                "\tActual seq start +/- 5: %s\n"
                "\tValidation seq: %s" % (
                        str(interval),
                        str(seq_record.seq[interval[0]:
                                min(interval[0] + 6, interval[1])]),
                        str(seq_record.seq[interval[0] - 5:interval[0]]) + '*' +
                                str(seq_record.seq[interval[0]:
                                        interval[0] + 5]),
                        validation_seq
        ))

    if validation_start_seq:
        assert validation_start_seq == str(seq_record.seq[
                interval[0]:interval[0] + len(validation_start_seq)])

    if validation_end_seq:
        assert validation_end_seq == str(seq_record.seq[
                interval[1] - len(validation_end_seq):interval[1]])

    interval_size = interval[1] - interval[0]

    # Change the underlying sequence.
    new_seq = (
            seq_record.seq[:interval[0]] +
            seq_record.seq[interval[0] + interval_size:]
    )
    seq_record.seq = new_seq

    # Update the features.
    updated_features = []
    for feature in seq_record.features:
        if feature.location.start >= interval[1]:
            # Entire feature is downstream of deletion. Move the whole thing.
            feature = feature._shift(-1 * interval_size)
        elif does_interval_overlap_feature(interval, feature):
            # Feature overlaps interval.
            feature = _shift_feature_overlapped_by_event(feature, interval, -1)
        if feature:
            updated_features.append(feature)
    seq_record.features = updated_features

    # Return the mutated seq_record.
    return seq_record


def _shift_feature_with_event_at_position(
        feature, event_position, shift_amount, extend_feature_ends=False):
    """Mutates the feature by shifting its end the given amount.

    The Bio.SeqFeature._shift() method doesn't update the sub-features
    correctly which is probably why it is underscore-indicated-private.

    So the effect of this method is:
        1) For each sub-feature:
            * Shift the whole thing, if after the event.
            * Shift just the end, if overlapping the event.
            * No shift, if the sub_feature lies entirely before the event.
        2) For the whole feature:
            * Shift the end of the whole feature.
    """
    # Save original length for final validation.
    # NOTE: When a SeqFeature has sub_features, SeqFeature.__len__() returns
    # the sum of the lengths of the sub-features rather than the size of
    # the feature.location interval, so we are explicit here with the assertion
    # check.
    orig_feature_len = len(feature.location)

    ### 1) Update the sub_features appropriately.
    updated_sub_features = []
    for sub_feature in feature.sub_features:
        if sub_feature.location.start > event_position:
            sub_feature = sub_feature._shift(shift_amount)
        elif ((sub_feature.location.end > event_position and
                    sub_feature.location.start <= event_position) or
                (extend_feature_ends and
                        sub_feature.location.end == event_position)):
            _shift_end_of_feature(sub_feature, shift_amount)
        updated_sub_features.append(sub_feature)
    feature.sub_features = updated_sub_features

    ### 2) Shift the end of the feature.
    _shift_end_of_feature(feature, shift_amount)

    # Assert the operation succeeded.
    expected_updated_feature_len = orig_feature_len + shift_amount
    assert expected_updated_feature_len == len(feature.location), (
            "Feature: %s, Expected: %d, actual: %d" % (
                    feature,
                    expected_updated_feature_len,
                    len(feature)))


def _shift_end_of_feature(feature, shift_amount):
    """Mutates the feature by shifting the end of its location.
    """
    feature.location = FeatureLocation(
            start=feature.location.start,
            end=feature.location.end._shift(shift_amount),
            strand=feature.strand)


def _shift_feature_overlapped_by_event(feature, event_interval, event_polarity):
    """Shifts the feature and its sub_features depending on which parts
    of the interval overlap it.

    NOTE: Clients should use the returned feature.

    Args:
        feature: SeqFeature that may be mutated, but clients should
            use the return value.
        event_interval: The interval of the event in the frame of the
            SeqRecord to which the feature belongs.
        event_polarity: 1 for insertion, -1 for deletion.

    Returns:
        An updated copy of the feature, or None.
    """
    assert event_polarity == -1 or event_polarity == 1

    event_interval_size = event_interval[1] - event_interval[0]

    # Determine how much to shift the feature
    if (event_interval[0] <= feature.location.start and
            event_interval[1] >= feature.location.end):
        # Event covers entire feature. Remove it.
        #   >>>>>>>>>
        # (...........)
        return None
    elif (event_interval[0] > feature.location.start and
            event_interval[1] < feature.location.end):
        # Event is inside feature. We'll shift the end
        # >>>>>>>>>>
        #   (....)
        shift_amount = event_interval_size
    elif (event_interval[0] > feature.location.start and
            event_interval[1] >= feature.location.end):
        # >>>>>>>>>
        #   (........)
        shift_amount = feature.location.end - event_interval[0]
    else:
        #     >>>>>>>>>
        # (........)
        shift_amount = event_interval[1] - feature.location.start
    _shift_end_of_feature(feature, shift_amount * event_polarity)

    # Update the sub_features appropriately.
    updated_sub_features = []
    for sub_feature in feature.sub_features:
        if sub_feature.location.start >= event_interval[1]:
            # Entire feature is downstream of deletion. Move the whole thing.
            sub_feature = sub_feature._shift(event_polarity * event_interval_size)
        elif does_interval_overlap_feature(event_interval, sub_feature):
            # Feature overlaps interval.
            sub_feature = _shift_feature_overlapped_by_event(
                    sub_feature, event_interval, event_polarity)
        if sub_feature:
            updated_sub_features.append(sub_feature)
    feature.sub_features = updated_sub_features

    return feature


def add_feature_to_seq_record(seq_record, new_seq_feature, ignore_ids=False):
    """Adds the feature to the SeqRecord list of features and updates
    the custom mapping that we use.

    Args:
        seq_record: The SeqRecord to modify.
        new_seq_feature: The SeqFeature to add.
        ignore_ids: If True, does not check for unique ids. Only works when
            not using an underlying map.

    Raises:
        AssertionError if a feature with that id already exists.
    """
    if (hasattr(seq_record, 'feature_id_to_index_map') and ignore_ids):
        raise AssertionError(
                "Cannot ignore ids when using feature_id_to_index_map in "
                "SeqRecord.")

    # HACK: Assuming our feature ids are safe, then this feature does not
    # need to be added again.
    if (hasattr(seq_record, 'feature_id_to_index_map') and
            new_seq_feature.id in seq_record.feature_id_to_index_map):
        return

    not_unique_err = "Feature id already exists %s" % (new_seq_feature.id,)
    if hasattr(seq_record, 'feature_id_to_index_map'):
        assert new_seq_feature.id not in seq_record.feature_id_to_index_map,\
                not_unique_err
        seq_record.features.append(new_seq_feature)
        last_index = len(seq_record.features) - 1
        seq_record.feature_id_to_index_map[new_seq_feature.id] = last_index
        num_keys = len(seq_record.feature_id_to_index_map.keys())
        num_features = len(seq_record.features)
        assert num_keys == num_features, ("Num keys: %d.  Num features: %d" %
                (num_keys, num_features))
    else:
        if ignore_ids:
            seq_record.features.append(new_seq_feature)
            return

        # Otherwise check uniqueness.
        unique = True
        for feature in seq_record.features:
            if feature.id == new_seq_feature.id:
                unique = False
                break
        if unique:
            seq_record.features.append(new_seq_feature)
        # else:
        #     raise AssertionError(not_unique_err)


def swap_feature_codon_at_position(
        genome_record, feature_id, position, previous, new):
    """Changes the value of the sequence at the given position, mutating
    genome_record to reflect this change.

    Args:
        genome_record: SeqRecord.
        feature_id: Id of the feature in the SeqRecord.
        position: Position relative to the start of the feature.
        previous: Previous sequence at that position to remove.
        new: New sequence to insert at that position.

    Returns:
        Handle to genome_record, which has been mutated.
    """
    len_previous = len(previous)
    len_new = len(new)

    feature = get_feature_by_id(genome_record, feature_id)
    orig_feature_seq = feature.extract(genome_record.seq)
    assert previous == str(orig_feature_seq[position:position + len_previous]), (
            "Expected: %s, actual at position: %s" % (
                    previous,
                    str(orig_feature_seq[position-5:position]).lower() +
                            str(orig_feature_seq[position:
                                    position + len_previous]) +
                            str(orig_feature_seq[position + len_previous:
                                    position + len_previous + 5]).lower()))

    new_feature_seq = (
            orig_feature_seq[:position] +
            new +
            orig_feature_seq[position + len_previous:]
    )
    assert len(orig_feature_seq) + len_new - len_previous == len(new_feature_seq)

    return update_feature_seq(
            genome_record,
            feature_id,
            new_feature_seq,
            allow_translation_change=True)


def get_feature_by_id(seq_record, feature_id):
    """Returns the feature with the given id, or None if not found.
    """
    if hasattr(seq_record, 'feature_id_to_index_map'):
        feature_index = seq_record.feature_id_to_index_map[feature_id]
        feature = seq_record.features[feature_index]
        error_msg = "Bug in our custom SeqRecord feature id map."
        assert feature_id == feature.id, error_msg
        return feature
    else:
        for feature in seq_record.features:
            if feature.id == feature_id:
                return feature

    return None


def replace_feature(seq_record, feature_id, new_feature):
    """Replace the feature in the given seq_record having the given
    feature_id with the new feature.
    """
    if hasattr(seq_record, 'feature_id_to_index_map'):
        feature_index = seq_record.feature_id_to_index_map[feature_id]
        old_feature = seq_record.features[feature_index]
        assert old_feature.id == new_feature.id
        seq_record.features[feature_index] = new_feature
    else:
        new_features = []
        for feature in seq_record.features:
            if feature.id == feature_id:
                new_features.append(new_feature)
            else:
                new_features.append(feature)
        seq_record.features = new_features


def get_feature_gene(feature):
    """Returns the gene for the feature.

    NOTE: Assumes each feature only has one gene.
    """
    gene_list = feature.qualifiers.get('gene')
    if gene_list:
        return gene_list[0]
    return None


def add_custom_seq_record_methods(seq_record):
    """Adds custom methods to the SeqRecord object.

    NOTE: We have not confirmed that these extensions work with all endpoints
    so use these methods with caution.
    """
    # Populate the map that maintains a dict mapping from feature id
    # to index in the feature list of the SeqRecord.
    feature_id_to_index_map = {}
    for feature_index in range(len(seq_record.features)):
        feature = seq_record.features[feature_index]
        feature_id_to_index_map[feature.id] = feature_index
    setattr(seq_record, 'feature_id_to_index_map', feature_id_to_index_map)


def update_seq_record_feature(seq_record, feature_id, result_obj):
    """Modifies the feature with the given id, including swapping
    out the underlying sequence in the global sequence record.

    Args:
        seq_record: The mutable SeqRecord object.
        feature_id: The id of the feature to update.
        result_obj: Object returned from refactoring code, with keys:
            * feature_id
            * is_success
            * exception_string
            * orig_feature_seq
            * new_feature_seq
            * profiler_stats

    NOTE: result_obj['new_feature_seq'] should be the actual 5'-to-3' sequence.
    If the feature is transcribed on the negative strand, then this method will
    properly insert the reverse complement of new_seq into the global underlying
    sequence which is kept as a forward strand.

    Returns:
        A handle to the updated copy of seq_record.
    """
    feature = get_feature_by_id(seq_record, feature_id)

    # Replace the underlying sequence.
    new_seq = result_obj['new_feature_seq']
    update_feature_seq(seq_record, feature_id, new_seq)

    # Add the new features.
    if 'new_feature_list' in result_obj:
        for new_feature in result_obj['new_feature_list']:
            add_feature_to_seq_record(seq_record, new_feature)

    # Add metadata to feature.qualifiers. This dictionary is then written
    # as metadata to a genbank output file, which may then be useful while
    # viewing in something like Geneious.
    for profiler_name, stats in result_obj['profiler_stats'].iteritems():
        feature.qualifiers['error_tolerance_' + profiler_name] = (
                stats['error_tolerance'])

    return seq_record


def update_feature_seq(
        seq_record, feature_id, new_seq, allow_translation_change=False):
    """Updates the underlying sequence of the feature.

    The new sequence should be the actual 5'-to-3' sequence for the feature.
    This method will handle inserting into the seq_record in the correct
    polarity.

    Args:
        seq_record: Mutated SeqRecord object.
        feature_id: Id of the feature.
        new_seq: The new sequence.
        allow_translation_change: Allows the caller to acknowledge that the
            translation may change with the new sequence.
    """
    feature = get_feature_by_id(seq_record, feature_id)
    return update_feature_seq_given_feature(
        seq_record, feature, new_seq,
        allow_translation_change=allow_translation_change)


def update_feature_seq_given_feature(
        seq_record, feature, new_seq, allow_translation_change=False):
    """Delegate of update_feature_seq in act of (pseudo) method overloading.
    """
    if isinstance(new_seq, str):
        new_seq_obj = Seq(new_seq, generic_dna)
    else:
        new_seq_obj = new_seq

    # Assert that translation is the same.
    if not allow_translation_change:
        original_feature_seq = feature.extract(seq_record.seq)
        original_translation = translate_custom(str(original_feature_seq))
        refactored_translation = translate_custom(str(new_seq_obj))
        error_msg = "Translation mismatch for feature %s" % feature.id
        assert original_translation == refactored_translation, error_msg

    polarity_asserted_seq = copy.deepcopy(new_seq_obj)
    if feature.strand == -1:
        try:
            polarity_asserted_seq = polarity_asserted_seq.reverse_complement()
        except AttributeError:
            assert isinstance(polarity_asserted_seq, str)
            polarity_asserted_seq = reverse_complement(polarity_asserted_seq)

    updated_seq = (
            seq_record.seq[:feature.location.start] +
            polarity_asserted_seq +
            seq_record.seq[feature.location.start + len(new_seq):]
    )
    assert len(updated_seq) == len(seq_record.seq)
    seq_record.seq = updated_seq

    return seq_record


def update_feature_seq_at_position(seq_record, position, ref_seq, alt_seq):
    """Updates the feature seq at the given position from the ref to the alt.

    All sequences provided in the forward direction.
    """
    # Make sure the ref is what we think it is.
    actual_seq = str(seq_record.seq[position:position + len(ref_seq)]).upper()
    assert ref_seq.upper() == actual_seq, "Expected/Actual: \n" + ref_seq.upper() + "\n" + actual_seq


    # Make the change.
    updated_seq = (
            seq_record.seq[:position] +
            alt_seq +
            seq_record.seq[position + len(ref_seq):]
    )
    seq_record.seq = updated_seq
    return seq_record



def translate_custom(seq_str):
    """Custom translation method that accounts for the right species.

    For now we just use MDS42. Eventually we may want to allow different
    cofigs depending on the genetic code.

    Returns a string of the amino acid sequence.
    """
    return translate_mds42(seq_str)


def translate_mds42(seq_str):
    """Translation method for mds42 which uses the right translation table
    and handles exceptions appropriately.

    This method was created in response to needing to properly handle
    the start codon during translation checking.

    Returns a string of the amino acid sequence.
    """
    try:
        result = translate(seq_str, table=11, cds=True).strip('*')
    except TranslationError as e:
        # For now the cases where we care about the start codon and the
        # cases where we are okay with there being a stop codon in frame
        # are exclusive.
        # TODO: Can we do something less likely to bite us?
        if 'stop' in e.message:
            result = translate(seq_str)
        else:
            raise e
    return result


def show_all_gbk_annotation_keys(record_or_source):
    """Useful utility method for seeing all the annotation keys available
    in a gbk file.
    """
    if isinstance(record_or_source, SeqRecord):
        record = record_or_source
    else:
        with open(record_or_source) as fh:
            record = SeqIO.read(fh, 'genbank')

    counts = {}
    for feature in record.features:
        if feature.type in counts:
            counts[feature.type] += 1
        else:
            counts[feature.type] = 1

    return counts


def find_features_starting_at(
        position, genome_record, feature_type_restrict=[]):
    """Returns a list of features with the given start position.

    Args:
        position: Integer (or castable to Integer) position that is the
            start of feature we're looking for.
        genome_record: SeqRecord object to search.
        feature_type_restrict: Feature types to limit search and results to.

    Returns:
        A list of features ending at the requested position.
    """
    return _find_features_verbing_at(
            'start', position, genome_record, feature_type_restrict)


def find_features_ending_at(
        end_position, genome_record, feature_type_restrict=[]):
    """Returns a list of features with the given end position.

    Args:
        position: Integer (or castable to Integer) position that is the end
            of feature we're looking for.
        genome_record: SeqRecord object to search.
        feature_type_restrict: Feature types to limit search and results to.

    Returns:
        A list of features ending at the requested position.
    """
    return _find_features_verbing_at(
            'end', end_position, genome_record, feature_type_restrict)


VALID_FIND_FEATURE_VERBS = set(['start', 'end'])

def _find_features_verbing_at(
        verb, position, genome_record, feature_type_restrict=[]):
    """Returns a list of features with the location descriptor verb.

    Args:
        verb: Either 'start' or 'end'.
        position: Integer (or castable to Integer) position that is either
            the end or start of the features we're looking for.
        genome_record: SeqRecord object to search.
        feature_type_restrict: Feature types to limit search and results to.

    Returns:
        A list of features ending at the requested position.
    """
    assert verb in VALID_FIND_FEATURE_VERBS

    position = int(position)

    # Maybe filter the features we're searching.
    features_to_search = genome_record.features
    if feature_type_restrict:
        features_to_search = filter(
                lambda feature: feature.type in feature_type_restrict,
                features_to_search
        )

    # Loop through the remaining candidate features.
    results = []
    for feature in features_to_search:
        if position == int(feature.location.__getattribute__(verb)):
            results.append(feature)
    return results


def does_feature_overlap_other_feature(query_feature, other_feature):
    """Determines whether the features overlap.
    """
    return does_interval_overlap_feature(
            (query_feature.location.start, query_feature.location.end),
            other_feature)


def does_interval_overlap_feature(interval, feature):
    """Checks whether the given interval overlaps the feature's location.

    Args:
        interval: A two-tuple of integers (start, end).
        feature: A SeqFeature.

    Returns:
        A boolean indicating whether the features overlap.
    """
    interval_start = interval[0]
    interval_end = interval[1]

    if feature.location.start == interval_start:
        # >>>>>>>>>
        # (....)
        return interval_end - interval_start > 0

    elif feature.location.start < interval_start:
        # >>>>>>>>>
        #    (..........)
        return feature.location.end > interval_start

    else:
        #      >>>>>>>>>
        # (........)
        return feature.location.start < interval_end


def get_region_codon_indeces_in_feature(feature, region):
    """Returns a list of codon indeces in the feature that are overlapped
    by the region.

    Handles strand polarity correctly.

    Args:
        feature: A SeqFeature object.
        region: Tuple pair (start, end), relative to the overall genome record.
    """
    num_codons = len(feature) / 3
    assert num_codons > 0

    # Determine the region to substitute over (handle polarity below).
    first_codon_index = max(
        (region[0] - feature.location.start) / 3,
        0
    )
    last_codon_index = min(
            (region[1] - 1 - feature.location.start) / 3,
            num_codons - 1
    )
    codon_indeces_to_change = range(first_codon_index, last_codon_index + 1)
    if not codon_indeces_to_change:
        return []

    # Account for polarity.
    if feature.strand == -1:
        codon_indeces_to_change = [
                num_codons - 1 - idx for idx in codon_indeces_to_change]

    return sorted(codon_indeces_to_change)


def calc_interval_list_to_features_overlapped(interval_list, feature_list):
    """Returns a list with one-to-one correspondence with the interval_list
    that indicates which of the features in feature_list the interval
    overlaps.

    Args:
        interval_list: List of tuples representing intervals (start, end).
        feature_list: List of features to look for overlaps against.

    Returns:
        List of lists where each index corresponds to the index of the
        interval passed in and each list is the list of features that interval
        overlaps.
    """
    # Make sure the list is sorted by start.
    # This allows us to generally do a little better than a full n*m
    # comparison below.
    all(interval_list[i][0] <= interval_list[i + 1][0]
            for i in xrange(len(interval_list) -1))

    # Initialize the list to return.
    idx_to_features_overlapped_list = [[] for i in xrange(len(interval_list))]

    # For each feature, check which intervals it overlaps with.
    for feature in feature_list:
        # Iterate through the sorted list of intervals.
        for idx, interval in enumerate(interval_list):
            if does_interval_overlap_feature(interval, feature):
                idx_to_features_overlapped_list[idx].append(feature)
            if interval[0] > feature.location.end:
                # All remaining intervals are after the feature.
                break

    return idx_to_features_overlapped_list


def validate_start_seq(genome_record, start_position, validation_start_seq):
    """Utility method for checking that the sub-sequence starting at the
    given start_position is what we think it is.

    Args:
        genome_record: SeqRecord containing the sequence.
        start_position: Integer position that is the first base in the
            validation sequence.
        validation_start_seq: The starting sub-sequence that we are validating.

    Raises:
        AssertionError if validation fails.
    """
    orig_seq = str(genome_record.seq)
    actual_start_seq = orig_seq[
            start_position:start_position + len(validation_start_seq)]
    assert validation_start_seq == actual_start_seq, (
            "Pos: %d, Expected: %s, Actual: %s" % (
                    start_position, validation_start_seq, actual_start_seq))


def validate_end_seq(genome_record, end_position, validation_end_seq):
    """Utility method for checking that the sub-sequence ending at the given
    given start_pos is what we think it is.

    Args:
        genome_record: SeqRecord containing the sequence.
        end_position: Integer position that is the pythoning upper bound on the
            validation sequence. That is, it is the first position past
            the validation sequence end.
        validation_end_seq: The end sub-sequence that we are validating.

    Raises:
        AssertionError if validation fails.
    """
    orig_seq = str(genome_record.seq)
    actual_end_seq = orig_seq[
            end_position - len(validation_end_seq):end_position]
    assert validation_end_seq == actual_end_seq, (
            "Actual: %s" % actual_end_seq)


def swap_region_seq(
        refactor_context, region_start, original_region_seq, new_region_seq,
        fuzzy=False):
    """Swaps the new sequence at the position defined by the original sequence.

    This might be useful, for example, when manually fixing GC content issues
    in a region.

    For now, we require that the new sequence be the same length as the
    original sequence. This can be changed if we come up with appropriate
    use-cases.

    Args:
        refactor_context: The RefactorContext. The contained SeqRecord
            gets mutated. Also okay to pass SeqRecord directly.
        region_start: Pythonic position of the starting region.
        original_region_seq: Original sequence in that region.
        new_region_seq: New sequence in that region.
        fuzzy: If True, slide the target window around a bit to find the
            target. Experts only!

    Returns:
        Handle to the SeqRecord in the RefactorContext.
    """
    if isinstance(refactor_context, SeqRecord):
        updated_genome_record = refactor_context
    else:
        updated_genome_record = refactor_context.get_genome_record()

    already_ref = False
    try:
        validate_start_seq(updated_genome_record, region_start, new_region_seq)
        already_ref = True
    except:
        pass
    if already_ref:
        return updated_genome_record

    # Initial validation.
    assert len(original_region_seq) == len(new_region_seq), (
            "The lengths of the regions must be the same.")

    if fuzzy:
        try:
            validate_start_seq(updated_genome_record, region_start,
                    original_region_seq)
            target_found = True
        except:
            target_found = False
            for start_pos in range(region_start - 5,
                    min(region_start + 5, len(updated_genome_record.seq))):
                try:
                    validate_start_seq(updated_genome_record, start_pos,
                            original_region_seq)
                    region_start = start_pos
                    target_found = True
                except:
                    continue
        assert target_found
    else:
        validate_start_seq(updated_genome_record, region_start, original_region_seq)

    # Perform the swap.
    swap_region_size = len(original_region_seq)

    original_full_seq = updated_genome_record.seq
    new_full_seq = (
            original_full_seq[:region_start] +
            new_region_seq +
            original_full_seq[region_start + swap_region_size:])
    assert len(original_full_seq) == len(new_full_seq)

    updated_genome_record.seq = new_full_seq
    return updated_genome_record


def swap_unique_seq(genome_record, original_seq, replacement_seq,
        fuzzy_region=None, fuzzy_epsilon=100):
    """Finds the exact seq in the genome_record and swaps it out
    with the new one.

    For now, we require that replacement_seq be the same length as
    original_seq.

    Args:
        fuzzy_region: Optional tuple describing interval where we limit
            our attention to. Useful when original_seq is not unique
            genome-wide.
        fuzzy_epsilon: Amount of padding on either side of fuzzy_region.

    Raises:
        AssertionError if sequence not found or if more than one.
    """
    assert len(original_seq) == len(replacement_seq)
    seq_len = len(original_seq)

    # Find start positions (followed by fuzzy filter).
    start_positions = [m.start() for m in
            re.finditer(original_seq, str(genome_record.seq))]

    # Filter out matches that aren't in fuzzy region.
    if fuzzy_region:
        fuzzy_start = fuzzy_region[0] - fuzzy_epsilon
        fuzzy_end = fuzzy_region[1] + fuzzy_epsilon
        start_positions = [s for s in start_positions
                if s >= fuzzy_start and s < fuzzy_end]

    if not len(start_positions):
        raise AssertionError(
                "Sequence %s not found." % original_seq[:10] + '...')
    if len(start_positions) > 1:
        raise AssertionError(
                "Sequence %s not unique." % original_seq[:10] + '...')

    seq_start = start_positions[0]
    genome_record.seq = (
            genome_record.seq[:seq_start] +
            replacement_seq +
            genome_record.seq[seq_start + seq_len:])


def get_codon_feature_location(feature, codon_index):
    """Returns the feature location that could be used for annotating
    a codon inside of the given feature, at the given index.
    """
    if feature.strand == 1:
        start = feature.location.start + codon_index * 3
        end = start + 3
        return FeatureLocation(start, end)
    elif feature.strand == -1:
        end = feature.location.end - codon_index * 3
        start = end - 3
        return FeatureLocation(start, end)
    else:
        raise AssertionError("Unknown feature strand.")


def remove_gene_features(genome_record, gene_dict):
    """Removes features by gene and type in the given gene_list.

    Args:
        genome_record: SeqRecord that will have its list of features mutated.
        gene_dict: Dictionary with keys being gene names and values being
            dictionaries with at least the keys:
                'types': Feature types for that gene to remove.
    """
    updated_feature_list = []
    for feature in genome_record.features:
        if ('gene' in feature.qualifiers and
                feature.qualifiers['gene'][0] in gene_dict and
                feature.type in gene_dict[feature.qualifiers['gene'][0]]['types']):
            continue # skip appending it below, thus removing it.
        # Otherwise keep the feature.
        updated_feature_list.append(feature)
    genome_record.features = updated_feature_list


def make_gene_misc_feature(genome_record, gene_list, note=None):
    """Converts the passed in list of genes into misc_feature types.

    We use the simplest implementation for now which involves deleting
    the 'gene' feature and converting the feature with type 'CDS' to
    have type 'misc_feature'.
    """
    updated_features = []
    for feature in genome_record.features:

        if ('gene' in feature.qualifiers and
                feature.qualifiers['gene'][0] in gene_list):
            if feature.type == 'CDS':
                feature.type = 'misc_feature'
                feature.qualifiers = {}
                feature.qualifiers['note'] = note if note is not None else (
                    "Annotated as CDS in parent strain.")
            elif feature.type == 'gene':
                # Skip directly to next iteration, causing the feature to
                # be dropped.
                continue # Getting rid of this annotation.
            else:
                raise AssertionError("Unexpected type %s for gene %s" % (
                        feature.type, feature.qualifiers['gene'][0]))
        updated_features.append(feature)
    genome_record.features = updated_features


def update_feature_note(feature, update_text):
    """Updates the note feature qualifier value.

    If no note exists, create one. If it does exist, append the text
    separated from existing test by a semicolon.
    """
    if not 'note' in feature.qualifiers:
        feature.qualifiers['note'] = [update_text]
    else:
        current_note = feature.qualifiers['note'][0]
        feature.qualifiers['note'] = [current_note + '; ' + update_text]


def does_feature_end_with_valid_stop_codon(feature, seq_record,
        stop_codons=set(['TAA', 'TGA', 'TAG'])):
    """Checks whether a feature ends with an in-frame stop codon.

    Args:
        feature: A SeqFeature.
        seq_record: SeqRecord that feature comes from. Used to get the
            feature sequence.

    Returns:
        Boolean indicating whether the feature ends with a stop codon.
    """
    if len(feature) % 3 != 0:
        return False
    seq = str(feature.extract(seq_record.seq))
    last_codon = seq[-3:]
    return last_codon in stop_codons


def does_feature_start_with_valid_start_codon(feature, seq_record,
        start_codons=set(['ATG', 'GTG', 'TTG', 'ATT', 'CTG'])):
    """Checks whether a feature ends with an in-frame stop codon.

    Args:
        feature: A SeqFeature.
        seq_record: SeqRecord that feature comes from. Used to get the
            feature sequence.

    Returns:
        Boolean indicating whether the feature ends with a stop codon.
    """
    seq = str(feature.extract(seq_record.seq))
    return seq[0:3] in start_codons


def maybe_get_feature_gene(feature):
    """Utitliy functoin that returns the string name for the gene in the
    feature if it has one. Otherwise returns None.

    # DEPRECATED: Use get_feature_gene() instead.
    """
    return get_feature_gene(feature)


def get_features_with_gene(seq_record, gene):
    """Returns the features in seq_record that have gene in their qualifiers.
    """
    matches = []
    for feature in seq_record.features:
        maybe_gene = get_feature_gene(feature)
        if maybe_gene and maybe_gene == gene:
            matches.append(feature)
    return matches


def build_gene_to_CDS_feature_map(genome_record):
    """Builds a map from gene name to CDS feature for the given genome record.
    """
    gene_to_feature_map = {}
    cds_features = [feature for feature in genome_record.features
            if feature.type == 'CDS']
    for feature in cds_features:
        maybe_gene = get_feature_gene(feature)
        if not maybe_gene:
            continue
        gene_to_feature_map[maybe_gene] = feature
    return gene_to_feature_map


def build_gene_synonym_map(genome_record):
    gene_to_synonym_map = {}
    cds_features = [feature for feature in genome_record.features
            if feature.type in set(['CDS', 'gene'])]
    for feature in cds_features:
        maybe_gene = get_feature_gene(feature)
        if not maybe_gene:
            continue

        synonym_list = [maybe_gene]

        # Add the locus tag as a synonym.
        if 'locus_tag' in feature.qualifiers:
            synonym_list.append(feature.qualifiers['locus_tag'][0])

        if 'gene_synonym' in feature.qualifiers:
            # Build a list containing all synonyms which will serve as the
            # value of the bimap.

            # Check each to see if it can be split.
            for synonym_phrase in feature.qualifiers['gene_synonym']:
                split_phrase = synonym_phrase.split(';')
                for synonym in split_phrase:
                    clean_synonym = synonym.strip()
                    if len(clean_synonym) > 0:
                        synonym_list.append(clean_synonym)

        # Remove duplicates.
        synonym_list = list(set(synonym_list))

        # Add a key for each synonym to this list.
        for synonym in synonym_list:
            gene_to_synonym_map[synonym] = synonym_list

    return gene_to_synonym_map


def build_linked_cds_feature_map(genome_record, gene_synonym_map=None):
    """Returns a map linkig cds feature gene names in spatial
    order to the feature following that gene.
    """
    linked_cds_map = {}

    if not gene_synonym_map:
        gene_synonym_map = build_gene_synonym_map(genome_record)

    # Get CDS features sorted by start position.
    cds_features = [feature for feature in genome_record.features
            if feature.type == 'CDS']
    cds_features = sorted(cds_features,
            key=lambda feature: feature.location.start)

    for idx in range(len(cds_features) - 1):
        this_feature = cds_features[idx]
        next_feature = cds_features[idx + 1]
        this_feature_gene = get_feature_gene(this_feature)
        next_feature_gene = get_feature_gene(next_feature)

        # Cover synonyms.
        for synonym in gene_synonym_map[this_feature_gene]:
            linked_cds_map[synonym] = next_feature

    return linked_cds_map


def fuzzy_search(genome_record_or_seq, target_seq, search_window_center,
        search_window_size=200, compare_step=1):
    """Performs search across the search window for the given sequence.

    Returns:
        An integer indicating the start position of the sequence if it is found
        in the window. Returns None if not found.
    """
    assert compare_step > 0

    if isinstance(genome_record_or_seq, SeqRecord):
        genome_seq = str(genome_record_or_seq.seq).upper()
    else:
        genome_seq = genome_record_or_seq

    target_seq = target_seq.upper()

    window_range = (search_window_center - search_window_size / 2,
            search_window_center + search_window_size / 2)

    # We wiggle the window_range outward from the center.
    window_start = search_window_center
    next_delta = 1

    while window_range[0] <= window_start and window_start <= window_range[1]:
        match = True
        for target_idx in range(0, len(target_seq), compare_step):
            if target_seq[target_idx] == genome_seq[window_start + target_idx]:
                continue
            else:
                match = False
                break
        if match:
            return window_start

        # Perform wiggle and prepare next wiggle.
        window_start = search_window_center + next_delta
        if next_delta < 0:
            next_delta *= -1
            next_delta += 1
        else:
            next_delta *= -1

    return None


def get_feature_label(feature, default=None):
    """Returns the feature label, or None if it doesn't exist.
    """
    try:
        return feature.qualifiers['label'][0]
    except:
        return default


def set_feature_label(feature, label):
    """Sets the label qualifier for a SeqFeature.
    """
    feature.qualifiers['label'] = [label]


def get_feature_by_gene_name(gene_name, gene_type, genome_record):
    """Returns feature with gene name and type, or None if not found.

    Raises:
        AssertionError if more than one match found.
    """
    assert gene_name

    matches = []
    for f in genome_record.features:
        if not f.type == gene_type:
            continue
        if gene_name == get_feature_gene(f):
            matches.append(f)

    assert len(matches) <= 1, (
            "More than one feature with type %s and gene_name %s" % (
                    gene_type, gene_name))

    if not len(matches):
        return None

    return matches[0]
