#include "thread.h"

//constructor
TransrateThread::TransrateThread() {
  // needs to store a pointer to the array object from `bam-read`
  scale = 0.65;
}

int TransrateThread::add(const BamAlignment& alignment) {
  queue.push_back(alignment);
  return 1;
}

int TransrateThread::process() {
  // if queue is not empty
  if (!list.empty()) {
  // remove item from front of the queue
    alignment = queue.pop_front();
    refid = alignment.RefID;

    ++array[refid].reads_mapped;

    array[refid].addAlignment(alignment);

    // store edit distance for sequence accuracy calculation
    if (alignment.HasTag("NM")) {
      if (alignment.GetTag("NM", nm_tag)) {
        len = alignment.Length;
        scale = (len - 35)/(double)len;
        seq_true = (len - nm_tag)/(double)len;
        seq_true = (seq_true - scale) * (1 / (1 - scale));
        array[refid].p_seq_true += seq_true;
      }
    }

    // count fragments where either or both mates mapped
    if (alignment.IsFirstMate() ||
      (alignment.IsSecondMate() && !alignment.IsMateMapped())) {
      ++array[refid].fragments_mapped;
    }

    // from now on ignore fragments unless both mates mapped
    if (!(alignment.IsFirstMate() && alignment.IsMateMapped())) {
      continue;
    }

    ++array[refid].both_mapped;

    // // count proper pairs, although we have our own definition because
    // // not all aligners use the same definition
    if (alignment.IsProperPair()) {
      ++array[refid].properpair;
    }

    // // mates must align to same contig, otherwise we record a bridge
    if (refid != alignment.MateRefID) {
      ++array[refid].bridges;
      continue;
    }

    ldist = max(alignment.Position-alignment.MatePosition,
                alignment.MatePosition-alignment.Position);
    if (ldist > realistic_distance) {
      // mates are too far apart
      continue;
    }

    // read orientation must match the generated library
    // in this case we only test for FR/RF orientation,
    // that is - we expect mates to be on opposite strands
    bool is_reversed = alignment.IsReverseStrand();
    bool is_mate_reversed = alignment.IsMateReverseStrand();

    if (!is_reversed && is_mate_reversed) {
      // in FR orientation, first read must start
      // before second read
      if (alignment.Position < alignment.MatePosition) {
        array[refid].good++;
      }
    } else if (is_reversed && !is_mate_reversed) {
      // in RF orientation, second read must start
      // before first read
      if (alignment.MatePosition < alignment.Position) {
        array[refid].good++;
      }
    }
  }

  return 1;

}