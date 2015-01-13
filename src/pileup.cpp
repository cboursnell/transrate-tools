#include "pileup.h"

using namespace BamTools;

using namespace std;

//constructor
TransratePileup::TransratePileup() {
  ref_length = 100;
  bases_mapped = 0;
  p_seq_true = 0.0;
  bridges = 0;
  name = "unknown";
  reads_mapped = 0;
  fragments_mapped = 0;
  both_mapped = 0;
  properpair = 0;
  good = 0;
  bases_uncovered = 0;
  p_unique = 0;
  p_not_segmented = 0;
  coverage.resize(ref_length);
}

void TransratePileup::setName(std::string _name) {
  name = _name;
}

void TransratePileup::setLength(int len) {
  ref_length = len;
  coverage.resize(len);
}

int TransratePileup::getCoverage(int i) {
  return coverage[i];
}

vector<int> TransratePileup::getCoverageArray() {
  return coverage;
}

void TransratePileup::calculateUncoveredBases() {
  for (int i = 0; i < ref_length; ++i) {
    if (coverage[i]==0) {
      ++bases_uncovered;
    }
  }
}

void TransratePileup::setPNotSegmented() {
  vector<int> states(30);
  int bin_width = (ref_length / 30)+1; // ceil
  int counter=0;
  int pos=0;
  int total=0;
  for (int i = 0; i < ref_length; ++i) {
    total += coverage[i];
    counter++;
    if (counter==bin_width || i==ref_length-1) {
      states[pos] = (int)min(24.0, max(0.0,log2(total/counter)));
      pos++;
      counter=0;
      total=0;
    }
  }
  //
  Segmenter segmenter(states, nullprior);
  p_not_segmented = segmenter.prob_k_given_R(0); // prob of 0 change points (1 segment)
}

void TransratePileup::addAlignment(const BamAlignment& alignment) {
  int numCigarOps = alignment.CigarData.size();
  int pos = alignment.Position;
  bases_mapped += alignment.Length;
  for (int i = 0; i < numCigarOps; ++i) {
    const CigarOp& op = alignment.CigarData.at(i);
    if (op.Type == 'M') {
      // increase coverage by 1 for each position that is a match
      for(p = 0; p < (int)op.Length; ++p) {
        if (pos < ref_length) {
          ++coverage[pos];
        }
        ++pos;
      }
    }
    if (op.Type == 'I') {
      // I means insertion - gaps in reference
    }
    if (op.Type == 'D' || op.Type == 'S') {
      // D means deletion - gaps in read
      pos += (int)op.Length;
    }
  }
}
