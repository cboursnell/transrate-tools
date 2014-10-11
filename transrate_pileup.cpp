#include "transrate_pileup.h"

using namespace BamTools;

using namespace std;

TransratePileup::TransratePileup() {
  ref_length = 0;
  // ref_length = rl;
  // coverage.reserve(ref_length);
}

void TransratePileup::clearCoverage(int rl) {
  ref_length = rl;
  coverage.clear();
  coverage.reserve(ref_length);
  for (int i = 0; i < ref_length; ++i) {
    coverage[i] = 0;
  }

  mapq.clear();
  mapq.reserve(ref_length);
  for (int i = 0; i < ref_length; ++i) {
    mapq[i] =0;
  }

}

int TransratePileup::getCoverage(int i) {
  return coverage[i];
}

int TransratePileup::getMapq(int i) {
  return mapq[i];
}

vector<int> TransratePileup::getCoverageArray() {
  return coverage;
}

vector<int> TransratePileup::getBinnedCoverage() {
  return coverage;
}

int TransratePileup::getBasesUncovered() {
  int bases_uncovered=0;
  for (int i = 0; i < ref_length; i++) {
    if (coverage[i]==0) {
      bases_uncovered++;
    }
  }
  return bases_uncovered;
}

double TransratePileup::getUniqueBases() {
  double total = 0;
  for (int i = 0; i < ref_length; i++) {
    if (coverage[i] > 0) {
      // cout << sqrt (mapq[i] / (double)coverage[i]) << " " << endl;
      total += log (sqrt (mapq[i] / (double)coverage[i]));
    }
  }
  return exp (total / (double)ref_length);
  // return pow (exp(total), 1 / (double)ref_length);
}

double TransratePileup::p_not_segmented() {
  vector<int> states(30);
  int bin_width = (ref_length / 30)+1; // ceil
  int counter=0;
  int pos=0;
  int total=0;
  for (int i = 0; i < ref_length; ++i) {
    total += coverage[i];
    counter++;
    if (counter==bin_width || i==ref_length-1) {
      states[pos] = max(0.0,log2(total/counter));
      pos++;
      counter=0;
      total=0;
    }
  }
  //
  Segmenter segmenter(states);
  return segmenter.prob_k_given_R(0); // prob of 0 change points (1 segment)
}

bool TransratePileup::addAlignment(const BamAlignment& alignment) {
  if ( !alignment.IsMapped()) {
    return false;
  }
  int numCigarOps = alignment.CigarData.size();
  int pos = alignment.Position;
  for (int i = 0; i < numCigarOps; ++i) {
    const CigarOp& op = alignment.CigarData.at(i);
    if (op.Type == 'M') {
      for(p = 0; p < (int)op.Length; ++p) {
        pos++;
        coverage[pos-1]++;
        mapq[pos-1] += pow(alignment.MapQuality,2);
      }
    }
    if (op.Type == 'I') {
      // I means insertion - gaps in reference
    }
    if (op.Type == 'D') {
      // D means deletion - gaps in read
      for(p = 0; p < (int)op.Length; ++p) {
        pos++;
      }
    }
  }
  return true;
}

