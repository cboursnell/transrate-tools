#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include "api/BamAlignment.h"
#include "segmenter.h"

extern double nullprior;

using namespace BamTools;

using namespace std;

class TransratePileup {

  private:
    int p;
    int ref_length;
    vector<int> coverage;
    vector<long> mapq;

  public:
    TransratePileup(int maxL);

    void clearCoverage(int rl);
    int getCoverage(int i);
    int getMapq(int i);
    vector<int> getCoverageArray();
    vector<int> getBinnedCoverage();
    int getBasesUncovered();
    double getUniqueBases();
    double p_not_segmented();
    bool addAlignment(const BamAlignment& alignment);

};
