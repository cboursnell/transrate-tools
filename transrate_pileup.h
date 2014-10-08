#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include "api/BamAlignment.h"

using namespace BamTools;

using namespace std;

class TransratePileup {

  private:
    int p;
    int ref_length;
    vector<int> coverage;
    vector<long> mapq;

  public:
    TransratePileup();

    void clearCoverage(int rl);
    int getCoverage(int i);
    int getMapq(int i);
    vector<int> getCoverageArray();
    vector<int> getBinnedCoverage();
    int getBasesUncovered();
    double getUniqueBases();
    bool addAlignment(const BamAlignment& alignment);

};