#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include "api/BamAlignment.h"
#include "segmenter.h"

extern double nullprior;
extern int fragment_length;

using namespace BamTools;

using namespace std;

class TransratePileup {

  private:
    // instance variables
    int p;
    vector<int> coverage;


  public:
    // constructor
    TransratePileup();

    // instance variables i can't be bothered to write accessor methods for
    int ref_length;
    int bases_mapped;
    double p_seq_true;
    int bridges;
    int length;
    string name;
    int reads_mapped;
    int fragments_mapped;
    int both_mapped;
    int properpair;
    int good;
    int bases_uncovered;
    int right_overhang;
    int left_overhang;
    double p_unique;
    double p_not_segmented;

    // methods
    void clearCoverage(int rl);
    void setName(std::string);
    void setLength(int);
    int getCoverage(int i);
    vector<int> getCoverageArray();
    void calculateUncoveredBases();
    void setPNotSegmented();
    void addAlignment(const BamAlignment& alignment);

};
