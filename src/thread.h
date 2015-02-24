#include <list>
#include "pileup.h"
#include "api/BamAlignment.h"

using namespace BamTools;

extern int realistic_distance;

class TransrateThread {

  public:
    TransrateThread();

    // instance variables
    double seq_true;
    double scale;
    BamAlignment alignment;
    std::list<BamAlignment> queue;
    int refid;
    int len;
    int ldist;
    int realistic_distance;
    uint32_t nm_tag;

    std::vector<TransratePileup> array;

    //methods
    int add(const BamAlignment& a);
    int process();

};