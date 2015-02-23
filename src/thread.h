#include <list>
#include "api/BamAlignment.h"

using namespace BamTools;

class TransrateThread {

  private:
    list<BamAlignment> queue;

  public:
    TransrateThread();

    // instance variables
    double seq_true;
    double scale;
    BamAlignment alignment;

    //methods
    int add(const BamAlignment& a);
    int process();

};