#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "api/BamReader.h"
#include "pileup.h"

using namespace BamTools;

using namespace std;

double nullprior = 0.7;

class BamRead
{
public:
  int estimate_fragment_size(std::string file);
  int load_bam(std::string);
  int bar;
  int seq_count;
  uint32_t nm_tag;
  int ldist;
  int realistic_distance;
  int refid;
  BamReader reader;
  BamAlignment alignment;
  std::vector<TransratePileup> array;
};
