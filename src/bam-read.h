#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "api/BamReader.h"
#include "pileup.h"
#include "thread.h"

using namespace BamTools;

using namespace std;

double nullprior = 0.7;

class BamRead
{
public:
  int estimate_fragment_size(std::string file);
  int load_bam(std::string, int);
  int thread_id;
  int bar;
  int seq_count;
  int realistic_distance;
  BamReader reader;
  BamAlignment alignment;
  std::vector<TransratePileup> array;   // data about each contig
  std::vector<TransrateThread> threads;
};
