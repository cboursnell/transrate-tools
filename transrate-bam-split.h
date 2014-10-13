#include <math.h> /* modf, tgamma */
#include <float.h> /* FLT_MIN */
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace BamTools;

using namespace std;

int max_indel_size = 10;

// counters
int unmapped = 0;
int not_proper = 0;
int mate_not_mapped = 0;
int mate_diff_contig = 0;
int both_same_orient = 0;
int inversed_left = 0;
int inversed_right = 0;
int oversize_ins = 0;
int oversize_del = 0;

/**
 * The Indel struct stores the information for a single insertion or deletion.
 *  @author    Adam Roberts
 *  @date      2012
 *  @copyright Artistic License 2.0
 **/
struct Indel {
  /**
   * A public size_t for the position of the Indel in the read. 0-based.
   */
  size_t pos;
  /**
   * A public size_t for the length of the Indel in the read.
   */
  size_t len;
  /**
   * Indel constructor
   */
  Indel(size_t p, size_t l) : pos(p), len(l) {}
};

class BamSplitter {

  // ctor & dtor
  public:
    BamSplitter(string& input)
      : m_input(input)
    { }

    ~BamSplitter(void) {
      m_reader.Close();
    }

  // 'public' interface
  public:
    bool Run(void);

  private:
     // open BamWriter
    bool OpenWriters(void);
    // close & delete BamWriters
    void CloseWriters(void);
    // open our BamReader
    bool OpenReader(void);
    // split alignments in BAM file based on eXpress's filtering
    bool Split(void);
    // check if an alignment is valid
    bool IsAlignmentValid(BamAlignment& a);

    // data members
    string m_input;
    string m_output;
    BamReader m_reader;
    BamWriter m_validwriter;
    BamWriter m_invalidwriter;
    string m_header;
    RefVector m_references;

};
