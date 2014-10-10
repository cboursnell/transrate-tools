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
