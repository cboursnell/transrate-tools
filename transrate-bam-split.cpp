#include "transrate-bam-split.h"

// Transrate BAM splitter
//
// Takes a path to a BAM file as input and split the file
// into valid and invalid alignments, according to the criteria
// used by eXpress.
//
// By Richard Smith-Unna 2014
//
// MIT license
//
// Uses code from BamTools and eXpress

string RemoveFilenameExtension(const string& filename) {
    size_t found = filename.rfind(".");
    return filename.substr(0, found);
}

bool BamSplitter::Run(void) {

  m_output = RemoveFilenameExtension(m_input);

  // open up BamReader
  if (!OpenReader()) {
    return false;
  }

  // open up BamWriters
  if (!OpenWriters()) {
    return false;
  }

  // do the splitting
  Split();

  return true;
}

bool BamSplitter::OpenReader(void) {

  // attempt to open BAM file
  if (!m_reader.Open(m_input)) {
    cerr << "bam-split ERROR: could not open BAM file: "
         << m_input << endl;
    return false;
  }

  // save file 'metadata' & return success
  m_header     = m_reader.GetHeaderText();
  m_references = m_reader.GetReferenceData();
  return true;
}

bool BamSplitter::OpenWriters(void) {

  const string validFilename = m_output + ".valid.bam";
  if (!m_validwriter.Open(validFilename, m_header, m_references)) {
    cerr << "bam-split ERROR: could not open "
         << validFilename << " for writing." << endl;
    return false;
  }

  const string invalidFilename = m_output + ".invalid.bam";
  if (!m_invalidwriter.Open(invalidFilename, m_header, m_references)) {
    cerr << "bam-split ERROR: could not open "
         << invalidFilename << " for writing." << endl;
    return false;
  }

  return true;
}

void BamSplitter::CloseWriters(void) {

  m_validwriter.Close();
  m_invalidwriter.Close();

}

bool BamSplitter::Split(void) {

  // iterate through alignments
  BamAlignment a;
  BamWriter *writer;
  bool valid;

  while (m_reader.GetNextAlignment(a)) {

    if (IsAlignmentValid(a)) {
      m_validwriter.SaveAlignment(a);
    } else {
      m_invalidwriter.SaveAlignment(a);
    }

  }

  // clean up BamWriters
  CloseWriters();

  return true;
}

bool BamSplitter::IsAlignmentValid(BamAlignment& a) {

  if (!a.IsMapped()) {
    return 0;
  }

  bool is_paired = a.IsPaired();

  if (is_paired && !a.IsProperPair()) {
    return 0;
  }

  bool is_reversed = a.IsReverseStrand();
  bool is_mate_reversed = a.IsMateReverseStrand();
  if (is_paired && (!a.IsMateMapped()
                || a.RefID != a.MateRefID
                || is_reversed == is_mate_reversed
                || (is_reversed && a.MatePosition > a.Position)
                || (is_mate_reversed && a.MatePosition < a.Position))) {
    return 0;
  }

  return 1;
}

int main (int argc, char* argv[]) {

  if (argc == 2) {
    string infile(argv[1]);

    BamSplitter splitter(infile);
    
    if (splitter.Run()) {
      return 0;
    } else {
      return 1;
    }

  } else {
    cout << "bam-split version 0.1\n"
            "Usage:\n"
            "bam-split <bam_file>" << endl;
    return 1;
  }

}
