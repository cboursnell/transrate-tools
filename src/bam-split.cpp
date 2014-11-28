#include "bam-split.h"

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

/**
 * A helper function that strips the extension from a filename.
 */
string strip_extension(const string& filename) {
    size_t found = filename.rfind(".");
    return filename.substr(0, found);
}

/**
 * A helper functon that calculates the length of the reference spanned by the
 * read and populates the indel vectors (for BAM input).
 * @param cigar_vec a vector containing the split cigar string.
 * @param inserts an empty Indel vector into which to add inserts.
 * @param deletes an empty Indel vector into which to add deletions.
 */
size_t cigar_length(vector<BamTools::CigarOp>& cigar_vec,
                    vector<Indel>& inserts, vector<Indel>& deletes) {
  inserts.clear();
  deletes.clear();
  size_t i = 0; // read index
  size_t j = 0; // genomic index
  for (size_t k = 0; k < cigar_vec.size(); ++k) {
    char op_char = cigar_vec[k].Type;
    size_t op_len = cigar_vec[k].Length;
    switch(op_char) {
      case 'I':
      case 'S':
        inserts.push_back(Indel(i, op_len));
        i += op_len;
        break;
      case 'D':
        deletes.push_back(Indel(i, op_len));
        j += op_len;
        break;
      case 'M':
      case 'N':
        i += op_len;
        j += op_len;
        break;
    }
  }
  return j;
}

bool BamSplitter::Run(void) {

  m_output = strip_extension(m_input);

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
    ++unmapped;
    return 0;
  }

  bool is_paired = a.IsPaired();

  if (is_paired && !a.IsProperPair()) {
    ++not_proper;
    return 0;
  }

  bool is_reversed = a.IsReverseStrand();
  bool is_mate_reversed = a.IsMateReverseStrand();
  if (is_paired) {
    if (!a.IsMateMapped()) {
      ++mate_not_mapped;
      return 0;
    } else if (a.RefID != a.MateRefID) {
      ++mate_diff_contig;
      return 0;
    } else if (is_reversed == is_mate_reversed) {
      ++both_same_orient;
      return 0;
    } else if (is_reversed && a.MatePosition > a.Position) {
      ++inversed_left;
      return 0;
    } else if (is_mate_reversed && a.MatePosition < a.Position) {
      ++inversed_right;
      return 0;
    }
  }

  // loop through all cigar operations and check they are not S
  int numCigarOps = a.CigarData.size();
  bool check = false;
  for (int i = 0; i < numCigarOps; ++i) {
    const CigarOp& op = a.CigarData.at(i);
    if (op.Type != 'S') {
      check = true;
    }
  }
  if (!check) {
    return 0;
  }

  std::vector<Indel> inserts;
  std::vector<Indel> deletes;

  cigar_length(a.CigarData, inserts, deletes);

  for (Indel& indel : inserts) {
    if (indel.len > max_indel_size) {
      ++oversize_ins;
      return 0;
    }
  }

  for (Indel& indel : deletes) {
    if (indel.len > max_indel_size) {
      ++oversize_del;
      return 0;
    }
  }

  return 1;
}

int main (int argc, char* argv[]) {

  if (argc == 2) {
    string infile(argv[1]);

    BamSplitter splitter(infile);

    if (splitter.Run()) {
      cout << "Unmapped:\t\t\t" << to_string(unmapped)
           << "\nNot proper:\t\t" << to_string(not_proper)
           << "\nMate not mapped:\t" << to_string(mate_not_mapped)
           << "\nMate diff contig:\t" << to_string(mate_diff_contig)
           << "\nBoth same orient:\t" << to_string(both_same_orient)
           << "\nInversed left:\t\t" << to_string(inversed_left)
           << "\nInversed right:\t\t" << to_string(inversed_right)
           << "\nOversize insert\t\t" << to_string(oversize_ins)
           << "\nOversize delete\t\t" << to_string(oversize_del)
           << endl;
      return 0;
    } else {
      return 1;
    }

  } else {
    cout << "bam-split version 1.0.0.beta1\n"
            "Usage:\n"
            "bam-split <bam_file>" << endl;
    return 1;
  }

}
