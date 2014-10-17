
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "api/BamReader.h"
#include "transrate_pileup.h"

using namespace BamTools;

using namespace std;

struct ContigRecord {
  int bases_mapped;
  int p_seq_true;
  int bridges;
  int length;
  string name;
  int fragments_mapped;
  int both_mapped;
  int properpair;
  int good;
  int bases_uncovered;
  double p_unique;
  double p_not_segmented;
};

class BetterBam {
    int realistic_distance;
    int seq_count;
    int i,j;
    uint32_t nm_tag;
    int ldist;
    int rdist;
    int maxL=0;
    std::string file;
    BamReader reader;
    BamAlignment alignment;
  public:
    std::vector<ContigRecord> array;
    BetterBam (std::string);

    void set_fragment_size(int size, int sd) {
      realistic_distance = size + 3 * sd;
    }

    int load_bam() {
      if (!reader.Open(file)) {
        cerr << "Could not open BAM file" << endl;
        return 1;
      }
      // get sam header
      SamSequenceDictionary dictionary = reader.GetHeader().Sequences;
      seq_count = dictionary.Size();
      array.resize(seq_count);
      // get an iterator for looping over the sequences in the header
      std::vector<SamSequence>::iterator it = dictionary.Begin();
      // fill the vector with intial values
      for (i = 0; i < seq_count; i++) {
        array[i].bases_mapped = 0;
        array[i].p_seq_true = 0;
        array[i].bridges = 0;
        if (it[i].HasLength()) {
          array[i].length = atoi(it[i].Length.c_str());
          if (array[i].length > maxL) {
            maxL = array[i].length;
          }
        } else {
          array[i].length = 0;
        }
        array[i].name = it[i].Name;
        array[i].fragments_mapped = 0;
        array[i].both_mapped = 0;
        array[i].properpair = 0;
        array[i].good = 0;
        array[i].bases_uncovered = 0;
        array[i].p_unique = 0;
        array[i].p_not_segmented = 1;
      }
      // loop through bam file
      i = -2;
      TransratePileup pileup(maxL);
      int ref_length = -1;
      while (reader.GetNextAlignment(alignment)) {
        if (alignment.IsMapped()) {
          // new contig
          if (alignment.RefID != i) {
            if (i>=0) {
              array[i].bases_uncovered = pileup.getBasesUncovered();
              array[i].p_unique = pileup.getUniqueBases();
              array[i].p_not_segmented = pileup.p_not_segmented();
            }
            i = alignment.RefID;
            ref_length = array[i].length;
            pileup.clearCoverage(ref_length);

          }
          if (alignment.IsPrimaryAlignment()) {
            pileup.addAlignment(alignment);
          }

          array[i].bases_mapped += alignment.Length;
          if (alignment.GetTag("NM", nm_tag)) {
            array[i].p_seq_true += nm_tag;
          }
          if (alignment.IsFirstMate() ||
            (alignment.IsSecondMate() && !alignment.IsMateMapped())) {
            array[i].fragments_mapped++;
          }
          if (alignment.IsFirstMate() && alignment.IsMateMapped()) {
            array[i].both_mapped++;
            if (alignment.IsProperPair() && alignment.RefID==alignment.MateRefID) {
              array[i].properpair++;
              // check orientation

              ldist = max(alignment.Position-alignment.MatePosition,
                          alignment.MatePosition-alignment.Position);
              if (ldist < realistic_distance) {
                if (!alignment.IsReverseStrand() && alignment.IsMateReverseStrand()) {
                  if (alignment.GetEndPosition() < alignment.MatePosition) {
                    array[i].good++;
                  }
                } else if (alignment.IsReverseStrand() && !alignment.IsMateReverseStrand()) {
                  if (alignment.GetEndPosition() > alignment.MatePosition) {
                    array[i].good++;
                  }
                }
              }
            } else {
              if (i != alignment.MateRefID) {
                ldist = min(alignment.Position,
                            array[i].length - alignment.Position);
                rdist = min(alignment.MatePosition,
                            array[alignment.MateRefID].length - alignment.MatePosition);
                if (ldist + rdist <= realistic_distance) {
                  array[i].bridges++;
                }
              }
            }
          }
        }

      }
      array[i].bases_uncovered = pileup.getBasesUncovered();
      array[i].p_unique = pileup.getUniqueBases();
      array[i].p_not_segmented = pileup.p_not_segmented();

      reader.Close();
      return 0;
    }

    int get_seq_count() {
      return seq_count;
    }

    ContigRecord get_info(int i) {
      return array.at(i);
    }

    int get_bases_mapped(int i) {
      return array.at(i).bases_mapped;
    }
};

//constructor
BetterBam::BetterBam (std::string s) {
    file = s;
    realistic_distance = 350;
}

int main (int argc, char* argv[]) {
  if (argc == 3) {
    string infile = argv[1];
    BetterBam bam (infile);
    bam.load_bam();
    int i;

    std::ofstream output;
    output.open (argv[2]);
    output << "name,p_seq_true,bridges,length,fragments_mapped,"
              "both_mapped,properpair,good,bases_uncovered,p_unique,"
              "p_not_segmented\n";
    for (i = 0; i < bam.get_seq_count(); i++) {
      output << bam.get_info(i).name << ",";
      if (bam.get_info(i).bases_mapped>0) {
        output << 1-((double)bam.get_info(i).p_seq_true/bam.get_info(i).bases_mapped) << ",";
      } else {
        output << "1,";
      }
      output << bam.get_info(i).bridges << ",";
      output << bam.get_info(i).length << ",";
      output << bam.get_info(i).fragments_mapped << ",";
      output << bam.get_info(i).both_mapped << ",";
      output << bam.get_info(i).properpair << ",";
      output << bam.get_info(i).good << ",";
      output << bam.get_info(i).bases_uncovered << ",";
      output << bam.get_info(i).p_unique << ",";
      output << bam.get_info(i).p_not_segmented << endl;
    }
    output.close();
    return 0;
  } else {
    cout << "bam-read version 0.3.1\nUsage:\nbam-read <bam_file> <output_csv>" << endl;
    return 1;
  }
}
