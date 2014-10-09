
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
  int edit_distance;
  int bridges;
  int length;
  string name;
  int reads_mapped;
  int both_mapped;
  int properpair;
  int good;
  int bases_uncovered;
  double mapq;
  double p_not_segmented;
};

class BetterBam {
    int realistic_distance;
    int seq_count;
    int i,j;
    uint32_t nm_tag;
    int ldist;
    int rdist;
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
      // file the vector with intial values
      for (i = 0; i < seq_count; i++) {
        array[i].bases_mapped = 0;
        array[i].edit_distance = 0;
        array[i].bridges = 0;
        if (it[i].HasLength()) {
          array[i].length = atoi(it[i].Length.c_str());
        } else {
          array[i].length = 0;
        }
        array[i].name = it[i].Name;
        array[i].reads_mapped = 0;
        array[i].both_mapped = 0;
        array[i].properpair = 0;
        array[i].good = 0;
        array[i].bases_uncovered = 0;
        array[i].mapq = 0;
        array[i].p_not_segmented = 1;
      }
      // loop through bam file
      i = -2;
      TransratePileup pileup;
      int ref_length = -1;
      while (reader.GetNextAlignment(alignment)) {

        if (alignment.IsMapped()) {
          // new contig
          if (alignment.RefID != i) {
            // -1 for unaligned reads
            if (i>=0) {
              array[i].bases_uncovered = pileup.getBasesUncovered();
              array[i].mapq = pileup.getUniqueBases();
              array[i].p_not_segmented = pileup.p_not_segmented();
            }
            i = alignment.RefID;
            ref_length = array[i].length;
            pileup.clearCoverage(ref_length);

          }
          if (alignment.IsPrimaryAlignment()) {
            pileup.addAlignment(alignment);
          }
        }

        array[i].bases_mapped += alignment.Length;
        if (alignment.GetTag("NM", nm_tag)) {
          array[i].edit_distance += nm_tag;
        }
        if (alignment.IsFirstMate() ||
          (alignment.IsSecondMate() && !alignment.IsMateMapped())) {
          array[i].reads_mapped++;
        }
        if (alignment.IsFirstMate() && alignment.IsMateMapped()) {
          array[i].both_mapped++;
          if (alignment.IsProperPair()) {
            array[i].properpair++;
            // check orientation
            if ((!alignment.IsReverseStrand() && alignment.IsMateReverseStrand()) ||
              (alignment.IsReverseStrand() && !alignment.IsMateReverseStrand())) {
              array[i].good++;
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
      array[i].bases_uncovered = pileup.getBasesUncovered();
      array[i].mapq = pileup.getUniqueBases();
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
    // string outfile = argv[2];
    output.open (argv[2]);
    output << "name,bases,edit_distance,bridges,length,reads_mapped,"
              "both_mapped,properpair,good,bases_uncovered,mapq"
              ",p_not_segmented\n";
    for (i = 0; i < bam.get_seq_count(); i++) {
      output << bam.get_info(i).name << ",";
      output << bam.get_info(i).bases_mapped << ",";
      output << bam.get_info(i).edit_distance << ",";
      output << bam.get_info(i).bridges << ",";
      output << bam.get_info(i).length << ",";
      output << bam.get_info(i).reads_mapped << ",";
      output << bam.get_info(i).both_mapped << ",";
      output << bam.get_info(i).properpair << ",";
      output << bam.get_info(i).good << ",";
      output << bam.get_info(i).bases_uncovered << ",";
      output << bam.get_info(i).mapq << ",";
      output << bam.get_info(i).p_not_segmented << endl;
    }
    output.close();
    return 0;
  } else {
    cout << "bam-read version 0.3\nUsage:\nbam-read <bam_file> <output_csv>" << endl;
    return 1;
  }
}