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

    // set realistic distance between pairs to have 0.03% false positive rate
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

        // read must be mapped
        if (!alignment.IsMapped()) {
          continue;
        }

        // check this read comes from the currently loaded contig
        // if not, load the new contig
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

        pileup.addAlignment(alignment);

        array[i].bases_mapped += alignment.Length;

        // store edit distance for sequence accuracy calculation
        if (alignment.HasTag("NM")) {
          if (alignment.GetTag("NM", nm_tag)) {
            array[i].p_seq_true += nm_tag;
          }
        }

        // count fragments where either or both mates mapped
        if (alignment.IsFirstMate() ||
          (alignment.IsSecondMate() && !alignment.IsMateMapped())) {
          array[i].fragments_mapped++;
        }

        // from now on ignore fragments unless both mates mapped
        if (!(alignment.IsFirstMate() && alignment.IsMateMapped())) {
          continue;
        }

        array[i].both_mapped++;

        // count proper pairs, although we have our own definition because
        // not all aligners use the same definition
        if (alignment.IsProperPair()) {
          ++array[i].properpair;
        }

        // mates must align to same contig, otherwise we record a bridge
        if (alignment.RefID != alignment.MateRefID) {
          array[i].bridges++;
          continue;
        }

        // fragment length must be plausible
        ldist = max(alignment.Position-alignment.MatePosition,
                    alignment.MatePosition-alignment.Position);
        if (ldist > realistic_distance) {
          // mates are too far apart
          continue;
        }

        // read orientation must match the generated library
        // in this case we only test for FR/RF orientation,
        // that is - we expect mates to be on opposite strands
        bool is_reversed = alignment.IsReverseStrand();
        bool is_mate_reversed = alignment.IsMateReverseStrand();

        if (!is_reversed && is_mate_reversed) {
          // in FR orientation, first read must start
          // before second read
          if (alignment.Position < alignment.MatePosition) {
            array[i].good++;
          }
        } else if (is_reversed && !is_mate_reversed) {
          // in RF orientation, second read must start
          // before first read
          if (alignment.MatePosition < alignment.Position) {
            array[i].good++;
          }
        }
      }
      array[i].bases_uncovered = pileup.getBasesUncovered();
      array[i].p_unique = pileup.getUniqueBases();
      array[i].p_not_segmented = pileup.p_not_segmented();

      if (reader.IsOpen()) {
        reader.Close();
      }
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
  realistic_distance = 450;
}

int main (int argc, char* argv[]) {
  if (argc >= 3) {
    if (argc == 4) {
      // user has supplied a segmentation null prior
      nullprior = atof(argv[3]);
    }
    string infile = argv[1];
    BetterBam bam (infile);
    bam.load_bam();
    int i;

    std::ofstream output;
    output.open (argv[2]);
    output << "name,p_seq_true,bridges,length,fragments_mapped,"
              "both_mapped,properpair,good,bases_uncovered,p_unique,"
              "p_not_segmented" << endl;
    for (i = 0; i < bam.get_seq_count(); i++) {
      output << bam.get_info(i).name << ",";
      if (bam.get_info(i).bases_mapped > 0) {
        output <<
        1-((double)bam.get_info(i).p_seq_true/bam.get_info(i).bases_mapped) <<
        ",";
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
    cout << "bam-read version 0.3.6\n"
         << "Usage:\n"
         << "bam-read <bam_file> <output_csv> <nullprior (optional)>\n\n"
         << "example:\n"
         << "bam-read in.bam out.csv 0.95"
         << endl;
    return 1;
  }
}
