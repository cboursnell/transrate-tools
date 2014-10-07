// ***************************************************************************
// bamtools_pileup_engine.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011
// ---------------------------------------------------------------------------
// Provides pileup at position functionality for various tools.
// ***************************************************************************

// #include "utils/utils_global.h"

#include <api/BamAlignment.h>
#include <vector>

using namespace BamTools;

// contains auxiliary data about a single BamAlignment
// at current position considered
struct PileupAlignment {

    // data members
    BamAlignment Alignment;
    int32_t PositionInAlignment;
    bool IsCurrentDeletion;
    bool IsNextDeletion;
    bool IsNextInsertion;
    int DeletionLength;
    int InsertionLength;
    bool IsSegmentBegin;
    bool IsSegmentEnd;

    // ctor
    PileupAlignment(const BamAlignment& al)
        : Alignment(al)
        , PositionInAlignment(-1)
        , IsCurrentDeletion(false)
        , IsNextDeletion(false)
        , IsNextInsertion(false)
        , DeletionLength(0)
        , InsertionLength(0)
        , IsSegmentBegin(false)
        , IsSegmentEnd(false)
    { }
};

// contains all data at a position
struct PileupPosition {

    // data members
    int RefId;
    int Position;
    std::vector<PileupAlignment> PileupAlignments;

    // ctor
    PileupPosition(const int& refId = 0,
                   const int& position = 0,
                   const std::vector<PileupAlignment>& alignments = std::vector<PileupAlignment>())
        : RefId(refId)
        , Position(position)
        , PileupAlignments(alignments)
    { }
};

class TransratePileup {

    public:
        TransratePileup(void);
        ~TransratePileup(void);

    public:
        bool AddAlignment(const BamAlignment& al, int ref_length);
        std::vector<int> getCoverage(void);
        void Flush(void);

    private:
        struct TransratePileupPrivate;
        TransratePileupPrivate* d;
};

