#ifndef FILTER_BAM_H
#define FILTER_BAM_H
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <api/BamAux.h>
#include <api/BamWriter.h>
#include <queue>
#include <boost/thread.hpp>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <string.h>
#include <utility>
#include <map>
using namespace std;

struct ReadId
{            
    string name;
    uint hi_tag;

    bool operator == (const ReadId & ri) const { 
        return (name == ri.name and hi_tag == ri.hi_tag);
    }
};

struct ReadIdHasher
{
    size_t operator()(const ReadId & ri) const {

        return ((boost::hash<string>()(ri.name) ^ (boost::hash<uint>()(ri.hi_tag) << 1)) >> 1);
    }
};

struct PairInfo {
    
    uint pos;
    uint insert; // e_mpos
    string cigar;
    string m_cigar;

    bool operator == (const PairInfo & pi) const { 

        return (pos == pi.pos and insert == pi.insert and cigar == pi.cigar and m_cigar == pi.m_cigar);
    }
};

struct PairInfoHasher {
    
    size_t operator()(const PairInfo & pi) const {

        return (((((boost::hash<uint>()(pi.pos) ^ (boost::hash<uint>()(pi.insert) << 1)) >> 1) ^ (boost::hash<string>()(pi.cigar) >> 1)) << 1) ^ (boost::hash<string>()(pi.m_cigar) >> 1));
    }
};


typedef pair<BamTools::BamAlignment*, BamTools::BamAlignment*> ReadPair;
typedef tr1::unordered_map <ReadId, BamTools::BamAlignment*, ReadIdHasher> ReadIDs;

typedef map <uint, tr1::unordered_map <ReadId, BamTools::BamAlignment*, ReadIdHasher> > FirstReads;
typedef tr1::unordered_map <PairInfo, ReadPair, PairInfoHasher > ReadPairs;     

typedef multimap <uint, BamTools::BamAlignment*> UniqueReads;

class BAMfilter
{
public:
    BAMfilter(){};
    BAMfilter(string &inbam, string &inoutdir, int &inidx)
    {
        bamfile = inbam;
        outdir = inoutdir;
        //outbamfile = inoutbamfile;
        idx = inidx;
    };
    string bamfile;
    string outdir;
    int idx;
    void filter_bam();
    void markDuplicates(BamTools::BamAlignment & current_alignment, FirstReads * first_reads, ReadPairs * read_pairs);
    void addUniqueReads(UniqueReads *, ReadPairs *);
    double writeUniqueReads(BamTools::BamWriter *, UniqueReads *, FirstReads *);
    void generateBamIndex(string);
    string generateCigarString(vector<BamTools::CigarOp> &);
};

#endif