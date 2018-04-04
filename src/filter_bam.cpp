#include <iostream>
#include <api/BamWriter.h>
#include "filter_bam.h"
#include <string>
using namespace std;

void BAMfilter::filter_bam()
{
	BamTools::BamReader reader;
	//string bamfile = argv[1];
	//string outfile = argv[2];
    //stringstream outfile_stream;
    //outfile_stream << outdir << "/IntAPT_filtered" << idx <<".bam";
    string outfile = outdir;//outfile_stream.str();
    string outbamfile = outfile;
	assert(reader.Open(bamfile));
	assert(reader.IsOpen());

    const BamTools::SamHeader sam_header = reader.GetHeader();
    const BamTools::RefVector references = reader.GetReferenceData();

    BamTools::BamWriter unstranded_writer;
    //cout << "bamfile: " << bamfile << endl;
    //cout << "outfile: " << outfile << endl;
    assert(unstranded_writer.Open(outfile, sam_header, references));

    int current_position = -1;
    int cur_reference = -1;
    BamTools::BamAlignment current_alignment;
    int total_mapped_pair_reads = 0;
    int num_reads_paired = 0;
    double nd_nm_mapped_paired_end_reads = 0;

    FirstReads first_reads;
    ReadPairs read_pairs;
    UniqueReads unique_reads;
    while (reader.GetNextAlignment(current_alignment))
    {
    	uint nh_tag;
        uint8_t test_int;

        if(!current_alignment.GetTag("NH", test_int))
        {
            continue;
        }

        nh_tag = test_int;
        if (nh_tag > 1) {

            continue;
        }

        if(!current_alignment.IsPaired())
            continue;
        if(current_alignment.IsFailedQC())
            continue;

        if (current_alignment.IsMapped() and current_alignment.IsPrimaryAlignment() and current_alignment.IsMateMapped()) {

            total_mapped_pair_reads += 1;

            // Check that read pair is OK
            if ((current_alignment.IsReverseStrand() != current_alignment.IsMateReverseStrand()) and (current_alignment.RefID == current_alignment.MateRefID)) {

                // Same threshold as CEM
                if (abs(current_alignment.InsertSize) > 700000) {

                    continue;
                }

                unstranded_writer.SaveAlignment(current_alignment);
            }
        }
    }
    /*cout << "total_mapped_pair_reads: " << total_mapped_pair_reads << endl;
    cout << "nd_nm_mapped_paired_end_reads: " << nd_nm_mapped_paired_end_reads << endl;
    cout << "num_reads_paired: " << num_reads_paired << endl;*/

}

void BAMfilter::addUniqueReads(UniqueReads * unique_reads, ReadPairs * read_pairs) {

    for (ReadPairs::iterator it = read_pairs->begin(); it != read_pairs->end(); it++) {

        unique_reads->insert(unique_reads->end(), pair<uint, BamTools::BamAlignment*>(it->second.first->Position, it->second.first));
        unique_reads->insert(unique_reads->end(), pair<uint, BamTools::BamAlignment*>(it->second.second->Position, it->second.second));
    }

    read_pairs->clear();
}


/*void generateBamIndex(string bam_file_name) {

    stringstream idx_system_stream;
    idx_system_stream << "nice " << options_variables.samtools_path << " index " << bam_file_name;
    
    string idx_system_string = idx_system_stream.str();
    
    if (system(NULL) == 0) {
    
        cerr << "Command processor is not avaliable" << endl;
        exit(-1);
    
    } else {
    
        assert(system(idx_system_string.c_str()) == 0);
        assert(boost::filesystem::exists(bam_file_name));
    }
}*/

string BAMfilter::generateCigarString(vector<BamTools::CigarOp> & cigar_data) {

    stringstream read_cigar; 

    for (vector<BamTools::CigarOp>::iterator it = cigar_data.begin(); it != cigar_data.end(); it++) {

        if ((it->Type == 'M') or (it->Type == 'D')) {

            read_cigar << it->Length << "M";

        } else if (it->Type == 'N') {

            read_cigar << it->Length << "N";
        
        } else if (it->Type == 'S') {

            read_cigar << it->Length << "S";
        
        } else if (!(it->Type == 'I')) {

            cerr << "ERROR: Unhandled cigar string symbol '" << it->Type << "'!" << endl;
            exit(-1);
        }
   }

   return read_cigar.str();
}


double BAMfilter::writeUniqueReads(BamTools::BamWriter * writer, UniqueReads * unique_reads, FirstReads * first_reads) {

    double nd_nm_pe_reads = 0;

    uint cur_nh_tag;
    uint cur_hi_tag;

    UniqueReads::iterator uit = unique_reads->begin();

    uint left_most_position;

    if (first_reads->empty()) {

        left_most_position = numeric_limits<uint>::max();
    
    } else {

        left_most_position = first_reads->begin()->first;
    }
         
    while (uit->first < left_most_position and !unique_reads->empty()) {

        cur_hi_tag = 0;
        uint8_t test_int;
        assert(uit->second->GetTag("NH", test_int));
        cur_nh_tag = test_int;
        assert(uit->second->GetTag("HI", test_int) or (cur_nh_tag == 1));
        cur_hi_tag = test_int;

        if (cur_nh_tag == 1) {

            nd_nm_pe_reads += 1;
        }

        assert(writer->SaveAlignment(*(uit->second)));

        delete uit->second;
        unique_reads->erase(uit++);
    }
         
    return nd_nm_pe_reads; 
}


void BAMfilter::markDuplicates(BamTools::BamAlignment & current_alignment, FirstReads * first_reads, ReadPairs * read_pairs) 
{

    uint cur_hi_tag = 0;
    uint8_t test_int;
    current_alignment.GetTag("HI", test_int);
    cur_hi_tag = test_int;
    ReadId ri;
    ri.name = current_alignment.Name;
    ri.hi_tag = cur_hi_tag;

    FirstReads::iterator first_reads_it = first_reads->find(current_alignment.MatePosition);

    if (first_reads_it == first_reads->end()) {

        FirstReads::iterator cur_pos_first_reads_it = first_reads->find(current_alignment.Position);

        if (cur_pos_first_reads_it == first_reads->end()) {

            ReadIDs temp_read_ids;
            assert(temp_read_ids.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);
            assert(first_reads->insert(pair<uint, ReadIDs>(current_alignment.Position, temp_read_ids)).second);

        } else {

            assert(cur_pos_first_reads_it->second.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);
        }

    } else { 
        
        ReadIDs::iterator pos_first_reads_it = first_reads_it->second.find(ri);

        if (pos_first_reads_it == first_reads_it->second.end()) {

            //assert(current_alignment.Position == current_alignment.MatePosition);
            assert(first_reads_it->second.insert(pair<ReadId, BamTools::BamAlignment*>(ri, new BamTools::BamAlignment(current_alignment))).second);            
        
        } else {

            PairInfo pi;
            pi.pos = pos_first_reads_it->second->Position;
            pi.insert = abs(pos_first_reads_it->second->InsertSize);

            pi.cigar = generateCigarString(pos_first_reads_it->second->CigarData);
            pi.m_cigar = generateCigarString(current_alignment.CigarData);

            assert(current_alignment.MatePosition == pi.pos);

            ReadPairs::iterator read_pairs_it = read_pairs->find(pi);

            if (read_pairs_it == read_pairs->end()) {

                assert(read_pairs->insert(pair<PairInfo, ReadPair>(pi, ReadPair(pos_first_reads_it->second, new BamTools::BamAlignment(current_alignment)))).second);
            
            } else {

                uint cur_first_nh_tag;
                uint8_t test_int;
                assert(pos_first_reads_it->second->GetTag("NH", test_int));
                cur_first_nh_tag = test_int;
        
                uint cur_second_nh_tag;
                assert(current_alignment.GetTag("NH", test_int));
                cur_second_nh_tag = test_int;
        
                assert(cur_first_nh_tag == cur_second_nh_tag);

                uint prev_nh_tag;
                assert(read_pairs_it->second.first->GetTag("NH", test_int));
                prev_nh_tag = test_int;

                if (cur_first_nh_tag == 1 and prev_nh_tag > 1) {

                    delete read_pairs_it->second.first;
                    delete read_pairs_it->second.second;
                    
                    read_pairs_it->second = ReadPair(pos_first_reads_it->second, new BamTools::BamAlignment(current_alignment));
                
                } else {

                    delete pos_first_reads_it->second;
                }
            }
            
            first_reads_it->second.erase(pos_first_reads_it);

            if (first_reads_it->second.empty()) {

                first_reads->erase(first_reads_it);
            }
        }
    }
}
