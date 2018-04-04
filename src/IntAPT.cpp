#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <list>
#include <math.h>
#include <algorithm>
#include <vector>
#include <sstream>
#include <set>
#include <map>
#include <limits.h>
#include "utility.h"
#include "readinstance.h"
#include "options.h"
#include "Info.h"
#include "filter_bam.h"
#include <boost/threadpool.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include "merge_instances.h"
#include "process_junc.h"
#include <assert.h>
#include "mergerange.h"
using namespace std;





class SAMprocess
{
public:
	SAMprocess(){};
	SAMprocess(string &inbamfile, string &inoutdir, int &inidx, string &ininstname, string &inreadtype, string &intoolpath)
	{
		bamfile = inbamfile;
		outdir = inoutdir;
		idx = inidx;
		instname = ininstname;
		readtype = inreadtype;
		toolpath = intoolpath;
	}

	SAMprocess(string &inbamfile, string &inoutdir, int &inidx, string &ininstname, string &inbedsfile, string &inboundsfile, string &intoolpath)
	{
		bamfile = inbamfile;
		outdir = inoutdir;
		idx = inidx;
		instname = ininstname;
		bedsfile = inbedsfile;
		boundsfile = inboundsfile;
		toolpath = intoolpath;
	}

	void run()
	{
		/*if (readtype == "p")
        {
            stringstream outbamfilesorteds;
            outbamfilesorteds << toolpath << "samtools sort " << bamfile << " " << bamfile <<".sorted";
            assert(system(outbamfilesorteds.str().c_str())==0);
            stringstream outbamfilesorted_str;
            outbamfilesorted_str << bamfile << ".sorted.bam";
            string outbamfilesorted = outbamfilesorted_str.str();
            bamfile = outbamfilesorted;
        }*/
		stringstream samtools_stream, instance_stream;

		stringstream samfile;
		samfile << outdir << "/IntAPT_sam" << idx << ".sam";

		samtools_stream << toolpath << "samtools view " << bamfile << " > " << samfile.str();
		if (system(NULL) == 0)
		{
			cerr << "Processor not available" << endl;
			exit(-1);
		}
		else
			assert(system(samtools_stream.str().c_str())==0);

		instance_stream << "nice " << toolpath << "processsamS -g 10 -k 100000 -c 4 -i -o " << instname << " " << samfile.str() << " 1>/dev/null "; 
		
		if (system(NULL) == 0)
		{
			cerr << "Processor not available" << endl;
			exit(-1);
		}
		else
			assert(system(instance_stream.str().c_str())==0);
	};

	void run2nd()
	{
		stringstream instance_stream;

		stringstream samfile;
		samfile << outdir << "/IntAPT_sam" << idx << ".sam";

		instance_stream << "nice " << toolpath << "processsamS -g 10 -k 100000 -c 4 -i -o " << instname << " -r " << bedsfile << " -e " << boundsfile << " " << samfile.str() << " 1>/dev/null "; 
		if (system(NULL) == 0)
		{
			cerr << "Processor not available" << endl;
			exit(-1);
		}
		else
			assert(system(instance_stream.str().c_str())==0);
		stringstream inputinstfilestream, intronfilestream;
		inputinstfilestream << instname << ".instance";
		string inputinstfile = inputinstfilestream.str();
		intronfilestream << outdir <<"/intronlist" << idx;
		string intronfile = intronfilestream.str();
		filter_2nd(inputinstfile,intronfile);
	};

	string bamfile;
	string outdir;
	int idx;
	string instname;
	string bedsfile;
	string boundsfile;
	string readtype;
	string toolpath;
};

class Job
{
public:
	Job(Info* info_in, int idx_in) 
	{
		info = info_in;
		idx = idx_in;
	}
	~Job() {}
	void run()
	{
		if (info->single)
		{
			info->run_single();
		}
		else
		{
			info->run();
		}
	}
private:
	Info* info;
	int idx;
};



int main(int argc, char* argv[])
{
	//get info list in read instance class
	int complicated = 1;
	Options opt;
	if (opt.parse_options(argc,argv))
		exit(-1);

	if (!boost::filesystem::exists(opt.outdir))
		boost::filesystem::create_directory(opt.outdir);
	string finalInstFile;
	if (!opt.use_inst)
	{
		string outbamfile = opt.bamfile;

		if (!boost::filesystem::exists(opt.bamfile))
		{
			cerr << "The input bamfile list doesn't exist" << endl;
			exit(-1);
		}

		ifstream infile;
		infile.open(outbamfile.c_str());
		string line;
		vector<string> bamlist;
		if (infile.is_open())
		{
			while (infile.good())
			{
				getline(infile,line);
				if (!line.empty())
					bamlist.push_back(line);
			}
		}
		else
		{
			cerr << "Failed to open input file" << endl;
			exit(-1);
		}
	 	int J = bamlist.size();
	 	cout << "Found " << J << " input bam files in the provided list..." << endl;
		vector<string> bamlist_filtered;

		// filter reads for paired end data
		
		if (opt.readtype == "p")
		{
			for (int i = 0; i < bamlist.size(); i++)
			{
				stringstream outfile_stream;
		    	outfile_stream << opt.outdir << "/IntAPT_filtered" << i <<".bam";
		    	bamlist_filtered.push_back(outfile_stream.str());
			}
			cout << "Input data is paired-end data ..." << endl;
			
			if (opt.N_cores > 1)
			{
				boost::threadpool::pool multi_tp(opt.N_cores);
				for (int i = 0; i < bamlist.size(); i++)
				{
					boost::shared_ptr<BAMfilter> filter(new BAMfilter(bamlist[i],bamlist_filtered[i],i));
					multi_tp.schedule(boost::bind(&BAMfilter::filter_bam,filter));
				}
				multi_tp.wait();
			}
			else
			{
				for (int i = 0; i < bamlist.size(); i++)
				{
					BAMfilter filter(bamlist[i],bamlist_filtered[i],i);
					filter.filter_bam();
				}
			}
		}
		else
			bamlist_filtered = bamlist;

		cout << "Processing bam file" << endl;
		if (system(NULL) == 0)
		{
			cerr << "Processor not available" << endl;
			exit(-1);
		}

		vector<string> first_instlist;

		for (int i = 0; i < bamlist_filtered.size(); i++)
		{
			stringstream outinstance_stream;
			outinstance_stream << opt.outdir << "/sample_" << i << ".1st";
			first_instlist.push_back(outinstance_stream.str());
		}

		if (opt.N_cores > 1)
		{
			boost::threadpool::pool multi_tp(opt.N_cores);
			for (int i = 0; i < J; i++)
			{
				boost::shared_ptr<SAMprocess> samproc(new SAMprocess(bamlist_filtered[i],opt.outdir,i,first_instlist[i],opt.readtype,opt.cempath));
				multi_tp.schedule(boost::bind(&SAMprocess::run,samproc));
			}
			multi_tp.wait();
		}
		else
		{
			for (int i = 0; i < J; i++)
			{
				SAMprocess samproc(bamlist_filtered[i],opt.outdir,i,first_instlist[i],opt.readtype,opt.cempath);
				samproc.run();
			}
		}

		cout << "Jointly analyzing the junction..." << endl;
		stringstream boundsfilestream, bedsfilestream, intronliststream;
		boundsfilestream << opt.outdir << "/allbounds.txt";
		bedsfilestream << opt.outdir << "/allbeds.bed";
		intronliststream << opt.outdir << "/intronlist";
		string boundsfile = boundsfilestream.str();
		string bedsfile = bedsfilestream.str();
		string intronfileprefix = intronliststream.str();
		process_junc(first_instlist,boundsfile,bedsfile, intronfileprefix);

		stringstream bedssortedfilestream;
		bedssortedfilestream << opt.outdir << "/allbeds_sorted.bed";
		string bedssortedfile = bedssortedfilestream.str();
		stringstream sortbedcmd;
		sortbedcmd << "sort -k 1,1 -k 2,2n " << bedsfile << " > " << bedssortedfile;
		assert(system(sortbedcmd.str().c_str())==0);

		//filter old MATLAB code
		//fitlering large regions that highly possible to be false positives, code in process_junc
		stringstream filteredbedstream;
		filteredbedstream << opt.outdir << "/allbeds_sorted_filtered.bed";
		string filteredbed = filteredbedstream.str();
		filter_bed(bedssortedfile,filteredbed);
		
		stringstream filteredsortedbedstream;
		filteredsortedbedstream << opt.outdir << "/allbeds_sorted_filtered_sorted.bed";
		string filteredsortedbed = filteredsortedbedstream.str();
		stringstream sortfilteredbedcmd;
		sortfilteredbedcmd << "sort -k 1,1 -k 2,2n " << filteredbed << " > " << filteredsortedbed; 
		assert(system(sortfilteredbedcmd.str().c_str())==0);

		stringstream mergebedstream;
		mergebedstream << opt.outdir << "/merged_sorted.bed";
		string mergebed_str = mergebedstream.str();

		mergebed(filteredsortedbed,mergebed_str);
		

		// test if graph is complex
		ifstream ifcountjunc,ifcountregion;
		double countjunc = 0, countregion = 0;
		ifcountjunc.open(boundsfile.c_str());
		string countline = "";
		if (ifcountjunc.is_open())
		{
			while (ifcountjunc.good())
			{
				getline(ifcountjunc,countline);
				countjunc++;
			}
		}
		ifcountjunc.close();
		
		ifcountregion.open(mergebed_str.c_str());
                countline = "";
                if (ifcountregion.is_open())
                {
                        while (ifcountregion.good())
                        {
                                getline(ifcountregion,countline);
                                countregion++;
                        }
                }
		cout << "Total num of junctions: " << countjunc << endl;
		cout << "Total num of regions: " << countregion << endl;
		if (countjunc/countregion>35)
			complicated = 1;
		else
			complicated = 0;
		
		// 2nd round 
		cout << "Re-analyzing the junctions" << endl;
		vector<string> second_instlist;
		for (int i = 0; i < bamlist_filtered.size(); i++)
		{
			stringstream outinstance_stream;
			outinstance_stream << opt.outdir << "/sample_" << i << ".2nd";
			second_instlist.push_back(outinstance_stream.str());
		}

		vector<string> second_instlist_filtered;
		for (int i = 0; i < bamlist_filtered.size(); i++)
		{
			stringstream outinstance_stream;
			outinstance_stream << opt.outdir << "/sample_" << i << ".2nd.filtered";
			second_instlist_filtered.push_back(outinstance_stream.str());
		}
		
		if (opt.readtype == "p")
		{
			for (int i = 0; i < bamlist_filtered.size();i++)
			{
				stringstream modbam;
				modbam << bamlist_filtered[i] << ".sorted.bam";
				bamlist_filtered[i] = modbam.str();
			}
		}
		if (opt.N_cores > 1) //2nd inst filtering is integrated in the multiprocess
		{
			boost::threadpool::pool multi_tp(opt.N_cores);
			for (int i = 0; i < J; i++)
			{
				boost::shared_ptr<SAMprocess> samproc(new SAMprocess(bamlist_filtered[i],opt.outdir,i,second_instlist[i],mergebed_str,boundsfile,opt.cempath));
				multi_tp.schedule(boost::bind(&SAMprocess::run2nd,samproc));
			}
			multi_tp.wait();
		}
		else
		{
			for (int i = 0; i < J; i++)
			{
				SAMprocess samproc(bamlist_filtered[i],opt.outdir,i,second_instlist[i],mergebed_str,boundsfile,opt.cempath);
				samproc.run2nd();
			}
		}

		// merge all ranges
		stringstream mergeinststream;
		mergeinststream << opt.outdir << "/IntAPT_proc.ins";
		string mergeinst = mergeinststream.str();
		MergeInst mergeinst_c(second_instlist_filtered,mergeinst);
		mergeinst_c.merge_instances();
		finalInstFile = mergeinst;
	}
	else
	{
		finalInstFile = opt.InstFile;
		complicated = 1;
	}

	if (opt.InstOnly)
		return 0;


	//start from inst file
	stringstream outputfile_stream;
	outputfile_stream << opt.outdir << "/IntAPT_enum_inst";
	string outputfile = outputfile_stream.str();
	cout << "Constructing splicing graph" << endl;
	readinstance loadData(opt);
	loadData.complicated = complicated;
	loadData.readinstance_p(finalInstFile,outputfile);
	cout << "Finished reading" << endl;
	ofstream outfile1;
	outfile1.open(outputfile.c_str(),std::ios_base::app);
	outfile1 << "totalNumReads:\t" << endl;
	for (int i = 0; i < loadData.TOTALNUMREADS.size(); i++)
	{
		outfile1 << loadData.TOTALNUMREADS[i] << "\t";
	}
	outfile1 << endl;

	vector<int> totalReads;
	double sumreads = 0;
	for (int i = 0; i < loadData.TOTALNUMREADS.size(); i++)
	{
		totalReads.push_back(round(loadData.TOTALNUMREADS[i]/2));
		sumreads += loadData.TOTALNUMREADS[i]/2;
	}
	
	cout << "Start sampling process\n" << endl;

	if (opt.N_cores > 1)
	{
		boost::threadpool::pool multi_tp(opt.N_cores);
		for (int i = 0; i < loadData.infolist.size(); i++)
		{
			if (!loadData.infolist[i].valid)
				continue;
			if (loadData.infolist[i].valid)
			{
				loadData.infolist[i].infoidx = i;
				loadData.infolist[i].TOTALREADSJ = totalReads;
				loadData.infolist[i].totalReads= round(sumreads / loadData.TOTALNUMREADS.size());
				boost::shared_ptr<Job> job(new Job(&loadData.infolist[i],i));
				multi_tp.schedule(boost::bind(&Job::run,job));
			}
			
		}
		multi_tp.wait();
	}
	else
	{
		for (int i = 0; i < loadData.infolist.size(); i++)
		{
			if (!loadData.infolist[i].valid)
				continue;
			loadData.infolist[i].infoidx = i;
			loadData.infolist[i].TOTALREADSJ = totalReads;
			loadData.infolist[i].totalReads= round(sumreads / loadData.TOTALNUMREADS.size());
			
			if (loadData.infolist[i].single)
			{
				loadData.infolist[i].run_single();
			}
			else
			{
				loadData.infolist[i].run();
			}
		}
	}

	cout << "Start writing..." << endl;
	for (int i = 0; i < loadData.infolist.size(); i++)
	{
		if (!loadData.infolist[i].valid)
			continue;
		loadData.infolist[i].output_all_bool = opt.output_all;
		if (loadData.infolist[i].single)
		{
			if (loadData.infolist[i].output_bool)
			{
				loadData.infolist[i].write(opt.outdir,i+1);
			}
		}
		else
		{
			loadData.infolist[i].write(opt.outdir,i+1);
		}
		
	}

	return 0;
}
