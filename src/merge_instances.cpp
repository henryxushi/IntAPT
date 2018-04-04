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
#include <assert.h>
#include "merge_instances.h"
using namespace std;




void MergeInst::merge_instances()
{
	vector<ifstream *> inst_stream;
	string line = "";
	int NofSample = 0;
	for (int i = 0; i < instancelist.size(); i++)
	{
		stringstream linestream;
		line = instancelist[i];
		linestream << line << ".instance";
		line = linestream.str();
		inst_stream.push_back(new ifstream());
		inst_stream[inst_stream.size()-1]->open(line.c_str());
		if(inst_stream[inst_stream.size()-1]->is_open())
		{
			NofSample++;
		}
		else
			inst_stream.pop_back();			
	}

	vector<int> totalReads;
	countTotalReads(inst_stream,NofSample,totalReads);

	inst_stream.clear();
	line = "";
	NofSample = 0;
	for (int i = 0; i < instancelist.size(); i++)
	{			
		stringstream linestream;
		line = instancelist[i];
		linestream << line << ".instance";
		line = linestream.str();
		inst_stream.push_back(new ifstream());
		inst_stream[inst_stream.size()-1]->open(line.c_str());
		if(inst_stream[inst_stream.size()-1]->is_open())
		{
			NofSample++;
		}
		else
			inst_stream.pop_back();			
		
	}

	ofstream outfile;
	outfile.open(outputfilename.c_str());
	int inst_count = 0;
	while(true)
	{
		string tag = getInstTag(inst_stream,NofSample);
		cout << tag << endl;
		if (tag.empty())
			break;
		outfile << "Instance\t" << ++inst_count << endl;


		vector<double> ReadLength;
		getReadLength(inst_stream,NofSample,ReadLength);
				
		
		vector<ExonInfo> exonInfo;
		getExonBound(inst_stream,NofSample,exonInfo);
		outfile << "segment\t" << exonInfo.size() << endl;
		for(int i = 0; i < exonInfo.size(); i++)
			outfile << "!" << i << "\t" << exonInfo[i].s1 << "\t" << exonInfo[i].s2 << endl;

		outfile << "segcount\t" << exonInfo.size() << endl;
		for(int i = 0; i < exonInfo.size(); i++)
		{
			outfile << "!" << i << "\t";
			for (int j = 0; j < exonInfo[i].NofReads.size();j++)
			{
				outfile << exonInfo[i].NofReads[j] << "\t";
			} 
			outfile << endl;
		}

		outfile << "segcvg\t" << exonInfo.size() << endl;
		for(int i = 0; i < exonInfo.size(); i++)
		{
			outfile << "!" << i << "\t";
			for (int j = 0; j < exonInfo[i].cov.size();j++)
			{
				outfile << exonInfo[i].cov[j] << "\t";
			} 
			outfile << endl;
		}

		outfile << "segp0\t" << exonInfo.size() << endl;
		for(int i = 0; i < exonInfo.size(); i++)
		{
			outfile << "!" << i << "\t";
			for (int j = 0; j < exonInfo[i].p.size();j++)
			{
				outfile << exonInfo[i].p[j] << "\t";
			} 
			outfile << endl;
		}

		outfile << "segmaxcvg\t" << exonInfo.size() << endl;
		for(int i = 0; i < exonInfo.size(); i++)
		{
			outfile << "!" << i << "\t";
			for (int j = 0; j < exonInfo[i].maxcov.size();j++)
			{
				outfile << exonInfo[i].maxcov[j] << "\t";
			} 
			outfile << endl;
		}


		map<vector<int> , vector<double> > segInfo;
		getSegInfo(inst_stream, NofSample, segInfo);
		outfile << "sgtype\t" << segInfo.size() << endl;


		int sgcount = 0;
		for (std::map<vector<int> , vector<double> >::iterator it = segInfo.begin();it != segInfo.end(); ++it)
		{
			outfile << "!" << sgcount++ << "\t";
			for(int i = 0; i < it->first.size(); i++)
			{
				outfile << it->first[i] << " ";
			}
			outfile << "0 " << endl;
		}

		sgcount = 0;
		outfile << "sgcount\t" << segInfo.size() << endl;
		for (std::map<vector<int> , vector<double> >::iterator it = segInfo.begin();it != segInfo.end(); ++it)
		{
			outfile << "!" << sgcount++ << "\t";
			for(int i = 0; i < it->second.size(); i++)
			{
				outfile << it->second[i] << "\t";
			}
			outfile << endl;
		}

		outfile << "petype\t" << 0 << endl;
		outfile << "pecount\t" << 0 << endl;
		outfile << "readlen\t" << 1 << endl;
		for (int i = 0; i < ReadLength.size(); i++)
		{
			outfile << ReadLength[i] << "\t";
		}
		outfile << endl;

		outfile << "totalreads\t" << 1 << endl;
		for (int i = 0; i < totalReads.size(); i++)
		{
			outfile << totalReads[i] << "\t";
		}
		outfile << endl;

		outfile << "otherinfo\t" << 2 << endl;
		int idx = tag.find_first_of("\t");
		outfile << "Chr\t" << tag.substr(0,idx) << endl;
		outfile << "Dir\t" << tag[tag.size()-1] << endl;
	}
	
}

void MergeInst::countTotalReads(vector<ifstream *> & inst_stream, int NofSample, vector<int> &totalReads)
{
	totalReads.resize(NofSample,0);
	for (int i = 0; i < NofSample; i++)
	{
		string line;
		while(inst_stream[i]->good())
		{
			getline(*inst_stream[i],line);
			if (line[0] != 'R')
				continue;
			int idx = line.find_first_of("\t");
			if (idx > 0)
			{
				string label = line.substr(0,idx);
				if (label.compare("Reads") == 0)
					totalReads[i] += atoi(line.substr(idx+1).c_str());
			}
		}
	}
}

string MergeInst::getInstTag(vector<ifstream *> & inst_stream, int NofSample)
{
	string tag = "";
	for (int i = 0; i < NofSample; i++)
	{
		string line;
		while(inst_stream[i]->good())
		{
			getline(*inst_stream[i],line);
			int idx = line.find_first_of("\t");
			string label = line.substr(0,idx);
		
			if (label.compare("Instance") == 0)
			{
				getline(*inst_stream[i],line);
				idx = line.find_first_of("\t");
				label = line.substr(0,idx);
				if (label.compare("Boundary") == 0)
				{
					if (i == 0)
						tag = line.substr(idx+1);
					else
					{
						if (tag != line.substr(idx+1))
						{
							cout << tag << " vs " << line.substr(idx+1) << endl;
						}
						assert(tag == line.substr(idx+1));
					}
					break;
				}
			}
		}
	}
	return tag;
}

void MergeInst::getReadLength(vector<ifstream *> & inst_stream, int NofSample, vector<double> &ReadLength)
{
	ReadLength.clear();
	for (int i = 0; i < NofSample; i++)
	{
		string line;
		while(inst_stream[i]->good())
		{
			getline(*inst_stream[i],line);
			int idx = line.find_first_of("\t");
			string label = line.substr(0,idx);
			assert(!label.compare("ReadLen"));
			ReadLength.push_back(atof(line.substr(idx+1).c_str()));
			break;
		}
	}
}

void MergeInst::getExonBound(vector<ifstream*> & inst_stream, int NofSample,vector<ExonInfo> &exonInfo)
{
	for (int i = 0; i < NofSample; i++)
	{
		string line;
		while(inst_stream[i]->good())
		{
			getline(*inst_stream[i],line);
			int idx = line.find_first_of("\t");
			string label = line.substr(0,idx);
			if (i == 0)
			{
				while (label.compare("Refs")==1 and inst_stream[i]->good())
				{
					getline(*inst_stream[i],line);
					int idx = line.find_first_of("\t");
					if (line[0] != 'R')
					{
						ExonInfo info;
						info.NofReads.clear();
						info.p.clear();
						info.s1 = atoi(line.substr(0,idx).c_str());
						int idx1 = line.find_first_of("\t",idx+1);
						info.s2 = atoi(line.substr(idx+1,idx1-idx).c_str());
						idx1 = line.find_first_of("\t",idx1+1);
						int idx2 = line.find_first_of("\t",idx1+1);
						info.NofReads.push_back(atof(line.substr(idx1+1,idx2-idx1).c_str()));
						int idx3 = line.find_first_of("\t",idx2+1);
						info.maxcov.push_back(atof(line.substr(idx2+1,idx3-idx2).c_str()));

						idx = line.find_last_of("\t");
						info.cov.push_back(atof(line.substr(idx+1).c_str()));
						idx1 = line.find_last_of("\t",idx-1);
						info.p.push_back(atof(line.substr(idx1+1,idx-idx1).c_str()));
						exonInfo.push_back(info);
					}
					else
						break;

				}
			}
			else
			{
				while(inst_stream[i]->good())
				{
					for (int j = 0; j < exonInfo.size();j++)
					{
						getline(*inst_stream[i],line);
						int idx = line.find_first_of("\t");	
						assert(exonInfo[j].s1 == atoi(line.substr(0,idx).c_str()));
						int idx1 = line.find_first_of("\t",idx+1);
						assert(exonInfo[j].s2 == atoi(line.substr(idx+1,idx1-idx).c_str()));
						idx1 = line.find_first_of("\t",idx1+1);
						int idx2 = line.find_first_of("\t",idx1+1);
						exonInfo[j].NofReads.push_back(atof(line.substr(idx1+1,idx2-idx1).c_str()));
						int idx3 = line.find_first_of("\t",idx2+1);
						exonInfo[j].maxcov.push_back(atof(line.substr(idx2+1,idx3-idx2).c_str()));

						idx = line.find_last_of("\t");
						exonInfo[j].cov.push_back(atof(line.substr(idx+1).c_str()));
						idx1 = line.find_last_of("\t",idx-1);
						exonInfo[j].p.push_back(atof(line.substr(idx1+1,idx-idx1).c_str()));
					}
					getline(*inst_stream[i],line);
					break;
				}
			}
			break;
		}
	}
}

void MergeInst::getSegInfo(vector<ifstream*> & inst_stream, int NofSample, map<vector<int> , vector<double> > &segInfo)
{
	for (int i = 0; i < NofSample; i++)
	{
		string line;
		while(inst_stream[i]->good())
		{
			getline(*inst_stream[i],line);
			getline(*inst_stream[i],line);
			while (true)
			{
				getline(*inst_stream[i],line);
				if (line[0] == 'P')
					break;
				int idx = line.find_first_of("\t");
				string conf_str = line.substr(0,idx);
				vector<int> conf;
				double Nreads;
				sep_space(conf_str,conf,Nreads);
				if (i == 0)
					segInfo[conf].push_back(Nreads);
				else
				{
					std::map<vector<int> , vector<double> >::iterator it = segInfo.find(conf);
					if (it == segInfo.end())
					{
						for(int j = 0; j < i; j++)
						{
							segInfo[conf].push_back(0);
						}
						segInfo[conf].push_back(Nreads);
					}
					else
					{
						it->second.push_back(Nreads);
					}
				}
			}
			for (std::map<vector<int> , vector<double> >::iterator it = segInfo.begin();it != segInfo.end(); ++it)
			{
				if (it->second.size()<i+1)
				{
					assert(it->second.size()==i);
					it->second.push_back(0);
				}
			} 	
			break;
		}
	}
}

void MergeInst::sep_space(string conf_str, vector<int> & conf, double &Nreads)
{
	int idx = conf_str.find_first_of(" ");
	int preidx = -1;
	while(idx > 0)
	{
		conf.push_back(atoi(conf_str.substr(preidx+1,idx-preidx-1).c_str()));
		preidx = idx;
		idx = conf_str.find_first_of(" ", preidx+1);
	}
	Nreads = atof(conf_str.substr(preidx+1).c_str());
}
