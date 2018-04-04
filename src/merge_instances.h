#ifndef MERGE_INSTANCES_H
#define MERGE_INSTANCES_H
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
using namespace std;

struct ExonInfo
{
	int s1;
	int s2;
	vector<double> p;
	vector<double> NofReads;
	vector<double> cov;
	vector<double> maxcov;

};

class MergeInst
{
public:
	MergeInst(){};
	MergeInst(vector<string> &ininstancelist, string &inoutputfilename)
	{
		instancelist = ininstancelist;
		outputfilename = inoutputfilename;
	};

	string getInstTag(vector<ifstream *> & inst_stream, int NofSample);
	void getReadLength(vector<ifstream *> & inst_stream, int NofSample, vector<double> &ReadLength);
	void getExonBound(vector<ifstream *> & inst_stream, int NofSample,vector<ExonInfo> &exonInfo);
	void print_v(vector<double> v);
	void print_v(vector<int> v);
	void print_vv(vector<vector<double> > v);
	void getSegInfo(vector<ifstream*> & inst_stream, int NofSample, map<vector<int> , vector<double> > &segInfo);
	void countTotalReads(vector<ifstream *> & inst_stream, int NofSample, vector<int> &totalReads);
	void sep_space(string conf_str, vector<int> & conf, double &Nreads);
	void merge_instances();
	vector<string> instancelist;
	string outputfilename;

};


#endif
