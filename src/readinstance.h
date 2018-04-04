#ifndef READINSTANCE_H
#define READINSTANCE_H

#include "utility.h"
#include "options.h"
using namespace std;



class readinstance
{
public:
	readinstance(){};
	readinstance(Options &in_opt)
	{
		opt = in_opt;
	};
	Options opt;
	vector<Info> infolist;
	void readInstance_p(string inputfile);
	void readinstance_p(string inputfile, string outputfile);
	void getInfoList(vector<Info>& outinfolist);
	int totalNumofPaths;
	vector<double> TOTALNUMREADS;
	void runSparseIso_Threaded();
	int complicated;	
};


#endif
