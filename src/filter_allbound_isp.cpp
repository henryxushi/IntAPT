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


int main(int argc, char *argv[])
{
	string instancelist = argv[1];
	ifstream infile1;
	infile1.open(instancelist.c_str());
	string line = "";
	int NofSample = 0;
	stringstream outstream;
	outstream << instancelist << ".filtered.bed";
	ofstream s_outfile;
	s_outfile.open(outstream.str().c_str());
	if (infile1.is_open())
	{
		while(infile1.good())
		{
			getline(infile1,line);
			int idx = line.find_first_of(".");
			if (idx != -1)
				continue;
			if (line.find("seg") != std::string::npos)
				continue;
			s_outfile << line << endl;
		}
		
	}
	else
		exit(-1);
	
	infile1.close();
	s_outfile.close();
	return 0;
}

