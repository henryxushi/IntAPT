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
	string boundfile = argv[1];
	stringstream bound_str;
	bound_str << boundfile << "_filtered.bed";
	string outfilename = bound_str.str();

	//cout << outname1.str() << endl << outname2.str() << endl;
	//int aa;
	//cin >> aa;
	ifstream infile1;
	infile1.open(boundfile.c_str());

	string curchr = "";
	int curs1 = 0, curs2 = 0;


	string line;
	int count = 0;
	ofstream outfile(outfilename.c_str());
	if (infile1.is_open())
	{
		while (infile1.good())
		{
			getline(infile1,line);
			int idx = line.find_first_of("\t");
			int idx1 = line.find_first_of("\t",idx+1);
			int idx2 = line.find_first_of("\t",idx1+1);
			int idx3 = line.find_first_of("\t",idx2+1);
			int idx4 = line.find_first_of("\t",idx3+1);
			int idx5 = line.find_first_of("\t",idx4+1);
			string strand = line.substr(idx4+1,idx5-idx4);
			if (strand == ".")
				continue;
			string chrstr = line.substr(0,idx);
			string s1str = line.substr(idx+1,idx1-idx);
			string s2str = line.substr(idx1+1,idx2-idx1);
			int s1 = atoi(s1str.c_str());
			int s2 = atoi(s2str.c_str());
		

			if (chrstr.compare(curchr) == 0)
			{
				if (s1 >= curs1 and s2 <= curs2)
					continue;
				
			}
			string line1 = line;
			if (line1[line1.size()-1] == '\t')
				line1 = line1.substr(0,line1.size()-1);
			if (count++ > 0)
				outfile << line1 << endl;
			//cout << line << endl;
			//cout << chrstr << "\t" << s1 << "\t" << s2 << endl;
			//int aa;
			//cin >> aa;
			curchr = chrstr;
			curs1 = s1;
			curs2 = s2;
		}
	}
	infile1.close();
	
	


	return 0;
}
