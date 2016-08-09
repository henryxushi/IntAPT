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
	if (infile1.is_open())
	{
		while(infile1.good())
		{
			getline(infile1,line);
			ifstream s_infile;
			s_infile.open(line.c_str());
			stringstream outstream;
			outstream << line << ".filtered.instance";
			ofstream s_outfile;
			s_outfile.open(outstream.str().c_str());
			string s_line = "", s_line1 = "";
			bool outbool = false;
			cout << "Start " << line << endl;
			int count = 0;
			if (s_infile.is_open())
			{
				while(s_infile.good())
				{
					getline(s_infile,s_line);
					//cout << "line: " << count++ << endl;
					if (s_line[0] == 'I')
					{
						outbool = false;
						getline(s_infile,s_line1);
						//cout << "line: " << count++ << endl;
						if (s_line1[0] == 'B' and s_line1[s_line1.size()-1] != '.')
						{
							stringstream linetemp;
							linetemp << s_line << endl << s_line1;
							s_line = linetemp.str();
							outbool = true;
						}
					}
					if (outbool)
						s_outfile << s_line << endl;
				}
			}
			cout << "Finished " << line << endl;
			//s_infile.close();
			//s_outfile.close();
		}
		
	}
	else
		exit(-1);
	
	infile1.close();
	return 0;
}

