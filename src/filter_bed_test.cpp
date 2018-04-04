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
	//store the input bed into a map with chr1- as key
	string inputbed = argv[1];
	ifstream infile;
	infile.open(inputbed.c_str());
	string outputbed = argv[2];
	ofstream outfile;
	outfile.open(outputbed.c_str());

	string line = "";
	map<string,vector<vector<int> > > bedcontents;//chr->start,end,length
	if (infile.is_open())
	{
		while(infile.good())
		{
			getline(infile,line);
			istringstream ss(line);
			vector<string> contents;
			string temp;
			while(getline(ss,temp,'\t'))
				contents.push_back(temp);
			if(contents.size()!=12)
				continue;
			if(contents[5].length()!=1)
				continue;
			stringstream keystream;
			keystream << contents[0] << contents[5];
			string keychr = keystream.str();
			vector<int> insertitem; //start, end, length
			insertitem.push_back(atoi(contents[1].c_str()));
			insertitem.push_back(atoi(contents[2].c_str()));
			insertitem.push_back(insertitem[1]-insertitem[0]+1);
			bedcontents[keychr].push_back(insertitem);
		}
	}
	infile.close();
	//sliding window to test
	int K = 3;
	int F = 10;
	int instcount = 1;
	
	for(std::map<string,vector<vector<int> > >::iterator it = bedcontents.begin(); it != bedcontents.end(); it++)
	{
		vector<vector<int> > pos = it->second;
		if (pos.size()<10*K)
			continue;
		for (int posidx = 0; posidx < pos.size(); posidx++)
		{
			int startidx = posidx-K;
			if (startidx < 0)
				startidx = 0;
			int endidx = posidx+K;
			if (endidx >= pos.size()-1)
				endidx = pos.size()-1;

			int startsum = 0;
			int countsum = 0;
			for (int ii = startidx; ii < posidx; ii++)
			{
				startsum += pos[ii][2];
				countsum += 1;
			}
			double startmean = 0;
			if (countsum > 0)
				startmean = (double)startsum/(double)countsum;

			int endsum = 0;
			countsum = 0;
			for (int ii = posidx+1; ii <= endidx; ii++)
			{
				endsum += pos[ii][2];
				countsum += 1;
			}
			double endmean = 0;
			if (countsum > 0)
				endmean = (double)endsum/(double)countsum;

			if ((double)pos[posidx][2] > (double)F*startmean and (double)pos[posidx][2] > (double)F*endmean)
			{
				//skip
			}
			else
			{
				string chr = it->first.substr(0,it->first.size()-1);
				string strand = it->first.substr(it->first.size()-1,1);
				outfile << fixed << chr << "\t" << pos[posidx][0] << "\t" << pos[posidx][1] << "\tInst_" << instcount << "\t" << pos[posidx][2] << "\t" << strand << "\t" << pos[posidx][0] << "\t" << pos[posidx][1] << "\t255,0,0\t1" << pos[posidx][2] << "\t0" << endl;  
				instcount++;			
			}

		}
	}
	outfile.close();
}

