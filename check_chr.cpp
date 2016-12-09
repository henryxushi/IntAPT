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

	//cout << outname1.str() << endl << outname2.str() << endl;
	//int aa;
	//cin >> aa;
	ifstream infile1;
	infile1.open(boundfile.c_str());

	string line;
	int count = 1;
	if (infile1.is_open())
	{
		while (infile1.good())
		{
			getline(infile1,line);
			int idx = line.find_first_of("\t");
			string chr = line.substr(0,idx-1);
			if (chr == "")
			{
				cout << "!!!" << endl;

				cout << count << endl;
				cout << line << endl;

			}
			cout << count << endl;
			count++;
			
		}
	}
	infile1.close();
	
	


	return 0;
}
