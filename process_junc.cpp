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


void process_junc(vector<string> &instancelist, string &outfilename, string &regoutfilename)
{
	ofstream outfile;
	outfile.open(outfilename.c_str());
	ofstream regoutfile;
	regoutfile.open(regoutfilename.c_str());
	string line = "";
	int NofSample = 0;
	map<string,set<double> > junc;
	map<string,set<vector<int> > > regs;
	for (int i = 0; i < instancelist.size(); i++)
	{
		stringstream linestream;
		line = instancelist[i];
		linestream << line << ".instance";
		line = linestream.str();	
	
		ifstream s_infile;
		s_infile.open(line.c_str());
		
		string s_line = "";
		if (s_infile.is_open())
		{
			while(s_infile.good())
			{
				getline(s_infile,s_line);
				
				if (s_line[0] == 'I')
				{
					getline(s_infile,s_line);
					//cout << "line: " << count++ << endl;
					if (s_line[0] == 'B' and s_line[s_line.size()-1] != '.')
					{
						string strand = s_line.substr(s_line.size()-1,1);
						istringstream boundstream(s_line);
						string chr;
						getline(boundstream,chr,'\t');
						getline(boundstream,chr,'\t');
						string s_start_str,s_end_str;
						getline(boundstream,s_start_str,'\t');
						getline(boundstream,s_end_str,'\t');
						int s_start = atoi(s_start_str.c_str());
						int s_end = atoi(s_end_str.c_str());
						vector<int> s_reg;
						s_reg.push_back(s_start);
						s_reg.push_back(s_end);
						string chr1 = chr;
						chr1.append(strand);
						regs[chr1].insert(s_reg);

						vector<vector<double> > potential_al;
						getline(s_infile,s_line);
						if (s_line[0] == 'R')
						{
							getline(s_infile,s_line);
							if (s_line[0] == 'S')
							{
								int idx = s_line.find_first_of("\t");
								assert(idx > -1);
								int Nsegs = atoi(s_line.substr(idx+1).c_str());
								int count1 = 0;
								for (int i = 0; i < Nsegs; i++)
								{
									getline(s_infile,s_line);
									istringstream ss(s_line);
									vector<double> token_v;
									string token;
									while(getline(ss,token,'\t'))
										token_v.push_back(atof(token.c_str()));
									assert(token_v.size()==9);
									if (token_v[7] < 0.95)
									{
										if (count1 == 0)
											junc[chr].insert(token_v[0]-1);
										else
											junc[chr].insert(token_v[0]);
										junc[chr].insert(token_v[1]+1);
										count1++;
									}
									else
									{
										vector<double> s_al;
										s_al.push_back(i);
										if (count1 == 0)
											s_al.push_back(token_v[0]-1);
										else
											s_al.push_back(token_v[0]);
										s_al.push_back(token_v[1]+1);
										potential_al.push_back(s_al);
									}

								}
								if (potential_al.size() > 0)
								{
									getline(s_infile,s_line);getline(s_infile,s_line);getline(s_infile,s_line);
									assert(s_line[0] == 'S');
									int idx = s_line.find_first_of("\t");
									assert(idx > -1);
									int NSGT = atoi(s_line.substr(idx+1).c_str());
									vector<vector<double> >  SGType;
									for (int i = 0; i < NSGT; i++)
									{
										getline(s_infile,s_line);
										idx = s_line.find_first_of("\t");
										assert(idx > -1);
										string s_line1 = s_line.substr(0,idx);
										idx = s_line1.find_last_of(" ");
										s_line1 = s_line1.substr(0,idx);
										istringstream ss(s_line1);
										vector<double> token_v;
										string token;
										while(getline(ss,token,' '))
											token_v.push_back(atof(token.c_str()));
										SGType.push_back(token_v);
									}
									//check all s_al
									for (int i = 0; i < potential_al.size(); i++)
									{
										vector<double> s_al = potential_al[i];
										int pos = s_al[0];
										//check if pos is an end
										int large = 0;
										int small = 0;
										for (int j = 0; j < SGType.size(); j++)
										{
											if (SGType[j][pos] == 1)
											{
												for (int jj = 0; jj < SGType[j].size(); jj++)
												{
													if (SGType[j][jj] == 1)
													{
														if (jj < pos)
															small = 1;
														if (jj > pos)
															large = 1;
													}
												}
											}
										}
										if (small == 1 and large == 0)
										{
											junc[chr].insert(s_al[1]);
											double l = 20;
											if (s_al[2] - s_al[1] < l)
												l = s_al[2] - s_al[1] - 1;
											junc[chr].insert(s_al[1]+l);
										}

										if (small == 0 and large == 1)
										{
											junc[chr].insert(s_al[2]);
											double l = 20;
											if (s_al[2] - s_al[1] < l)
												l = s_al[2] - s_al[1] - 1;
											junc[chr].insert(s_al[2]-l);
										}
									}
								}
							}
								
						}
				
					}
				}
			}
		}
		
	}
	
	for (map<string,set<double> >::iterator it = junc.begin(); it != junc.end(); it++)
	{
		for (set<double>::iterator itset = it->second.begin(); itset != it->second.end(); itset++)
			outfile << fixed << it->first << "\t" << (int)*itset << "\t" << (int)*itset << "\t-" << endl;
	}
	
	int countins = 0;
	for (map<string,set<vector<int> > >::iterator it = regs.begin(); it != regs.end(); it++)
	{
		for (set<vector<int> >::iterator itset = it->second.begin(); itset != it->second.end(); itset++)
		{
			int s_start = (*itset)[0];
			int s_end = (*itset)[1];
			regoutfile << it->first.substr(0,it->first.size()-1) << "\t" << s_start << "\t" << s_end << "\tInst_" << ++countins << "\t" << s_end-s_start << "\t" << it->first[it->first.size()-1] << "\t" << s_start << "\t" << s_end <<"\t255,0,0\t1\t" << s_end-s_start << "\t0" << endl;
		}
	}

	outfile.close();
	regoutfile.close();
}

