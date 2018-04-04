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
#include "process_junc.h"
using namespace std;


void process_junc(vector<string> &instancelist, string &outfilename, string &regoutfilename, string &intronfileprefix)
{
	ofstream outfile;
	outfile.open(outfilename.c_str());
	ofstream regoutfile;
	regoutfile.open(regoutfilename.c_str());
	string line = "";
	int NofSample = 0;
	map<string,set<double> > junc;
	map<string,set<vector<int> > > regs;
	for (int fileidx = 0; fileidx < instancelist.size(); fileidx++)
	{
		stringstream linestream;
		line = instancelist[fileidx];
		linestream << line << ".instance";
		line = linestream.str();	
	
		ifstream s_infile;
		s_infile.open(line.c_str());
		int current_intron_point = 1;
		int prev_intron_point = 1;
		string s_line = "";
		stringstream temp;
		temp << intronfileprefix << fileidx;
		ofstream intronoutfile;
		intronoutfile.open(temp.str().c_str());
		if (s_infile.is_open())
		{
			while(s_infile.good())
			{
				getline(s_infile,s_line);
				
				if (s_line[0] == 'I')
				{
					getline(s_infile,s_line);
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
								int intron_s1 = -1;
								int intron_s2 = -1;
								int count2 = 0;
								for (int i = 0; i < Nsegs; i++)
								{
									getline(s_infile,s_line);
									istringstream ss(s_line);
									vector<double> token_v;
									string token;
									while(getline(ss,token,'\t'))
										token_v.push_back(atof(token.c_str()));
									assert(token_v.size()==9);
									if (i > 0)// and count2 == 1)
									{
										intron_s2 = token_v[0]-1;
										if (intron_s2 >= intron_s1)
										{
											intronoutfile << fixed << chr << " " << intron_s1 << " " << intron_s2 << endl;
										}
										intron_s1 = token_v[1]+1;
										intron_s2 = -1;
									}
									else
									{
										intron_s1 = token_v[1]+1;
									}
									if (token_v[7] < 0.95)
									{
										if (count1 == 0)
										{
											current_intron_point = token_v[0]-1;
											intronoutfile << fixed << chr << " " << prev_intron_point << " " << current_intron_point << endl;
										}
										prev_intron_point = token_v[1]+1;

										if (count1 == 0)
											junc[chr].insert(token_v[0]-1);
										else
											junc[chr].insert(token_v[0]);
										junc[chr].insert(token_v[1]+1);
										count1++;
									}
									else
									{
										intronoutfile << fixed << chr << " " << (int)token_v[0] << " " << (int)token_v[1] << endl;
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
		s_infile.close();
		intronoutfile.close();
		
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

void filter_bed(string &inputbed, string &outputbed)
{
	//store the input bed into a map with chr1- as key
	ifstream infile;
	infile.open(inputbed.c_str());
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
				outfile << fixed << chr << "\t" << pos[posidx][0] << "\t" << pos[posidx][1] << "\tInst_" << instcount << "\t" << pos[posidx][2] << "\t" << strand << "\t" << pos[posidx][0] << "\t" << pos[posidx][1] << "\t255,0,0\t1\t" << pos[posidx][2] << "\t0" << endl;  
				instcount++;			
			}

		}
	}
	outfile.close();
}

void filter_2nd(string &inputinst, string &outputinst)
{
	string instancelist = inputinst;
	ifstream infile1;
	infile1.open(instancelist.c_str());


	string intronlistfile= outputinst;
	ifstream infile2;
	infile2.open(intronlistfile.c_str());
	string line = "";
	int NofSample = 0;

	map<string, map<int,int> > intron_list;
	//read intron list
	if (infile2.is_open())
	{
		while(infile2.good())
		{
			getline(infile2,line);
			int idx = line.find_first_of(" ");
			if (idx == -1)
				continue;
			istringstream ss(line);
			vector<string> token_v;
			string token;
			while(getline(ss,token,' '))
				token_v.push_back(token);
			intron_list[token_v[0]][atoi(token_v[1].c_str())] = atoi(token_v[2].c_str());
		}
	}
	else
	{
		cout << "Failed to load the intronlist, file may not exist" << endl;
		exit(-1);
	}

	// gen file name
	line = inputinst;
	int idx = line.find_last_of(".");
	stringstream stream_outfilename;
	stream_outfilename << line.substr(0,idx) << ".filtered.instance";
	string s_outfilename = stream_outfilename.str();
	ofstream s_outfile;
	s_outfile.open(s_outfilename.c_str());

	ifstream s_infile;
	s_infile.open(line.c_str());
	
	string s_line = "";
	if (s_infile.is_open())
	{
		while(s_infile.good())
		{
			getline(s_infile,s_line);
			s_outfile << s_line << endl;
			
			if (s_line[0] == 'I')
			{
				getline(s_infile,s_line);
				s_outfile << s_line << endl;
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

					
					getline(s_infile,s_line);
					s_outfile << s_line << endl;
					if (s_line[0] == 'R')
					{
						getline(s_infile,s_line);
						s_outfile << s_line << endl;
						if (s_line[0] == 'S')
						{
							int idx = s_line.find_first_of("\t");
							assert(idx > -1);
							int Nsegs = atoi(s_line.substr(idx+1).c_str());
							int count1 = 0;
							int intron_s1 = -1;
							int intron_s2 = -1;
							for (int i = 0; i < Nsegs; i++)
							{
								getline(s_infile,s_line);
								istringstream ss(s_line);
								vector<int> token_v;
								string token;
								while(getline(ss,token,'\t'))
									token_v.push_back(atoi(token.c_str()));
								assert(token_v.size()==9);

								// test if the exon is actually an intron
								if (token_v[3] > 0)
								{
									map<int,int> s_intron_list = intron_list[chr1];

									int pos1 = token_v[0];
									int pos2 = token_v[1];
									map<int,int>::iterator it = s_intron_list.lower_bound(pos1);
							
									string repline = s_line;
									int rep = 0;
									//condition 1
									if (it == s_intron_list.end())
									{
					
									}
									else
									{
										if (it->first != pos1)
										{
											if (it->first >= pos2)
											{
												
											}
											else
											{
												rep = 1;
											}
										}
										else
											rep = 1;
										
									}
									
									if (it != s_intron_list.begin())
									{
										map<int,int>::iterator it1 = it;
										it1--;
										if (it1->second <= pos1)
										{
										}
										else
											rep = 1;
									}
									
									if (rep == 1)
									{
										s_outfile << fixed << token_v[0] << "\t" << token_v[1] << "\t" << token_v[2] << "\t0\t0\t0\t0\t1\t0" << endl;
									}
									else
										s_outfile << s_line << endl;
									
								}
								else
								{
									s_outfile << s_line << endl;
								}

							}
							
						}
							
					}
			
				}
			}
		}
	}
	cout << "Finished " << line << endl;
	s_infile.close();
	s_outfile.close();
}
