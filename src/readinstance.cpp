#include "readinstance.h"
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


void readinstance::readinstance_p(string inputfile, string outputfile)
{
	//cout << "Processing\t" << inputfile << "\toutput\t" << outputfile << endl;
	ifstream infile1;
	infile1.open(inputfile.c_str());
	std::remove(outputfile.c_str());
	
	string line;
	Instance * inst = new Instance();
	Info info;
	inst->maxExonCov = -1;
	int instcount = 0;
	int outexon = 0, outexon_abun = 0, outexon_cvg = 0, outSGT = 0, outPET = 0, outSGT_abun = 0, outp0 = 0;
	vector<vector<int> > ExonBound;
	vector<vector<double> > exonAbundance;
	vector<vector<int> > SegConf;
	vector<vector<int> > SegAbundance;
	vector<vector<double> > exoncvg;
	vector<int> totalNumOfReads;
	vector<double> exonp0;
	int sumReadstemp = 0;
	string minloc, maxloc;
	int totalNumOfPaths = 0;
	int totalNumOfGenes = 0;
	int countlinetest = 0;
	if (infile1.is_open())
	{
		while (infile1.good())
		{
			getline(infile1,line);
			//cout << "line: " << ++countlinetest << endl;
			int idx; int a;
			idx = line.find_first_of("\t");
			string substring = line.substr(0,idx);
			vector<int> singleExonBound;
			vector<int> singleSegConf;
			vector<vector<int> > singlePETreads1;
			vector<int> singlePETreads;
			vector<int> singlePETconn;
			vector<double> singleexonAbundance;
			vector<int> singleSegAbundance;
			vector<double> singleexoncvg;
			
			if (!substring.compare("Instance") || line.empty())
			{
				if (instcount > 0)
				{
					//process....
					if (instcount % 1 == 0)
						//cout << "process\t" << instcount << "\t" << inst->location << endl;//"\t" << inst->NumofSegs << endl; 

					if (!inst->strand.compare("-"))
					{
						std::reverse(inst->exonbound.begin(), inst->exonbound.end());
						for (int i = 0; i < inst->exonAbundance.size(); i++)
							std::reverse(inst->exonAbundance[i].begin(), inst->exonAbundance[i].end());
						std::reverse(inst->exonAbundance.begin(), inst->exonAbundance.end());

						for (int i = 0; i < inst->exoncvg.size(); i++)
							std::reverse(inst->exoncvg[i].begin(), inst->exoncvg[i].end());
						std::reverse(inst->exoncvg.begin(), inst->exoncvg.end());

						for (int i = 0; i < inst->SegConf.size(); i++)
							std::reverse(inst->SegConf[i].begin(), inst->SegConf[i].end());
						std::reverse(inst->SegConf.begin(), inst->SegConf.end());
						std::reverse(inst->SegAbundance.begin(), inst->SegAbundance.end());
						std::reverse(inst->connSeg.begin(), inst->connSeg.end());
						std::reverse(inst->exonp0.begin(), inst->exonp0.end());
						for (int i = 0; i < inst->PETconnect.size(); i++)
						{
							int temp = 0, temp1 = 0;
							temp = inst->NumofSegs - inst->PETconnect[i][0] + 1;
							temp1 = inst->NumofSegs - inst->PETconnect[i][1] + 1;
							inst->PETconnect[i][0] = temp1;
							inst->PETconnect[i][1] = temp;
						}
					}
					//if (inst->NumofSegs < 200 && inst->strand.compare("."))
					if (inst->strand.compare("."))
					{
						//cout << "Start processing" << endl;
						if (opt.readtype == "p")
						{
							//inst->process_inst(info, outputfile);
							inst->process_inst_se(info, outputfile);
						}
						else if (opt.readtype == "s")
							inst->process_inst_se(info, outputfile);
						else
							exit(-1);
						if (info.valid)
							infolist.push_back(info);
						//process_inst(inst,outputfile, totalNumOfPaths);
						//cout << "Current Num of Genes: " << ++totalNumOfGenes << endl;
						if (totalNumOfReads.empty())
						{
							for (int i = 0; i < inst->totalreads.size();i++)
							{
								totalNumOfReads.push_back(0);
							}
							for (int i = 0; i < inst->totalreads.size();i++)
							{
								totalNumOfReads[i] += inst->totalreads[i];
							}
						}
					}
					else
					{
						//cout << "pass\t" << inst->NumofSegs << endl;
					}
				
					for (int i = 0; i < inst->exonbound.size(); i++)
					{
						vector<int>().swap(inst->exonbound[i]);
					}
					for (int i = 0; i < inst->SegConf.size(); i++)
					{
						vector<int>().swap(inst->SegConf[i]);
						inst->SegConf[i].clear();
					}
					vector<vector<int> >().swap(inst->exonbound);
					vector<vector<int> >().swap(inst->SegConf);

					for (int i = 0; i < inst->SegAbundance.size(); i++)
					{
						vector<int>().swap(inst->SegAbundance[i]);
						inst->SegAbundance[i].clear();
					}
					vector<vector<int> >().swap(inst->SegAbundance);

					inst->exonp0.clear();
					vector<double>().swap(inst->exonp0);

					vector<int>().swap(inst->connSeg);

					for (int i = 0; i < inst->exonAbundance.size(); i++)
					{
						vector<double>().swap(inst->exonAbundance[i]);
						inst->exonAbundance[i].clear();
					}

					for (int i = 0; i < inst->exoncvg.size(); i++)
					{
						vector<double>().swap(inst->exoncvg[i]);
						inst->exoncvg[i].clear();
					}
					vector<vector<double> >().swap(inst->exonAbundance);
					vector<vector<double> >().swap(inst->exoncvg);
					inst->exonbound.clear();
					inst->SegConf.clear();
					inst->SegAbundance.clear();
					inst->connSeg.clear();
					inst->exonAbundance.clear();
					inst->exoncvg.clear();
					delete inst;
					inst = NULL;
					inst = new Instance();
					Info info;
					//inst->frag_mu = frag_mu;
					//inst->frag_sigma = frag_sigma;
					outexon = 0; outexon = 1;
					for (int i = 0; i < ExonBound.size(); i++)
					{
						vector<int>().swap(ExonBound[i]);
						ExonBound[i].clear();
					}
					for (int i = 0; i < SegConf.size(); i++)
					{
						vector<int>().swap(SegConf[i]);
						SegConf[i].clear();
					}
					vector<vector<int> >().swap(ExonBound);
					ExonBound.clear();
					vector<vector<int> >().swap(SegConf);
					SegConf.clear();
					for (int i = 0; i < SegAbundance.size(); i++)
					{
						vector<int>().swap(SegAbundance[i]);
						SegAbundance[i].clear();
					}
					vector<vector<int> >().swap(SegAbundance);
					SegAbundance.clear();
					for (int i = 0; i < exonAbundance.size(); i++)
					{
						vector<double>().swap(exonAbundance[i]);
						exonAbundance[i].clear();
					}
					vector<vector<double> >().swap(exonAbundance);
					exonAbundance.clear();

					for (int i = 0; i < exoncvg.size(); i++)
					{
						vector<double>().swap(exoncvg[i]);
						exoncvg[i].clear();
					}
					vector<vector<double> >().swap(exoncvg);
					exoncvg.clear();

					exonp0.clear();
					vector<double>().swap(exonp0);
					//cout << "Finished instance" << endl;
				}
				instcount = instcount + 1;
			}
			else if (!substring.compare("segment"))
			{
				inst->NumofExon = atoi(line.substr(idx+1).c_str());
				//cout << "segment" << inst->NumofExon << endl;
				outexon = 1;
			}
			else if (outexon > 0)
			{
				if (!substring.compare("segcount"))
				{
					inst->exonbound = ExonBound;
					outexon = 0;
					outexon_abun = 1;
					//cout << "segcount" << endl;
				}
				else
				{
					int idx1 = line.find_last_of("\t");
					int idx2 = line.find_first_of("\t");
					/*
					double exon_abun = (double)(atof(line.substr(idx1+1).c_str()));
					if (exon_abun > inst->maxExonCov)
					{
						inst->maxExonCov = exon_abun;
					}
					exonAbundance.push_back(exon_abun);*/
					singleExonBound.push_back(atoi(line.substr(idx2+1,idx1-idx2).c_str()));
					singleExonBound.push_back(atoi(line.substr(idx1+1).c_str()));
					if (outexon == 1)
						minloc = line.substr(idx2+1,idx1-idx2-1);
					maxloc = line.substr(idx1+1);
					ExonBound.push_back(singleExonBound);
					//cout << singleExonBound[0] << "\t" << singleExonBound[1] << endl;
					singleExonBound.clear();
					outexon++;
				}

			}
			else if (outexon_abun == 1)
			{
				if (!substring.compare("segcvg"))
				{
					inst->exonAbundance = exonAbundance;
					inst->NumofSamples = exonAbundance[0].size();
					exonAbundance.clear();
					outexon_abun = 0;
					outexon_cvg = 1;
				}
				else
				{
					//double exon_abun  = 0;
					int idx1 = line.find_first_of("\t");
					int preidx = idx1;
					idx1 = line.find_first_of("\t",preidx+1);
					while (idx1 != -1)
					{
						singleexonAbundance.push_back(atof(line.substr(preidx+1,idx1-preidx).c_str()));
						//cout << atof(line.substr(preidx+1,idx1-preidx).c_str()) << "\t";
						//exon_abun += atof(line.substr(preidx+1,idx1-preidx).c_str());
						preidx = idx1;
						idx1 = line.find_first_of("\t",idx1+1);
					}
					exonAbundance.push_back(singleexonAbundance);
					//cout << endl;
					singleexonAbundance.clear();
					//if (exon_abun > inst->maxExonCov)
					//{
					//	inst->maxExonCov = exon_abun;
					//}
				}

			}
			else if (outexon_cvg == 1)
			{
				if (!substring.compare("segp0"))
				{
					inst->exoncvg = exoncvg;
					exoncvg.clear();
					outp0 = 1;
					outexon_cvg = 0;
				}
				else
				{
					int idx1 = line.find_first_of("\t");
					int preidx = idx1;
					idx1 = line.find_first_of("\t",preidx+1);
					while (idx1 != -1)
					{
						singleexoncvg.push_back(atof(line.substr(preidx+1,idx1-preidx).c_str()));
						preidx = idx1;
						idx1 = line.find_first_of("\t",idx1+1);
					}
					exoncvg.push_back(singleexoncvg);
					singleexoncvg.clear();
				}

			}
			else if (outp0 == 1)
			{
				if (!substring.compare("segmaxcvg"))
				{
					inst->exonp0 = exonp0;
					exonp0.clear();
					vector<double>().swap(exonp0);
					outp0 = 0;
				}
				else
				{
					int idx2 = line.find_first_of("\t");
					string confstring = line.substr(idx2+1);
					int idx1 = confstring.find_first_of("\t");
					int preidx = -1;
					int count_exon = 0;
					double max_p0 = -1;
					double count_p0 = 0;
					while (idx1 != -1)
					{
						double getp0 = atof(confstring.substr(preidx+1,idx1-preidx).c_str());
						double p0_th = 0.85;
						if (getp0 > p0_th)
							count_p0++;
						if (max_p0 < getp0)
							max_p0 = getp0;
						preidx = idx1;
						idx1 = confstring.find_first_of("\t",idx1+1);
					//	cout << getp0 << " ";
					}
					//cout << endl;
					//exonp0.push_back(max_p0);
					exonp0.push_back(count_p0);
					//cout << max_p0 << endl;
				}

			}
			else if (!substring.compare("sgtype"))
			{
				inst->NumofSegs = atoi(line.substr(idx+1).c_str());
				//cout << "sgtype\t" << inst->NumofSegs << endl;
				outSGT = 1;
			}
			else if (outSGT == 1)
			{
				if (!substring.compare("sgcount"))
				{
					inst->SegConf = SegConf;
					outSGT = 0;
					outSGT_abun = 1;
					//cout << "sgcount" << endl;
				}
				else
				{
					int idx2 = line.find_first_of("\t");
					string confstring = line.substr(idx2+1);
					int idx1 = confstring.find_first_of(" ");
					int preidx = -1;
					int count_exon = 0;
					while (idx1 != -1)
					{
						singleSegConf.push_back(atoi(confstring.substr(preidx+1,idx1-preidx).c_str()));
						//cout << atoi(confstring.substr(preidx+1,idx1-preidx).c_str()) << "\t";
						count_exon += atoi(confstring.substr(preidx+1,idx1-preidx).c_str());
						preidx = idx1;
						idx1 = confstring.find_first_of(" ",idx1+1);
					}
					//cout << endl;
					singleSegConf.pop_back();
					SegConf.push_back(singleSegConf);
					inst->connSeg.push_back(count_exon);
				}
			}
			else if (outSGT_abun == 1)
			{
				//cout << "process SegAbundance" << endl;
				if (!substring.compare("petype"))
				{
					inst->SegAbundance = SegAbundance;
					outSGT_abun = 0;
				}
				else
				{
					int idx2 = line.find_first_of("\t");
					string confstring = line.substr(idx2+1);
					//cout << confstring << endl;
					int idx1 = confstring.find_first_of("\t");
					int preidx = -1;
					while (idx1 != -1)
					{
						singleSegAbundance.push_back(atoi(confstring.substr(preidx+1,idx1-preidx).c_str()));
						//cout << atoi(confstring.substr(preidx+1,idx1-preidx).c_str()) << "\t";
						preidx = idx1;
						idx1 = confstring.find_first_of("\t",idx1+1);
					}
					//cout << endl;
					SegAbundance.push_back(singleSegAbundance);
				}
			}
			else if (!substring.compare("Chr"))
			{
				inst->location = line.substr(idx+1) + "\t" + minloc + "\t" + maxloc + "\t";
				getline(infile1,line);
				int idx1 = line.find_last_of("\t");
				inst->strand = line.substr(idx1+1);
				inst->location = inst->location + inst->strand;
				//cout << inst->location << endl;
			}
			else if (!substring.compare("readlen"))
			{
				getline(infile1,line);
				int idx1 = line.find_first_of("\t");
				int preidx = -1;
				//cout << "readlen" << endl;
				while (idx1 != -1)
				{
					inst->ReadLen.push_back(atoi(line.substr(preidx+1,idx1-preidx).c_str()));
					//cout << atoi(line.substr(preidx+1,idx1-preidx).c_str()) << "\t";
					preidx = idx1;
					idx1 = line.find_first_of("\t",idx1+1);
				}
				//cout << endl;
			}
			else if (!substring.compare("totalreads"))
			{
				getline(infile1,line);
				int idx1 = line.find_first_of("\t");
				int preidx = -1;
				//cout << "totalreads" << endl;
				while (idx1 != -1)
				{
					inst->totalreads.push_back((double)atoi(line.substr(preidx+1,idx1-preidx).c_str()));
					//cout << atoi(line.substr(preidx+1,idx1-preidx).c_str()) << endl;
					//cout << atoi(line.substr(preidx+1,idx1-preidx).c_str()) << "\t";
					preidx = idx1;
					idx1 = line.find_first_of("\t",idx1+1);
				}
				TOTALNUMREADS = inst->totalreads;
				//cout << endl;
			}

		}
		
					
	}
}
