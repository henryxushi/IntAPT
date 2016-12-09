#include "utility.h"
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
using namespace std;

bool ifcontain_vec(vector<int> &a, vector<int> &b);
bool ifcontain_arr(int a[], vector<int> &b);
double normpdf(double x, double mu, double sigma); // not exactly the probability, it is relative one
double phi(double x);
double phi_range(double x, double y);
double effectL(int l1, int l2, int l3, int r, double mu, double sigma, int bias);
PathVar::PathVar()
{
	minjuncfrac=0.01;
	minparentjuncfrac = 0.01;
	minadjjuncfrac=2;
	minintronretentioncvgfrac=0.4;
	minintronretentionlen=10;
}

void PathVar::setsinglepar()
{
	minjuncfrac=0.2;
	minparentjuncfrac = 0.2;
	minadjjuncfrac=2;
	minintronretentioncvgfrac=0.4;
	minintronretentionlen=10;
}

void PathVar::setdefaultpar()
{
	minjuncfrac=0.01;
	minparentjuncfrac = 0.01;
	minadjjuncfrac=2;
	minintronretentioncvgfrac=0.4;
	minintronretentionlen=10;
}

void Instance::process_inst_se(Info& info, string outputfile)
{
	//enumberate path
	info.clear();
	if (SegAbundance.size() == 0)
		return;
	
	
	if (true) 
	{
		L.clear();
		for (int i = 0; i < NumofExon; i++)
		{
			L.push_back((exonbound[i][1]) - (exonbound[i][0]) + 1);
		}
		int output = 1;
		ofstream outfile1;
		outfile1.open(outputfile.c_str(),std::ios_base::app);
		remove_low();
		L.clear();
		for (int i = 0; i < NumofExon; i++)
		{
			L.push_back((exonbound[i][1]) - (exonbound[i][0]) + 1);
		}
		cal_junccount();
		vector<vector<double> > proc_abun;
		vector<vector<double> > proc_R;
		vector<vector<double> > proc_L;
		map<int, vector<int> > iscontained;
		findcontain(iscontained);
		outfile1 << "Instance\t" << location << endl;
		info.label = location;
		outfile1 << "ExonBound\t" << exonbound.size() << endl;
		for (int i = 0; i < exonbound.size(); i++)
		{
			outfile1 << exonbound[i][0] << "\t" << exonbound[i][1] << endl;
		}
		info.modexonbound(exonbound);
		
		double th1 = 0;
		double min_th1 = 0;

		mod_y_se(iscontained, proc_abun, proc_R, proc_L, outfile1, output, th1, min_th1);

		map<int,set<int> > exonmap; 
		getmap(exonmap, th1); // exonmap minvalue 1
		/*for (map<int,set<int> >::iterator mapit = exonmap.begin();mapit!=exonmap.end();++mapit)
		{
			cout << mapit->first << ": ";
			for (set<int>::iterator setit = mapit->second.begin(); setit != mapit->second.end();setit++)
				cout << *setit <<" ";
			cout << endl;
		}*/

		//if(!exonmap.empty())
		//{
			set<int>::iterator it;
			int not_source_list[NumofExon];
			int not_sink_list[NumofExon];
			for (int i = 0; i < NumofExon; i++)
			{
				not_source_list[i] = 0;
				not_sink_list[i] = 0;
			}

			for (map<int,set<int> >::iterator ii = exonmap.begin(); ii != exonmap.end(); ++ii)
			{
				int key = (*ii).first;
				set<int> temp = (*ii).second;
				if (!temp.empty())
				{
					not_sink_list[key - 1] = 1;
					for (it = temp.begin(); it != temp.end(); ++it)
					{
						not_source_list[*it - 1] = 1; // not a source node
					}
				}
			}
			
			vector<vector<int> > paths;
			double mod_y_th = 0;
			int quit = enumerate_path(exonmap, not_source_list, not_sink_list, outfile1, paths);
			while(quit == -1)
			{
				for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
				{
					vector<int>().swap(it->second);
					//exonmap.erase(it);
				}

				map<int,vector<int> >().swap(iscontained);
				iscontained.clear();
				findcontain(iscontained);
				th1 = min_th1;
				if (th1 == 10000)
				{
					break;
				}
				//cout << "new threshold\t" << th1 << endl;
				output = 0;
				proc_abun.clear();
				vector<vector<double> >().swap(proc_abun);
				proc_R.clear();
				vector<vector<double> >().swap(proc_R);
				proc_L.clear();
				vector<vector<double> >().swap(proc_L);
				mod_y_se(iscontained, proc_abun, proc_R, proc_L, outfile1, output, th1, min_th1);
				//cout << "current number of segs in inst\t" << SegConf.size() << "\t" << NumofSegs << "\t" << SegAbundance.size() << endl;
				for (map<int,set<int> >::iterator it = exonmap.begin();it!=exonmap.end();++it)
				{
					set<int>().swap(it->second);
					//exonmap.erase(it);
				}

				map<int,set<int> >().swap(exonmap);
				exonmap.clear();
				map<int,set<int> > exonmap; 
				getmap(exonmap, th1); // exonmap minvalue 1
				set<int>::iterator it;
				int not_source_list1[NumofExon];
				int not_sink_list1[NumofExon];
				for (int i = 0; i < NumofExon; i++)
				{
					not_source_list1[i] = 0;
					not_sink_list1[i] = 0;
				}
				for (map<int,set<int> >::iterator ii = exonmap.begin(); ii != exonmap.end(); ++ii)
				{
					int key = (*ii).first;
					set<int> temp = (*ii).second;
					if (!temp.empty())
					{
						//cout << inst->NumofExon << endl;
						//cout << key << endl;
						not_sink_list1[key - 1] = 1;
						for (it = temp.begin(); it != temp.end(); ++it)
						{
							not_source_list1[*it - 1] = 1; // not a source node
						}
					}
				}
				quit = enumerate_path(exonmap, not_source_list1, not_sink_list1, outfile1, paths);
			}
			proc_abun.clear();
			vector<vector<double> >().swap(proc_abun);
			proc_R.clear();
			vector<vector<double> >().swap(proc_R);
			proc_L.clear();
			vector<vector<double> >().swap(proc_L);
			output = 1;
			map<int,vector<int> >().swap(iscontained);
			iscontained.clear();
			findcontain(iscontained);
			mod_y_se(iscontained, proc_abun, proc_R, proc_L, outfile1, output, th1, min_th1);
			
			//cout << "final paths\t" << paths.size() << endl;
			for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
			{
				vector<int>().swap(it->second);
			}
			map<int,vector<int> >().swap(iscontained);
			iscontained.clear();
			findcontain(iscontained);
			
			vector<vector<int> > paths1;
			for (int i = 0; i < paths.size(); i++)
			{
				if (paths[i].size() > 2)
				{
					vector<int> singlepath1;
					int idxpath = 1;
					for (int j = 0; j < NumofExon; j++)
					{
						if (j+1 == paths[i][idxpath])
						{
							singlepath1.push_back(1);
							idxpath++;
						}
						else
							singlepath1.push_back(0);
					}
					paths1.push_back(singlepath1);
				}

			}
			

			if (paths.size() == 0)
			{
				info.valid = false;
				outfile1 << "Abundance" << endl << -1 << "\t" << -1 <<"\t" << -1 << endl;
				outfile1 << "Paths" << endl;
				outfile1 << -1 << endl << -1 << endl;
			}
			else
			{
				info.valid = true;
				info.mody(proc_abun);
				//cout << "R processed" << endl;
				/*for (int ii = 0; ii < proc_R.size(); ii++)
				{
					for (int jj = 0; jj < proc_R[ii].size(); jj++)
						cout << proc_R[ii][jj] << " ";
					cout << endl;
				}*/
				info.modR(proc_R);
				vector<double> push_L;
				for (int ii = 0; ii < NumofSegs; ii++)
					push_L.push_back(proc_L[0][ii]);
				info.modL(push_L);
				outfile1 << "Abundance" << endl;
				for (int ii = 0; ii < NumofSegs; ii++)
				{
					outfile1 << proc_L[0][ii] << "\t";
					for (int i = 0; i < proc_abun.size(); i++)
						outfile1 << proc_R[i][ii] << "\t";
					outfile1 << endl;
				}
				

				outfile1 << "Paths" << endl;
				vector<vector<int> > segpath_all;
				vector<vector<int> > paths_exon;
				for (int i = 0; i < paths.size(); i++)
				{
					vector<int> segpath;
					vector<int> paths_vec;
					int countpath1 = 0;
					for (int ii = 0; ii < paths[i].size(); ii++)
						countpath1++;
					if (countpath1 > 0) //may skip single exon transcripts
					{
						mod_X(paths[i],segpath,paths_vec,outfile1);
					}
					if (segpath.size() > 0)
					{
						segpath_all.push_back(segpath);
						paths_exon.push_back(paths_vec);
					}
					
				}
				//cout << "modinfo " << endl;
				if (segpath_all.size() == 0)
					info.valid = false;
				else
				{
					info.valid = true;
					info.modX(segpath_all);
					info.modX_exon(paths_exon);
				}
			}
			

			for (int i = 0; i < paths.size(); i++)
			{
				paths[i].clear();
				vector<int>().swap(paths[i]);
			}
			paths.clear();
			vector<vector<int> >().swap(paths);
		//}
		for (map<int,set<int> >::iterator it = exonmap.begin();it!=exonmap.end();++it)
		{
			set<int>().swap(it->second);
			//exonmap.erase(it);
		}

		map<int,set<int> >().swap(exonmap);
		exonmap.clear();

		for (map<int,vector<int> >::iterator it = iscontained.begin();it!=iscontained.end();++it)
		{
			vector<int>().swap(it->second);
		}
		//cout << "Finish moding and assign single" << endl;
		map<int,vector<int> >().swap(iscontained);
		iscontained.clear();
		if (info.valid)
		{
			info.single = true;
			if (info.X.rows() > 0)
			{
				if (info.X.cols() > 1)
					info.single = false;
			}
		}
	}
}


void Instance::remove_low()
{
	cal_sumexoncvg();
	double th = 0.01;
	double abun_th = 15;
	double abun_th_low = 0;//2;
	vector<int> removed_exon;
	vector<double> mean_abun;
	double max_mean_abun = 0;
	int Nsample_sup[NumofExon];
	//int count = 0;

	for (int i = 0; i < NumofExon; i++)
	{
		double mean_abun1 = 0;
		Nsample_sup[i] = 0;
		for (int ii = 0; ii < NumofSamples; ii++)
		{
			//mean_abun1 += (exonAbundance[i][ii] / L[i]);
			mean_abun1 += exonAbundance[i][ii];
			if (exonAbundance[i][ii]/L[i] > 3)
				Nsample_sup[i]+=100;
			else if (exonAbundance[i][ii]/L[i] > 0)
				Nsample_sup[i]+=1;
		}
		mean_abun1 /= (double)NumofSamples;
		mean_abun.push_back(mean_abun1);
		if (max_mean_abun < mean_abun1)
			max_mean_abun = mean_abun1;
	}

	for (int i = 0; i < NumofExon; i++)
	{
		th = 0; //do not trust the th due to overdispersion.
		//bool ifremove = (mean_abun[i] < max_mean_abun * th) || (mean_abun[i] < -1 || Nsample_sup[i] < 0.2 * NumofSamples || exonp0[i] > NumofSamples - 1 || sumexoncvg[i] / NumofSamples < 3); //use coverage to cut not work for ENCSR310FIS, over-dispersion obvious 
		bool ifremove = (mean_abun[i] < max_mean_abun * th) || (mean_abun[i] < -1 ||(exonp0[i] > NumofSamples - 1) || sumexoncvg[i] / NumofSamples < 0 || Nsample_sup[i] < 0.2 * NumofSamples); //use coverage to cut not work for ENCSR310FIS, over-dispersion obvious 
		//cout << i << ": " << mean_abun[i] <<  " " << Nsample_sup[i] << " " << exonp0[i] << " " << sumexoncvg[i] << endl;
		if (ifremove)
			removed_exon.push_back(i+1);
	}

	mean_abun.clear();

	set<int> removed_seg;
	for (int i = 0; i < NumofSegs; i++)
	{
		for (int j = 0; j < removed_exon.size(); j++)
		{
			if(SegConf[i][removed_exon[j]-1] == 1)
			{
				removed_seg.insert(i+1);
				break;
			}
		}	
		int Nsample_sup_seg = 0;
		for (int j = 0; j < NumofSamples; j++)
		{
			if (SegAbundance[i][j] > 0)
				Nsample_sup_seg++;
		}
		if (Nsample_sup_seg < 0.2 * NumofSamples)// <= 1) NOTICE HERE!!!!!!!!!!!!!!
			removed_seg.insert(i+1);
	}

	for (int i = 0; i < removed_exon.size(); i++)
	{
		//cout << "removed\t" << removed_exon[i]-1 << endl;
		vector<int>().swap(exonbound[removed_exon[i]-1-i]);
		exonbound.erase(exonbound.begin()+removed_exon[i]-1-i);
		vector<double>().swap(exonAbundance[removed_exon[i]-1-i]);
		exonAbundance.erase(exonAbundance.begin()+removed_exon[i] - 1 -i);
		vector<double>().swap(exoncvg[removed_exon[i]-1-i]);
		exoncvg.erase(exoncvg.begin()+removed_exon[i] - 1 -i);

		for (int j = 0; j < NumofSegs; j++)
		{
			//cout << "i " << i << " j " << j << " SegConf[j].size() " << SegConf[j].size() << " removed_exon[i] " << removed_exon[i] << endl; 
			SegConf[j].erase(SegConf[j].begin()+removed_exon[i] - 1 - i);
		}
	}
	NumofExon = NumofExon - removed_exon.size();

	

	// remove duplicate segments
	/*for (int i = 0; i < SegConf.size(); i++)
	{
		for (int j = i + 1; j < SegConf.size(); j++)
		{
			int countsame = 0;
			for (int ii = 0; ii < SegConf[i].size(); ii++)
			{
				countsame += abs(SegConf[i][ii] - SegConf[j][ii]);
			}
			if (countsame == 0)
			{
				for(int jj = 0; jj < NumofSamples; jj++)
				{
					SegAbundance[i][jj] += SegAbundance[j][jj];
				}
				removed_seg.insert(removed_seg.end(),j+1);
			}
		}
	}*/
	
	// remove all 0's segments
	for (int j = 0; j < SegConf.size(); j++)
	{
		int count = 0;
		for (int ii = 0; ii < SegConf[j].size(); ii++)
		{
			if (SegConf[j][ii] == 1)
				count++;
		}
		if (count == 0)
		{
			removed_seg.insert(removed_seg.end(),j+1);
		}

	}
	
	set<int>::iterator it;
	int removetemp = 0;
	for(it = removed_seg.begin(); it != removed_seg.end(); ++it)
	{	
		//cout << "remove seg (remove_low)" << *it << endl; 
		vector<int>().swap(SegConf[*it-1-removetemp]);
		SegConf.erase(SegConf.begin() + *it - 1 - removetemp);
		vector<int>().swap(SegAbundance[*it-1-removetemp]);
		SegAbundance.erase(SegAbundance.begin() + *it - 1 - removetemp);
		removetemp++;
	}

	set<int>().swap(removed_seg);


	vector<int> connSeg1 (SegConf.size(),0);
	for (int j = 0; j < SegConf.size(); j++)
	{
		//cout << j << "\t" << SegConf.size() << endl;
		
		for (int ii = 0; ii < SegConf[j].size(); ii++)
		{
			if (SegConf[j][ii] == 1)
			{
				connSeg1[j]++;
			}
		}
	}
	
	//cout << connSeg.size() << endl;
	connSeg.swap(connSeg1);
	//cout << connSeg.size() << endl;
	vector<int>().swap(connSeg1);
	connSeg1.clear();
	NumofSegs = connSeg.size();
}

void Instance::cal_sumexoncvg()
{
	sumexoncvg.clear();
	sumexoncvg.resize(NumofExon,0);
	for (int i = 0; i < NumofExon; i++)
	{
		for (int j = 0; j < NumofSamples; j++)
			sumexoncvg[i] += exoncvg[i][j];
	}
}

void Instance::cal_junccount()
{
	junccount.clear();
	junccount.resize(NumofExon,vector<double> (NumofExon,0));
	juncsup.clear();
	juncsup.resize(NumofExon,vector<double> (NumofExon,0));
	for (int i = 0; i < NumofExon; i++)
	{
		for (int j = i+1; j < NumofExon; j++)
		{
			junccount[i][j] = cal_connreads(i,j);
			juncsup[i][j] = cal_connsup(i,j);
		}
	}
	sumexoncvg.clear();
	sumexoncvg.resize(NumofExon,0);
	for (int i = 0; i < NumofExon; i++)
	{
		for (int j = 0; j < NumofSamples; j++)
			sumexoncvg[i] += exoncvg[i][j];
	}

}

double Instance::cal_connreads(int i, int j)
{
	double count = 0;
	for(int n = 0; n < SegConf.size(); n++)
	{
		if (SegConf[n][i] == 1 and SegConf[n][j] == 1)
		{
			int flag = 0;
			for (int m = i+1; m <= j-1; m++)
				flag += SegConf[n][m];
			if (flag == 0)
			{
				for (int m = 0; m < NumofSamples; m++)
					count += SegAbundance[n][m];
			}
		}
	}
	return count;
}

double Instance::cal_connsup(int i, int j)
{
	set<int> count;
	for(int n = 0; n < SegConf.size(); n++)
	{
		if (SegConf[n][i] == 1 and SegConf[n][j] == 1)
		{
			int flag = 0;
			for (int m = i+1; m <= j-1; m++)
				flag += SegConf[n][m];
			if (flag == 0)
			{
				for (int m = 0; m < NumofSamples; m++)
				{
					if (SegAbundance[n][m] > 0)
						count.insert(m);
				}
			}
		}
	}
	return (double)count.size();
}

void Instance::findchildren(int c, vector<int> &nodes)
{
	nodes.clear();
	double maxj = 0; //maximum junction count of c
	double maxnonjj = 0; //maximum non-adjancent count of c
	int ncandidate = 0;
	// determine the value of maxj and maxnonjj
	for (int i = c+1; i < NumofExon; i++)
	{
		if (junccount[c][i] > maxj)
			maxj = junccount[c][i];
		if (i > c+1 and junccount[c][i] > maxnonjj)
			maxnonjj = junccount[c][i];
		if (junccount[c][i] > 0)
			ncandidate++;
	}
	//cout << c+1 << ": " << ncandidate << " candidates" << endl;
	bool jumpretention = false;
	for (int i = c+1; i < NumofExon; i++)
	{
		if (junccount[c][i] > 0)
		{
			pathvar.setdefaultpar();
			if (juncsup[c][i] < min((double)2,0.2*NumofSamples))
				pathvar.setsinglepar();
			//if (i == c+1) //here should consider genomic location
			bool isaccepted = false;
			if (abs(exonbound[i][0]-exonbound[c][1]) == 1 or abs(exonbound[i][1]-exonbound[c][0]) == 1)
			{
				//intron retention
				//cout << "retention " << c+1 << " " << i+1 << endl;
				double mfrac = pathvar.minjuncfrac * pathvar.minadjjuncfrac;
				if (mfrac>1)
					mfrac = 1;
				double cvgc = sumexoncvg[c]; //total coverage of c
				double cvgi = sumexoncvg[i]; //total coverage of i
				//cout << "coverage " << cvgc << " " << cvgi  << " " << cvgc*pathvar.minintronretentioncvgfrac << endl;
				int ccount = junccount[c][i];
				//bool isaccepted = false;
				if (ncandidate == 1) isaccepted = true;
				if (ccount >= maxj*mfrac)
				{
					if (L[i] < pathvar.minintronretentionlen) //see if the length is short
						isaccepted = true;
					else
					{
						if(cvgi >= cvgc*pathvar.minintronretentioncvgfrac) //and cvgi>2 //control the coverage difference
							isaccepted = true;
					}
				}
				if (isaccepted)
					nodes.push_back(i);
				else
					jumpretention = true;
			}
			//else
			if (isaccepted == false)
			{
				//junctions
				double tmaxj = maxj;
				if (jumpretention) tmaxj = maxnonjj;
				if (junccount[c][i] >= tmaxj*pathvar.minjuncfrac)
				{
					nodes.push_back(i);
					//cout << "push " << i << endl;
				}
				else
				{
					if (junccount[c][i] >= maxnonjj*pathvar.minjuncfrac)
					{
						nodes.push_back(i);
					//	cout << "push " << i << endl;
					}
					//else
					//	cout << "push " << i << " failed" << endl;
				}
			}
		}
		pathvar.setdefaultpar();
	}

}

void Instance::getmap(map<int,set<int> > &exonmap, double &th)
{
	cal_junccount();
	cal_sumexoncvg();
	//cout << "c graph:" << endl;
	map<int,set<int> > rev_exonmap;
	for (int i = 0; i < NumofExon; i++)
	{
		vector<int> nodes;
		findchildren(i,nodes);
		if (nodes.size()>0)
		{
			//cout << i+1 << ": ";
			for (int j = 0; j < nodes.size(); j++)
			{
				//cout << nodes[j]+1 << " ";
				exonmap[i+1].insert(nodes[j]+1);
				rev_exonmap[nodes[j]+1].insert(i+1);
			}
			//cout << endl;
		}
	}
	/*cout << "map before p:" << endl;
	for (map<int,set<int> >::iterator mapit = exonmap.begin();mapit!=exonmap.end();++mapit)
	{
		cout << mapit->first << ": ";
		for (set<int>::iterator setit = mapit->second.begin(); setit != mapit->second.end();setit++)
			cout << *setit <<" ";
		cout << endl;
	}*/
	refineparents(exonmap,rev_exonmap);
	/*cout << "map after p:" << endl;
	for (map<int,set<int> >::iterator mapit = exonmap.begin();mapit!=exonmap.end();++mapit)
	{
		cout << mapit->first << ": ";
		for (set<int>::iterator setit = mapit->second.begin(); setit != mapit->second.end();setit++)
			cout << *setit <<" ";
		cout << endl;
	}*/
}

void Instance::refineparents(map<int,set<int> > &exonmap, map<int,set<int> > &rev_exonmap)
{
	//cout << "remove parents:" << endl;
	for(map<int,set<int> >::iterator it = rev_exonmap.begin(); it != rev_exonmap.end(); it++)
	{
		int node_id = it->first;
		if (it->second.size()>1)
		{
			vector<int> remove_nodes;
			findparent(node_id-1,it->second, remove_nodes);
			//cout << "node *" << node_id << ": ";
			for(int i = 0; i < remove_nodes.size(); i++)
			{
				//cout << remove_nodes[i] << "-";
				map<int,set<int> >::iterator rmit = exonmap.find(remove_nodes[i]);
				set<int>::iterator rmsetit = rmit->second.find(node_id);
				//cout << *rmsetit << "-" << rmit->second.size() <<"-";
				rmit->second.erase(rmsetit);
				//cout << rmit->second.size() <<" ";
				if (rmit->second.size() == 0)
					exonmap.erase(rmit);
			}
			//cout << endl;
		}
		

	}
}

void Instance::findparent(int c, set<int> &parents, vector<int> &remove_nodes)
{
	// c must be index which is node -1
	remove_nodes.clear();
	double maxj = 0; //maximum junction count of c
	double maxnonjj = 0; //maximum non-adjancent count of c
	int ncandidate = 0;
	// determine the value of maxj and maxnonjj
	for (set<int>::iterator it = parents.begin(); it != parents.end(); it++)
	{
		if (junccount[*it-1][c] > maxj)
			maxj = junccount[*it-1][c];
	}
	maxnonjj = maxj;
	ncandidate = parents.size();

	for (set<int>::iterator it = parents.begin(); it != parents.end(); it++)
	{
		pathvar.setdefaultpar();
		if (juncsup[*it-1][c] < min((double)2,0.2*NumofSamples))
			pathvar.setsinglepar();

		if (junccount[*it-1][c] < maxj * pathvar.minparentjuncfrac)
			remove_nodes.push_back(*it);
	}
}

int Instance::enumerate_path(map<int,set<int> > &exonmap, int not_source_list[], int not_sink_list[], ofstream &outfile, vector<vector<int> > &paths)
{
	for (int i = 0; i < NumofExon; i++)
	{
		if (not_source_list[i] == 0)
		{
			//cout << "source " << i + 1 << endl;
			set<int>::iterator it = exonmap[0].end();
			exonmap[0].insert(it,i+1);
		}
	}
	for (int i = 0; i < NumofExon; i++)
	{
		//cout << not_sink_list[i] << endl;
		if (not_sink_list[i] == 0)
		{
			//cout << "sink " << i + 1 << endl;
			set<int>::iterator it = exonmap[i+1].end();
			exonmap[i+1].insert(it,NumofExon+1);
		}

	}

	//vector<vector<int> > paths;
	//vector<int> single_path;
	//apply depth first search
	list<Seg> stack;
	Seg* seg = new Seg();
	seg->node = 0;
	seg->visited.push_back(0);
	int goal = NumofExon + 1;
	stack.push_back(*seg);
	vector<int>().swap(seg->visited);
	delete seg;
	seg = NULL;
	int countpath = 0;
	while(!stack.empty())
	{
		Seg seg1 = stack.back();
		stack.pop_back();
		if (seg1.node == goal)
		{
			if (seg1.visited.size() > 3) // path with seg larger than 1
			{
				int path_L = 0;
				for (int i = 1; i < seg1.visited.size()-1; i++)
				{
					path_L += L[seg1.visited[i]-1];
					//cout << seg1.visited[i] << " ";
				}
				//cout << "L:" << path_L << endl;
				if (path_L > 100)
				{
					paths.push_back(seg1.visited);
					//mod_X(seg1.visited, SegConf, NumofSegs, NumofExon, outfile);
					countpath++;
					//path_all.push_back(seg1.visited);
				}
				
			}
		}
		else
		{
			set<int> conn = exonmap[seg1.node];
			for (set<int>::iterator it = conn.begin(); it != conn.end(); ++it)
			{
				Seg* push_seg = new Seg();
				push_seg->node = *it;
				push_seg->visited = seg1.visited;
				push_seg->visited.push_back(*it);
				stack.push_back(*push_seg);
				vector<int>().swap(push_seg->visited);
				delete push_seg;
				push_seg = NULL;

			}
		}
		if (countpath > MAX_NUM_PATH)
		{
			for (int i = 0; i < paths.size(); i++)
			{
				paths[i].clear();
				vector<int>().swap(paths[i]);
			}
			paths.clear();
			vector<vector<int> >().swap(paths);
			return -1;
		}
			
	}
	if (countpath <= MAX_NUM_PATH)
	{
		/*for (int i = 0; i < paths.size(); i++)
			mod_X(paths[i],SegConf,NumofSegs,NumofExon,outfile);
		for (int i = 0; i < paths.size(); i++)
		{
			paths[i].clear();
			vector<int>().swap(paths[i]);
		}
		paths.clear();
		vector<vector<int> >().swap(paths);	*/
		return 0;
	}
	else
	{
		for (int i = 0; i < paths.size(); i++)
		{
			paths[i].clear();
			vector<int>().swap(paths[i]);
		}
		paths.clear();
		vector<vector<int> >().swap(paths);
		return -1;
	}

	
} 

void Instance::mod_X(vector<int> &path, vector<int> &segpath, vector<int> &paths_vec, ofstream &outfile)
{
	//transpose of final X !!!!!!!!!!!
	//cout << NumofExon << endl;
	paths_vec.clear();
	paths_vec.resize(NumofExon,0);
	int path_vec[NumofExon];
	memset(path_vec,0,sizeof(path_vec));
	for (int i = 1; i < path.size() - 1; i++)
	{
		paths_vec[path[i]-1] = 1;
		path_vec[path[i]-1] = 1;
		//cout << path[i] << " ";
	}
	//cout << endl;

	int seg_map[NumofSegs];
	memset(seg_map,0,sizeof(seg_map));
	int countmap = 0; double abun = 0;
	for (int i = 0; i < NumofSegs; i++)
	{
		if(ifcontain_arr(path_vec,SegConf[i]))
		{
			seg_map[i] = 1;
			countmap++;
			//abun += proc_abun[i];
		}
		else
		{
			seg_map[i] = 0;
		}
	}
	//double path_mean_cov = abun / (double) countmap;
	//if (countmap > 0 && abun > 0.1)
	if (countmap > 0)
	{
		for (int i = 0; i < NumofExon; i++)
		{
			outfile << path_vec[i] << " ";
		}
		outfile << endl;
		for (int i = 0; i < NumofSegs; i++)
		{
			outfile << seg_map[i] << " ";
		}
		outfile << endl;
		for (int i = 0; i < NumofSegs; i++)
		{
			segpath.push_back(seg_map[i]);
		}
	}
	else
	{
		segpath.clear();
	}
	
}

void Instance::mod_y_se(map<int, vector<int> > &iscontained, vector<vector<double> > &SegAbundance1, vector<vector<double> > &proc_R, vector<vector<double> > &proc_L, ofstream &outfile, int &output, double &th1, double &min_abun)
{
	//cout << "here" << endl;
	//vector<vector<int> > &SegConf, vector<int> &SegAbundance, vector<vector<int> > &ExonBound, vector<int> &connSeg, vector<double> &SegAbundance1, int &NumofSegs, int NumofExon, int READLEN
	vector<int> READLEN = ReadLen;
	//vector<double> SegAbundance1;
	double max_abun = -1.0;
	double th = 0.01;
	int counttemp = 0;
	vector<vector<double> > sumy_v;
	vector<double> adjusted_length_v;
	vector<double> s_SegAbundance1;
	vector<double> s_R;
	vector<double> s_L;
	for (int n = 0; n < READLEN.size(); n++)
	{
		vector<int> SegAbundance11;
		//cout << "sample " << n << endl;
		for (int nn = 0; nn < SegAbundance.size(); nn++)
		{
			SegAbundance11.push_back(SegAbundance[nn][n]);
			//cout << SegAbundance[nn][n] << " ";
		}
		//cout << endl;
		int READLEN1 = READLEN[n];
		for (int i = 0; i < NumofSegs; i++)
		{
			int sumy = SegAbundance11[i];
			//cout << i << "\t" << sumy << "!!!" << endl;
			double abun = 0;
			int adjusted_length = 1;
			map<int,vector<int> >::iterator it;
			it = iscontained.find(i+1);
			if (it != iscontained.end())
			{
				for (int ii = 0; ii < iscontained[i+1].size(); ii++)
				{
					sumy += SegAbundance11[(iscontained[i+1])[ii]-1];
				}
				if (connSeg[i] == 1)
				{
					for(int j = 0; j < NumofExon; j++)
					{
						if((SegConf[i])[j] == 1)
						{
							adjusted_length = L[j];
							break;
						}
					}
					
				}
				else
				{
					for(int j = 0; j < NumofExon; j++)
					{
						if((SegConf[i])[j] == 1)
						{
							adjusted_length = L[j];
							break;
						}
					}
					adjusted_length = min(adjusted_length, READLEN1);
					
				}
			}
			else
			{
				if (connSeg[i] == 1)
				{
					for(int j = 0; j < NumofExon; j++)
					{
						if((SegConf[i])[j] == 1)
						{
							adjusted_length = L[j];
							break;
						}
					}
				}
				else
				{
					int sumlength = 0; int count = 0; int init_length = 0; int final_length = 0;
					for(int j = 0; j < NumofExon; j++)
					{
						if((SegConf[i])[j] == 1)
						{
							final_length = L[j];
							if (count > 0)
								sumlength += L[j];
							else
								init_length = L[j];
							count++;
						}
					}
					adjusted_length = READLEN1 - (sumlength - final_length);
					adjusted_length = min(init_length, adjusted_length);
					adjusted_length = min(final_length, adjusted_length);
				}
			}
			if (adjusted_length <= 0)
				adjusted_length = 1;
			/*cout << "SegConf " << i << ": " << endl;
			for (int kk = 0; kk < SegConf[i].size(); kk++)
			{
				if (SegConf[i][kk] == 1)
				{
					cout << kk << ":" << L[kk] << " ";
				}
			}
			cout << endl;*/
			
			//cout << "seg " << i+1 << " " << sumy << " " << adjusted_length << endl;
			abun = (double)sumy / ((double)adjusted_length);
			s_R.push_back(sumy);
			s_L.push_back(adjusted_length);
			s_SegAbundance1.push_back(abun);
			adjusted_length_v.push_back((double)adjusted_length);
			SegAbundance11.clear();
			
		}
		SegAbundance1.push_back(s_SegAbundance1);
		s_SegAbundance1.clear();
		proc_R.push_back(s_R);
		s_R.clear();
		proc_L.push_back(s_L);
		s_L.clear();
	}

	if (output == 0)
	{
		vector<int> removed_seg;
		vector<double> mean_abun;
		double max_mean_abun = 0;
		int Nsample_sup[NumofSegs];
		for (int i = 0; i < NumofSegs; i++)
		{
			double mean_abun1 = 0;
			Nsample_sup[i] = 0;
			for (int ii = 0; ii < READLEN.size(); ii++)
			{
				mean_abun1 += SegAbundance1[ii][i];
				//mean_abun1 += SegAbundance[ii][i];
				if (SegAbundance1[ii][i] > 2)
					Nsample_sup[i]+=100;
				else if (SegAbundance1[ii][i] > 0)
					Nsample_sup[i]+=1;

			}
			mean_abun1 /= (double)READLEN.size();
			mean_abun.push_back(mean_abun1);

			if (max_mean_abun < mean_abun1)
				max_mean_abun = mean_abun1;

		}

		min_abun = 10000; //due to overdispersion, th may not be a good choice; using the number of samples supported may be good
		for (int i = 0; i < NumofSegs; i++)
		{
			if (connSeg[i] > 1 && mean_abun[i] < th1+0.001 || Nsample_sup[i] < 0.2*READLEN.size())
			{
				removed_seg.push_back(i+1);
			}
			else
			{
				if (connSeg[i] > 1)
				{
					if (mean_abun[i] < min_abun)
						min_abun = mean_abun[i];
				}
			}
		}

		for (int j = 0; j < removed_seg.size(); j++)
		{
			//cout << "removed seg\t" << removed_seg[j] << endl;
			vector<int>().swap(SegConf[removed_seg[j]-1-j]);
			SegConf.erase(SegConf.begin() + removed_seg[j] - 1 - j);
			connSeg.erase(connSeg.begin() + removed_seg[j] - 1 - j);
			vector<int>().swap(SegAbundance[removed_seg[j]-1-j]);
			SegAbundance.erase(SegAbundance.begin() + removed_seg[j] - 1 - j);

		}
		NumofSegs = NumofSegs - removed_seg.size();

	}

}

void Instance::findcontain(map<int, vector<int> > &iscontained) // contain: index start from 1
{
	for (int i = 0; i < SegConf.size(); i++)
	{
		for (int j = 0; j < SegConf.size(); j++)
		{
			if (ifcontain_vec(SegConf[i],SegConf[j]) && i != j)
			{
				iscontained[j+1].push_back(i+1);
				//cout << j+1 << "\tcontain\t" << i+1 << endl;
			}
		}
	}
}

bool ifcontain_vec(vector<int> &a, vector<int> &b)
{
	int tempmax = INT_MAX;
	int ifcontain = 0; int idx1 = tempmax; int countj = 0;
	for (int ii = 0; ii < a.size(); ii++)
	{
		if (a[ii] - b[ii] == -1)
		{
			ifcontain = -1;
			return false;
		}
		else if (a[ii] - b[ii] == 1)
		{
			if (idx1 > ii)
			{
				idx1 = ii;
			}
		} 
		if (b[ii] == 1 && countj >= 0)
		{
			//cout << "here\t" << idx1 << "\t" << ii << endl;
			if (idx1 < ii)
			{
				ifcontain = -1;
				return false;
			}
			countj += 1;
		}
		
		ifcontain += a[ii] - b[ii];
	}
	if (ifcontain > 0)
	{
		return true;
	}
}

bool ifcontain_arr(int a[], vector<int> &b) //path to segs
{
	int tempmin = INT_MAX;
	int tempmax = -1;
	//int ifcontain = 0; int idx1 = tempmax; int countj = 0;
	for (int ii = 0; ii < b.size(); ii++)
	{
		if (a[ii] - b[ii] == -1)
			return false;
		if (b[ii] == 1)
		{
			//countj++;
			if (ii < tempmin)
				tempmin = ii;
			if (ii > tempmax)
				tempmax = ii;
		}
	}
	if (tempmin == tempmax)
		return true;
	else
	{
		for (int i = tempmin+1; i <=tempmax-1; i++)
		{
			if (a[i] - b[i] == 1)
				return false;
		}	
	}

	return true;
}

void print_vector(vector<int> aa)
{
	for (int i = 0; i < aa.size(); i++)
	{
		cout << aa[i] << " ";
	}
	cout << endl;
}
