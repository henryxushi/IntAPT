#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include "Info.h"
using namespace std;

const int MAX_NUM_GENE = 100000;
const int MAX_NUM_PATH = 20;
const double MAX_INTRON_COV = 0.95;

struct Seg
{
	int node;
	vector<int> visited;
};

struct PathVar
{
public:
	PathVar();
	double maxniso;
	double minjuncfrac;
	double minparentjuncfrac;
	double minadjjuncfrac;
	double minintronretentioncvgfrac;
	double minintronretentionlen;
};

class Instance
{
public:
	Instance(){}
	string location;
	string strand;
	vector<vector<int> > exonbound;
	vector<vector<double> > exonAbundance;
	vector<vector<int> > SegConf;
	vector<vector<int> > SegConfPET;
	vector<int> ReadLen;
	vector<int> L;
	int NumofExon;
	int NumofSegs;
	int NumofSamples;
	double maxExonCov;
	double frag_mu;
	double frag_sigma;
	PathVar pathvar;
	vector<vector<int> > SegAbundance;
	vector<int> connSeg; // find SGTypes with connection, count number of exons
	vector<vector<int> > PETconnect;
	vector<int> PET_switch;
	vector<vector<vector<int> > > PETreads;
	vector<int> SegAbundancePET;
	vector<set<int> > inc_exon;
	vector<double> totalreads;
	vector<double> exonp0;
	vector<vector<double> > exoncvg;
	vector<vector<double> > junccount;
	vector<double> sumexoncvg;
	
	int NumofPaths;
	
	void process_inst(Info& info, string outputfile);
	void process_inst_se(Info& info, string outputfile);
	void correct_y(vector<vector<int> > &paths, ofstream &outfile,Info &info);
	void remove_low();
	void getmap(map<int,set<int> > &exonmap, double &th);
	int enumerate_path(map<int,set<int> > &exonmap, int not_source_list[], int not_sink_list[], ofstream &outfile, vector<vector<int> > &paths);
	void mod_X(vector<int> &path, vector<int> &segpath, vector<int> &path_vec, ofstream &outfile);
	void mod_y(map<int, vector<int> > &iscontained, vector<double> &SegAbundance1, ofstream &outfile, int &output, double &th1, double &min_abun);
	//void mod_y_se(map<int, vector<int> > &iscontained, vector<double> &SegAbundance1, vector<double> &proc_R, vector<double> &proc_L, ofstream &outfile, int &output, double &th1, double &min_abun);
	void mod_y_se(map<int, vector<int> > &iscontained, vector<vector<double> > &SegAbundance1, vector<vector<double> > &proc_R, vector<vector<double> > &proc_L, ofstream &outfile, int &output, double &th1, double &min_abun);

	void findcontain(map<int, vector<int> > &iscontained); // contain: index start from 1
	int calintronlength(int &s1, int &s2, vector<int> &path, int &l1, int &l2, int &l3, int label, int &bias);
	void intronlength_helper(int &pos1, int &pos2, vector<int> &path, int &intron_l, int &l3, int &bias);
	void exonlength_plus(int &pos1, int &pos2, vector<int> &path, int &intron_l,int &l3);//calculate exon length between two points
	void exonlength_reverse(int &pos1, int &pos2, vector<int> &path, int &intron_l, int &l3);
	int segStart(int &s1, int &s2);
	void segStart(int &s1, int &s2, set<int> &s3);
	int segStart(int &s1);
	void cal_junccount();
	void cal_sumexoncvg();
	double cal_connreads(int i, int j);
	void findchildren(int c, vector<int> &nodes);
	void refineparents(map<int,set<int> > &exonmap,map<int,set<int> > &rev_exonmap);
	void findparent(int c, set<int> &parents, vector<int> &remove_nodes);
};

#endif
