#ifndef INFO_H
#define INFO_H

#include <vector>
#include <string>
#include <set>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

class Info
{
public:
	Info(){
		single = true;
	}
	// general parameters
	MatrixXd X;
	MatrixXd X_m; // X with length; preference of select: long segment to short segment
	VectorXd sumX_m; // sum X horizontally
	MatrixXd X_exon;
	MatrixXd y;
	MatrixXd y_m;
	MatrixXd R;
	VectorXd L;
	VectorXd bias;
	VectorXd transL;
	VectorXd exonL;
	VectorXd lw; //length weight of segment

	// parameters specific to each sample
	vector<MatrixXd> X_v;
	vector<VectorXd> y_v;
	vector<VectorXd> R_v;
	vector<VectorXd> L_v;
	


	vector<vector<int> > exonbound;
	string label;
	vector<int> final_isoidx;
	vector<double> final_FPKM;
	MatrixXd beta_frac_low_high;
	MatrixXd beta_frac_low;
	MatrixXd beta_frac_high;
	int totalReads;
	vector<int> TOTALREADSJ;
	int infoidx;
	void run();
	void run_single();
	void clear();
	void clear_info();
	bool single;
	bool bias_flag;
	void preprocess();
	void calsamplespecificpar(vector<MatrixXd> & XTX, vector<MatrixXd> & XT, vector<vector<int> > &slt_idx_v, vector<vector<int> > &beta_map_idx_v, vector<int> &all_slt_idx, VectorXd &trans_length);
	void modX_exon(vector<vector<int> > &paths);
	void modX(vector<vector<int> > &paths);
	void mody(vector<vector<double> > &s_y);
	void addy(vector<double> &s_y);
	void modR(vector<vector<double> > &s_R);
	void addR(vector<double> &s_R);
	void modL(vector<double> &s_L);
	void addL(double s_L);
	void modexonbound(vector<vector<int> > &s_bound);
	void modbias(vector<double> &s_bias);
	void getIsoIdxForX(vector<int> &current_idx, vector<int> &new_idx, MatrixXd &X_m, VectorXd &sumX_m, VectorXd &trans_length); //choose the isoform to explain all reads
	void calFPKM(VectorXd s_beta, VectorXd& FPKM, VectorXd& s_reads);
	void write(string outdir, int idx);
	void getExonRegion(int isoidx, set<vector<int> >& exon_region); //get the exon boundaries for isoform isoidx
	double get_muth(double p, int idx, VectorXd s_beta);
	bool output_bool;
	bool valid;
	bool output_all_bool;
};

#endif
