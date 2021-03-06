#include "Info.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <algorithm>
#include "sampling.h"
#include "basic_math.h"
using namespace std;

void print_v(vector<double> v);
void print_vv(vector<vector<double> > v);
void print_v_int(vector<int> v);
void print_vv_int(vector<vector<int> > v);
void vectorcptoMat(vector<vector<double> > &a, MatrixXd &b);
bool vector_matchX(MatrixXd X, int a, int b);
void connect_region_recur(set<vector<int> > region_in, set<vector<int> >& region_out);
void connect_region(set<vector<int> > region_in, set<vector<int> >& region_out);
void groupiso(vector<vector<int> > &transloc, vector<int> &transgroup, int &max_group_idx);
bool sort1( const vector<int>& v1, const vector<int>& v2 ) {
    return v1[0] < v2[0];
};
bool sort1r( const vector<int>& v1, const vector<int>& v2 ) {
    return v1[0] > v2[0];
};
bool sort2( const vector<int>& v1, const vector<int>& v2 ) {
    return v1[1] < v2[1];
};
bool sort3( const vector<int>& v1, const vector<int>& v2 ) {
    return v1[2] < v2[2];
};
bool comparator(const pair<double,int> &l, const pair<double,int> &r)
{
	return l.first < r.first;
};

void Info::run()
{
	int N = X.rows();
	int M = X.cols();
	if (M > 1 and N > 1)
	{
		preprocess();
		vector<MatrixXd> XTX_v;
		vector<MatrixXd> XT_v;
		vector<vector<int> > slt_idx_v;
		vector<vector<int> > beta_map_idx_v;
		vector<int> all_slt_idx;
		VectorXd trans_length;
		calsamplespecificpar(XTX_v, XT_v, slt_idx_v, beta_map_idx_v, all_slt_idx, trans_length);
		
		trans_length = transL;
		int niso = all_slt_idx.size();
		N = X.rows();
		M = X.cols();
		double s_w_total = niso / M;
		int J = y.cols();

		double v0 = pow(10,-6);
		Sampler sample;
		
		double a1 = 3, a2 = 4;
		double a_sigma0 = pow(10,2), b_sigma0 = pow(10,-2);
		int iter = 500;
		
		MatrixXd s_rho_all = y;
		VectorXd s_tau = VectorXd::Ones(M);
		VectorXd s_beta = VectorXd::Zero(M);
		VectorXd s_R = VectorXd::Zero(M);
		MatrixXd s_beta_all = MatrixXd::Zero(M,J);
		MatrixXd s_R_all = MatrixXd::Zero(M,J);
		VectorXd s_gamma = VectorXd::Ones(M);
		VectorXd s_z = v0 * VectorXd::Ones(M);
		
		double s_sigma2 = 1E-5;
		double s_w = niso / M;
		vector<MatrixXd > samples_beta(iter,MatrixXd(M,J));
		vector<MatrixXd > samples_R(iter,MatrixXd(M,J));
		vector<VectorXd > samples_z(iter,VectorXd(M));

		for (int i = 0; i < iter; i++)
		{
			for (int j = 0; j < J; j++)
			{
				MatrixXd Xt = X_v[j];
				VectorXd yt = y_v[j];
				VectorXd Rt = R_v[j];
				VectorXd Lt = L_v[j];

				vector<int> beta_map_idxt = beta_map_idx_v[j];
				vector<int> slt_idxt = slt_idx_v[j];
				int Nt = Xt.rows();
				int Mt = Xt.cols();
				VectorXd s_betat = VectorXd::Zero(Mt);
				VectorXd s_Rt = VectorXd::Zero(Mt);
				VectorXd s_rho = yt;
				VectorXd s_gammat = v0 * VectorXd::Ones(Mt);
				if (Mt > 1)
				{
					if (i == 0)
					{
						for (int ii = 0; ii < slt_idxt.size(); ii++)
							s_gammat(slt_idxt[ii]) = 1;
					}
					else
					{
						s_gammat = VectorXd::Zero(Mt);
						s_betat = VectorXd::Zero(Mt);
						for (int ii = 0; ii < beta_map_idxt.size(); ii++)
						{
							s_gammat(ii) = s_gamma(beta_map_idxt[ii]);
							s_betat(ii) = s_beta_all(beta_map_idxt[ii],j);
						}
					}

					MatrixXd A;
					MatrixXd diag_gamma;
					vector<double> diag_v(Mt,0);
					for (int k = 0; k < Mt; k++)
						diag_v[k] = s_sigma2/s_gammat(k);
					diag(diag_v,diag_gamma);
					
					A = XTX_v[j] + diag_gamma;
					MatrixXd invA = A.inverse();
					MatrixXd mu_betaM = invA * XT_v[j] * s_rho;
					assert(mu_betaM.cols()==1);

					VectorXd mu_beta(Map<VectorXd>(mu_betaM.data(),mu_betaM.rows()*mu_betaM.cols()));
					MatrixXd invsigma_beta = 1/s_sigma2 * A;
					vector<int> s_beta_idx = slt_idxt;
					for (int k = 0; k < Mt;k++)
					{
						int countMt = 0;
						for (int kk = 0; kk < Nt; kk++)
							countMt += Xt(kk,k);
						if (countMt > 0)
						{
							vector<int>::iterator it = find(s_beta_idx.begin(),s_beta_idx.end(),k);
							if (it == s_beta_idx.end())
								s_beta_idx.push_back(k);
						}
					}					

					double p = 0;
					for (int k = 0; k < 10; k++)
					{
						for (int ii = 0; ii < s_beta_idx.size(); ii++)
						{
							int idx = s_beta_idx[ii];
							//idx = ii;
							VectorXd s_betatemp = s_betat;
							s_betatemp(idx) = 0;
							double vi2 = 1/invsigma_beta(idx,idx);
							VectorXd btemp,ctemp;
							btemp = mu_beta - s_betat;
							btemp(idx) = 0;
							mattoVecR(invsigma_beta,ctemp,idx);
							ctemp(idx) = 0;
							double ui = mu_beta(idx) + vi2 * btemp.dot(ctemp);
							double mu_th = get_muth(p, idx, s_beta);
							double betai = sample.trnormrnd(s_betat(idx),ui,vi2,mu_th);
							s_betat(idx) = betai;
						}
					}

					s_Rt = VectorXd::Zero(s_betat.size());
					for (int ii = 0; ii < Nt; ii++)
					{
						vector<int> idx_include;
						
						for (int jj = 0; jj < Mt; jj++)
						{
							if (Xt(ii,jj) > 0)
							{
								idx_include.push_back(jj);
							}
						}
						VectorXd s_beta_include(idx_include.size());
						for (int jj = 0; jj < idx_include.size(); jj++)
						{
							s_beta_include(jj) = s_betat(idx_include[jj]);
						}
						
						if (idx_include.size()>0)
						{
							VectorXd s_beta_include_frac = s_beta_include / s_beta_include.sum();
							double ReadsAssign = Rt(ii);
							for (int kk = 0; kk < idx_include.size();kk++)
							{
								s_Rt(idx_include[kk]) += ReadsAssign * s_beta_include_frac(kk);
							}
						}
					}

				}
				else if (Mt == 1)
				{
					s_betat(0) = Rt.sum()/trans_length(0);
					s_Rt(0) = Rt.sum();
				}
				
				for (int ii = 0; ii < beta_map_idxt.size(); ii++)
				{
					s_beta_all(beta_map_idxt[ii],j) = s_betat(ii);
					s_R_all(beta_map_idxt[ii],j) = s_Rt(ii);
					
				}
			}
			samples_beta[i] = s_beta_all;
			samples_R[i] = s_R_all;

			int Nexp_total = sample.binornd(M,s_w_total);
			
			if (Nexp_total > niso + 2)
				Nexp_total = min(niso + 2,M);
			if (Nexp_total < niso)
				Nexp_total = niso;

			VectorXd meanbetaall(s_R_all.rows());
			for (int ii = 0; ii < meanbetaall.size(); ii++)
			{
				for (int jj = 0; jj < s_R_all.cols(); jj++)
					meanbetaall(ii) += s_R_all(ii,jj);
			}
			
			VectorXd p0_z = VectorXd::Ones(meanbetaall.size())-meanbetaall/meanbetaall.sum();
			vector<pair<double,int> > sort_p0;
			for (int ii = 0; ii < p0_z.size(); ii++)
			{
				pair<double, int> s_pair = std::make_pair(p0_z(ii),ii);
				sort_p0.push_back(s_pair);
			}

			std::sort(sort_p0.begin(),sort_p0.end(),comparator);

			s_z = v0 * VectorXd::Ones(M);
			int Nz_select = 0;
			int prev_Nz = 0;
			int freez = 0;
			while (Nz_select < Nexp_total)
			{
				for (int ii = 0; ii < sort_p0.size(); ii++)
				{
					double rand_p0 = sample.uniform_01();
					if (rand_p0 > sort_p0[ii].first)
					{
						s_z(sort_p0[ii].second) = 1;
						meanbetaall(sort_p0[ii].second) = 0;
						Nz_select++;
						if (Nz_select >= Nexp_total)
							break;
					}
				}
				if (prev_Nz == Nz_select)
					freez++;
				else
				{
					prev_Nz = Nz_select;
					p0_z = VectorXd::Ones(meanbetaall.size())-meanbetaall/meanbetaall.sum();
					sort_p0.clear();
					for (int ii = 0; ii < p0_z.size(); ii++)
					{
						pair<double, int> s_pair = std::make_pair(p0_z(ii),ii);
						sort_p0.push_back(s_pair);
					}
					std::sort(sort_p0.begin(),sort_p0.end(),comparator);
				}
				if (freez == 5)
					break;
			}
			while (Nz_select < Nexp_total)
			{
				for (int ii = 0; ii < M; ii++)
				{
					if (s_z(ii) == v0)
					{
						double rand_z = sample.uniform_01();
						if (rand_z > 0.5)
						{
							s_z(ii) = 1;
							Nz_select++;
						}
					}
				}
			}

			samples_z[i] = s_z;
				
			s_sigma2 = 1e-5; // implement sigma sampling later
			for (int m = 0; m < M; m++)
			{
				VectorXd s_beta_m = s_beta_all.row(m);
				s_tau(m) = sample.invgamrnd(a1+0.5*J,a2+s_beta_m.dot(s_beta_m)/2/s_z(m));
				if (s_tau(m) > 10)
					s_tau(m) = 10;
				s_gamma(m) = s_tau(m) * s_z(m);
			}
			

		}
				
		int burnin = iter - 200;
		vector<double> final_z(M,0);
		vector<double> final_beta(M,0);
		vector<double> final_R(M,0);
		beta_frac_low_high = MatrixXd::Zero(M,2);

		//esimate beta and R
		
		set<double> beta_helper;
		set<double> R_helper;
		for (int ii = 0; ii < M; ii++)
		{
			beta_helper.clear();
			R_helper.clear();
			for (int i = burnin; i < iter; i++)
			{
				double beta_sum = 0;
				double R_sum = 0;
				for (int j = 0; j < J; j++)
				{
					beta_sum += samples_beta[i](ii,j);
					R_sum += samples_R[i](ii,j);
				}
				beta_helper.insert(beta_sum/J);
				R_helper.insert(R_sum/J);
			}
			if (beta_helper.size()%2 == 1)
			{
				set<double>::iterator it = beta_helper.begin();
				std::advance(it,(beta_helper.size()-1)/2);
				final_beta[ii] = *it;
				it = R_helper.begin();
				std::advance(it,(R_helper.size()-1)/2);
				final_R[ii] = *it;
			}
			else
			{
				set<double>::iterator it = beta_helper.begin();
				std::advance(it,(beta_helper.size())/2 - 1);
				double beta1 = *it;
				std::advance(it,1);
				double beta2 = *it;
				final_beta[ii] = (beta1+beta2)/2;

				it = R_helper.begin();
				std::advance(it,(R_helper.size())/2 - 1);
				double R1 = *it;
				std::advance(it,1);
				double R2 = *it;
				final_R[ii]= (R1+R2)/2;
			}
			set<double>::iterator it = beta_helper.begin();
			std::advance(it,round(0.025*beta_helper.size()));
			beta_frac_low_high(ii,0) = (*it)/final_beta[ii];
			it = beta_helper.begin();
			std::advance(it,beta_helper.size()-1-round(0.025*beta_helper.size()));
			beta_frac_low_high(ii,1) = (*it)/final_beta[ii];
		}


		//estimate z
		for (int ii = 0; ii < M; ii++)
		{
			double sumz = 0;
			for (int i = burnin; i < iter; i++)
			{
				beta_helper.insert(samples_beta[i](ii));
				sumz+= samples_z[i](ii);
			}
			final_z[ii] = (double)sumz / (double) (iter-burnin);
		}

		VectorXd final_beta_v = VectorXd::Map(&final_beta[0],final_beta.size());
		VectorXd FPKM, s_reads;
		calFPKM(final_beta_v,FPKM,s_reads);
		for (int ii = 0; ii < FPKM.size(); ii++)
		{
			bool select = final_z[ii] >= 0.5 and final_R[ii] > 12;
			if (select)
				final_isoidx.push_back(1);
			else
				final_isoidx.push_back(0);
			final_FPKM.push_back(FPKM(ii));
		}
	}

}

void Info::calsamplespecificpar(vector<MatrixXd> & XTX, vector<MatrixXd> & XT, vector<vector<int> > &slt_idx_v, vector<vector<int> > &beta_map_idx_v, vector<int> &all_slt_idx, VectorXd &trans_length)
{
	slt_idx_v.clear();
	MatrixXd X_R = X;
	for (int i = 0; i < X.rows(); i++)
	{
		double sumR = 0;
		for (int j = 0; j < R.cols(); j++)
			sumR += R(i,j);
		for (int j = 0; j < X.cols(); j++)
			X_R(i,j) = X(i,j) * sumR;
	}
	VectorXd sum_X_R1;
	matSum(X_R,sum_X_R1,1);
	VectorXd all_slt_idxv = VectorXd::Zero(X.cols());

	//for each sample...
	for (int j = 0; j < R.cols(); j++)
	{
		vector<int> idx; // segment idx included in sample j

		for (int i = 0; i < R.rows(); i++)
		{
			if (R(i,j) > 0)
			{
				idx.push_back(i);				
			}
		}

		vector<int> beta_idx;
		//for each transcript
		for (int i = 0; i < X.cols(); i++)
		{
			double sum_idxc = 0;
			double sum_allc = 0;
			for (int ii = 0; ii < X.rows(); ii++)
				sum_allc += X_R(ii,i);
			for (int ii = 0; ii < idx.size(); ii++)
				sum_idxc += X_R(idx[ii],i);
			if (sum_idxc/sum_allc > 0.9999)
				beta_idx.push_back(i);
		}
		VectorXd sum_X_R_idx;

		MatrixXd X1 = MatrixXd::Zero(idx.size(),beta_idx.size());
		MatrixXd Xm1 = MatrixXd::Zero(idx.size(),beta_idx.size());
		VectorXd y1 = VectorXd::Zero(idx.size());
		VectorXd R1 = VectorXd::Zero(idx.size());
		VectorXd L1 = VectorXd::Zero(idx.size());
		VectorXd sum_Xm1 = VectorXd::Zero(idx.size());
		trans_length = VectorXd::Zero(beta_idx.size());

		for (int i = 0; i < idx.size(); i++)
		{
			for(int ii = 0; ii < beta_idx.size(); ii++)
			{
				X1(i,ii) = X(idx[i],beta_idx[ii]);
				Xm1(i,ii) = X1(i,ii) * R(idx[i],j);
				sum_Xm1(i) += Xm1(i,ii);
				trans_length(ii) += X1(i,ii) * L(idx[i]);		
			}
			R1(i) = R(idx[i],j);
			L1(i) = L(idx[i]);
			y1(i) = R(idx[i],j)/L(idx[i]);
			
		}
		
		MatrixXd XT1 = X1.transpose();
		XT.push_back(XT1);
		MatrixXd XTX1 = XT1 * X1;
		XTX.push_back(XTX1);
		X_v.push_back(X1);
		y_v.push_back(y1);
		R_v.push_back(R1);
		L_v.push_back(L1);
		beta_map_idx_v.push_back(beta_idx);

		vector<int> current_idx;
		current_idx.clear();
		vector<int> new_idx;
		new_idx.clear();
		
		getIsoIdxForX(current_idx, new_idx, Xm1, sum_Xm1, trans_length);
		
		
		for (int i = 0; i < new_idx.size(); i++)
		{
			all_slt_idxv(beta_idx[new_idx[i]]) = 1;
		}
		slt_idx_v.push_back(new_idx);
		
	}
	
	for (int i = 0; i < all_slt_idxv.size(); i++)
	{
		if (all_slt_idxv(i) == 1)
			all_slt_idx.push_back(i);
	}
	
		
}

void Info::preprocess() // calculate unique of X
{
	vector<int> X_label;
	vector<vector<double> > y1;
	vector<vector<double> > X1;
	vector<double> L1;
	vector<vector<double> > R1;
	vector<double> lwv;
	X_label.resize(X.rows(),0);

	for (int i = 0; i < X.rows(); i++)
	{
		double count = 0;
		for (int j = 0; j < X.cols(); j++)
		{
			count+=X(i,j);
		}
		if (count == 0)
			X_label[i] = 1;
	}
	for (int i = 0; i < X.rows(); i++)
	{
		if (X_label[i] == 0)			
		{
			X_label[i] = 1;
			vector<double> sumR;
			for (int j = 0; j < R.cols(); j++)
			 	sumR.push_back(R(i,j));
			double sumL = L(i);
			for (int j = i+1; j < X.rows(); j++)
			{
				if (X_label[j] == 0)
				{
					if (vector_matchX(X,i,j))
					{
						X_label[j] = 1;
						for (int jj = 0; jj < sumR.size(); jj++)
							sumR[jj] += R(j,jj);
						sumL+=L(j);
					}
					
				}
			}
			vector<double> s_X1(X.cols(),0);
			for (int j = 0; j < s_X1.size();j++)
				s_X1[j] = X(i,j);
			
			X1.push_back(s_X1);
			L1.push_back(sumL);
			R1.push_back(sumR);
			vector<double> s_y1;
			for (int j = 0; j < sumR.size(); j++)
				s_y1.push_back(sumR[j]/sumL);
			y1.push_back(s_y1);
			lwv.push_back(sumL/100);
		}
		
	}
	
	vectorcptoMat(X1,X);
	vectorcptoMat(y1,y);
	vectorcptoMat(R1,R);
	L = VectorXd::Map(&L1[0],L1.size());
	lw = VectorXd::Map(&lwv[0],lwv.size());
	X_m = MatrixXd::Zero(X.rows(),X.cols());
	sumX_m = VectorXd::Zero(X.rows());

	for (int i = 0; i < X.rows(); i++)
	{
		double sumRcount = 0;
		for (int j = 0; j < y.cols(); j++)
		{
			sumRcount += R(i,j);
		}	
		for (int j = 0; j < X.cols(); j++)
		{
			X_m(i,j) = sumRcount * X(i,j);
			sumX_m(i) += X_m(i,j);
		}
	}

	//get transL
	transL = VectorXd::Zero(X.cols());
	for (int i = 0; i < X.rows(); i++)
	{
		for (int j = 0; j < X.cols(); j++)
		{
			if (X(i,j) > 0)
				transL(j) += L(i);
		}
	}

}

void Info::getIsoIdxForX(vector<int> &current_idx, vector<int> &new_idx, MatrixXd &X_m, VectorXd &sumX_m, VectorXd &trans_length)
{
	vector<int> current_exon_idx(X_m.rows(),0);
	for (int i = 0; i < current_idx.size(); i++)
	{
		for (int j = 0; j < X_m.rows(); j++)
		{
			if (X_m(j,current_idx[i]) > 0)
				current_exon_idx[j] = 1;
		}
	}
	
	set<int> slt_exon_idx;
	for (int i = 0; i < X_m.rows(); i++)
	{
		if (current_exon_idx[i] == 0 and sumX_m(i) > 0)
			slt_exon_idx.insert(i);
	}
	set<int> slt_exon_idx_temp;
	for (int i = 0; i < X_m.cols(); i++)
		slt_exon_idx_temp.insert(i);
	while (!slt_exon_idx.empty())
	{
		vector<double> sumX_slt;
		matSum(X_m,slt_exon_idx,slt_exon_idx_temp,1,sumX_slt);
		int maxidx = maxidxV(sumX_slt);
		new_idx.push_back(maxidx);
		for (int i = 0; i < X_m.rows(); i++)
		{
			if (X_m(i,maxidx) > 0)
				slt_exon_idx.erase(i);
		}
	}


}


void vectorcptoMat(vector<vector<double> > &a, MatrixXd &b)
{
	b = MatrixXd::Zero(a.size(),a[0].size());
	for (int i = 0; i < a.size(); i++)
	{
		for (int j = 0; j < a[0].size(); j++)
		{
			b(i,j) = a[i][j];
		}
	}
}

void Info::modX_exon(vector<vector<int> > &paths)
{
	assert(paths.size()>0);
	X_exon = MatrixXd::Zero(paths[0].size(),paths.size());
	for (int i = 0; i < paths.size(); i++)
	{
		for (int j = 0; j < paths[i].size(); j++)
		{
			X_exon(j,i) = paths[i][j];
		}
	}
}

void Info::modX(vector<vector<int> > &paths)
{
	assert(paths.size()>0);
	if (paths.size()>1)
		single = false;
	X = MatrixXd::Zero(paths[0].size(),paths.size());
	for (int i = 0; i < paths.size(); i++)
	{
		for (int j = 0; j < paths[i].size(); j++)
		{
			X(j,i) = paths[i][j];
			if (paths[i][j] == 0)
				X(j,i) = 0;
		}
	}
}

void Info::mody(vector<vector<double> > &s_y)
{
	y = MatrixXd::Zero(s_y[0].size(),s_y.size());
	for (int i = 0; i < s_y[0].size(); i++)
	{
		for (int j = 0; j < s_y.size(); j++)
			y(i,j) = s_y[j][i];
	}
}

void Info::modR(vector<vector<double> > &s_R)
{
	R = MatrixXd::Zero(s_R[0].size(),s_R.size());
	for (int i = 0; i < s_R[0].size(); i++)
	{
		for (int j = 0; j < s_R.size(); j++)
			R(i,j) = s_R[j][i];
	}
}

void Info::modL(vector<double> &s_L)
{
	L = VectorXd::Map(&s_L[0],s_L.size());
}


void Info::clear()
{
	X.resize(0,0);
	X_m.resize(0,0); 
	sumX_m.resize(0); 
	X_exon.resize(0,0);
	y.resize(0,0);
	R.resize(0,0);
	L.resize(0);
	bias.resize(0);
	transL.resize(0);
	exonL.resize(0);
	lw.resize(0);
	label = "";
	single = true;
}

void Info::clear_info()
{
	X.resize(0,0);
	X_m.resize(0,0); 
	sumX_m.resize(0); 
	y.resize(0,0);
	R.resize(0,0);
	L.resize(0);
	bias.resize(0);
	transL.resize(0);
	exonL.resize(0);
	lw.resize(0);
	single = true;
}

double Info::get_muth(double p, int idx, VectorXd s_beta)
{
	MatrixXd mu_tempM = X * s_beta;
	VectorXd mu_temp;
	mattoVecC(mu_tempM,mu_temp,0);
	double min_v = 0;
	for (int i = 0; i < mu_temp.size(); i++)
	{
		if (X(i,idx) > 0)
		{
			if (min_v > mu_temp(i)-p)
				min_v = mu_temp(i)-p;
		}
	}
	return -min_v;
}

void Info::modexonbound(vector<vector<int> > &s_bound)
{
	exonbound = s_bound;
}

void Info::modbias(vector<double> &s_bias)
{
	bias = VectorXd::Map(&s_bias[0],s_bias.size());
}

void print_v(vector<double> v)
{
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}

void print_v_int(vector<int> v)
{
	for (int i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}

void print_vv(vector<vector<double> > v)
{
	for (int i = 0; i < v.size(); i++)
		print_v(v[i]);
}

void print_vv_int(vector<vector<int> > v)
{
	for (int i = 0; i < v.size(); i++)
		print_v_int(v[i]);
}

void Info::run_single()
{
	double reads = R.sum();
	double length = L.sum();
	if (reads > 12)
	{
		output_bool = true;
		final_FPKM.push_back(reads / length * pow(10,9) / totalReads);
		final_isoidx.push_back(0);
		beta_frac_low_high = MatrixXd::Zero(1,2);
		beta_frac_low_high(0,0) = 0.5;
		beta_frac_low_high(0,1) = 1.5;
	}
	else
		output_bool = false;

}

bool vector_matchX(MatrixXd X, int a, int b)
{
	for (int i = 0; i < X.cols(); i++)
	{
		if (X(a,i) != X(b,i))
			return false;
	}
	return true;
}

void Info::write(string outdir, int gene_idx)
{
	stringstream output_all, output;
	output << outdir << "/IntAPT.gtf";
	output_all << outdir << "/IntAPT_enumerate.gtf";
	string output_str, output_all_str;
	output_str = output.str();
	output_all_str = output_all.str();

	FILE *outfile;
	FILE *outfile_all;
	if (gene_idx == 1)
	{
		outfile = fopen(output_str.c_str(),"w");
		if (output_all_bool)
			outfile_all = fopen(output_all_str.c_str(),"w");
	}
	else
	{
		outfile = fopen(output_str.c_str(),"a+");
		if (output_all_bool)
			outfile_all = fopen(output_all_str.c_str(),"a+");
	}
	if (single) // or single
	{
		double singletransReads = R.sum()/R.cols();
		if (singletransReads > 12 and final_FPKM[0] > 1)
		{
			int idx = label.find_first_of("\t");
			string chr = label.substr(0,idx);
			int idx1 = label.find_last_of("\t");
			string strand = label.substr(idx1+1);
			stringstream outgene_stream;
			outgene_stream << "CUFF." << gene_idx;
			string outgene = outgene_stream.str();

			int i = 0;
			
			set<vector<int> > exon_region;
			getExonRegion(i,exon_region);
			set<vector<int> >::iterator it = exon_region.begin();
			int e_start = (*it)[0];
			it = exon_region.begin();
			for (int ii = 0; ii < exon_region.size()-1; ii++)
				std::advance(it,1);
			
			int e_end = (*it)[1];
			double translength_s = 0;
                        for (it = exon_region.begin();it!=exon_region.end();it++)
	                        translength_s += (*it)[1] - (*it)[0] + 1;
			
			stringstream trans_stream;
			trans_stream << outgene << "." << i+1;
			string outtrans = trans_stream.str();
			if (output_all_bool)
				fprintf(outfile_all, "%s\tIntAPT\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%.6f\"; frac \"%.6f\";\n",chr.c_str(), e_start, e_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), final_FPKM[i], 1.0);
			if (exon_region.size() > 1 and translength_s > 500)
				fprintf(outfile, "%s\tIntAPT\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%.6f\"; frac \"%.6f\";\n",chr.c_str(), e_start, e_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), final_FPKM[i], 1.0);

			int exonNum = 1;
			for (it = exon_region.begin(); it != exon_region.end(); it++)
			{
				int exon_start = (*it)[0];
				int exon_end = (*it)[1];
				if (output_all_bool)
					fprintf(outfile_all, "%s\tIntAPT\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%.6f\"; frac \"%.6f\";\n",chr.c_str(), exon_start, exon_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), exonNum, final_FPKM[i], 1.0);
				if (exon_region.size() > 1 and translength_s > 500)				
					fprintf(outfile, "%s\tIntAPT\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%.6f\"; frac \"%.6f\";\n",chr.c_str(), exon_start, exon_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), exonNum, final_FPKM[i], 1.0);

				exonNum++;
			}
		}
	}
	else
	{
		int idx = label.find_first_of("\t");
		string chr = label.substr(0,idx);
		int idx1 = label.find_last_of("\t");
		string strand = label.substr(idx1+1);
		stringstream outgene_stream;
		outgene_stream << "CUFF." << gene_idx;
		string outgene = outgene_stream.str();
		vector<int> final_isoidx_slt;
		
		for (int i = 0; i < final_isoidx.size(); i++)
		{
			if (final_isoidx[i] == 1)
				final_isoidx_slt.push_back(i);
		}

		if (final_isoidx_slt.size() > 5)
		{
			vector<vector<int> > transloc;
			vector<int> transgroup;
			int max_group_idx;
			vector<double> transFPM;
			vector<double> translength;
			//find exon start, end and transcript length
			for (int i = 0; i < final_isoidx_slt.size(); i++)
			{
				set<vector<int> > exon_region;
				getExonRegion(final_isoidx_slt[i],exon_region);
				set<vector<int> >::iterator it = exon_region.begin();
				set<vector<int> >::reverse_iterator itr = exon_region.rbegin();
				vector<int> transloc_v;
				transloc_v.push_back((*it)[0]);
				transloc_v.push_back((*itr)[1]);
				transloc.push_back(transloc_v);
				double translength_s = 0;
				for (it = exon_region.begin();it!=exon_region.end();it++)
					translength_s += (*it)[1] - (*it)[0] + 1;
				translength.push_back(translength_s);
				if (exon_region.size() <= 1 or translength_s < 500)
					final_isoidx[final_isoidx_slt[i]] = 0;
				transFPM.push_back(translength_s*final_FPKM[final_isoidx_slt[i]]);
			}
			
			groupiso(transloc,transgroup,max_group_idx);
			
			for (int i = 1; i <= max_group_idx; i++)
			{
				//search for the group idx
				set<vector<int> > sorted_list;
				for (int j = 0; j < transgroup.size(); j++)
				{
					if (transgroup[j] == i)
					{
						vector<int> sorted_list_s;
						sorted_list_s.push_back(transFPM[j]);
						sorted_list_s.push_back(j);
						sorted_list.insert(sorted_list_s);
					}
				}
				int groupsize = sorted_list.size();
				if (groupsize > 5)
				{
					set<vector<int> >::iterator it = sorted_list.begin();
					for (int j = 0; j < groupsize-5; j++)
					{
						final_isoidx[final_isoidx_slt[(*it)[1]]] = 0;
						std::advance(it,1);
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < final_isoidx_slt.size(); i++)
			{
				set<vector<int> > exon_region;
				getExonRegion(final_isoidx_slt[i],exon_region);
				set<vector<int> >::iterator it = exon_region.begin();
				double translength_s = 0;
				for (it = exon_region.begin();it!=exon_region.end();it++)
					translength_s += (*it)[1] - (*it)[0] + 1;
				if (exon_region.size() <= 1 or translength_s < 500)
					final_isoidx[final_isoidx_slt[i]] = 0;
			}
		}

		double sum_FPKM = 0;
		for (int i = 0; i < final_isoidx.size(); i++)
			if (final_isoidx[i] == 1)
				sum_FPKM += final_FPKM[i];
		vector<double> frac_FPKM;
		for (int i = 0; i < final_isoidx.size(); i++)
		{
            		if (final_isoidx[i] == 1)
				frac_FPKM.push_back(final_FPKM[i]/sum_FPKM);
			else
				frac_FPKM.push_back(0);
		}
		for (int i = 0; i < final_isoidx.size(); i++)
		{
			set<vector<int> > exon_region;
			getExonRegion(i,exon_region);
			set<vector<int> >::iterator it = exon_region.begin();
			int e_start = (*it)[0];
			it = exon_region.begin();
			for (int ii = 0; ii < exon_region.size()-1; ii++)
				std::advance(it,1);
			
			int e_end = (*it)[1];
			stringstream trans_stream;
			trans_stream << outgene << "." << i+1;

			string outtrans = trans_stream.str();
			if (output_all_bool)
				fprintf(outfile_all, "%s\tIntAPT\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%.6f\"; frac \"%.6f\";\n",chr.c_str(), e_start, e_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), final_FPKM[i], frac_FPKM[i]);
			if (final_isoidx[i] == 1)
				fprintf(outfile, "%s\tIntAPT\ttranscript\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%.6f\"; frac \"%.6f\";\n",chr.c_str(), e_start, e_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), final_FPKM[i], frac_FPKM[i]);

			int exonNum = 1;
			for (it = exon_region.begin(); it != exon_region.end(); it++)
			{
				int exon_start = (*it)[0];
				int exon_end = (*it)[1];
				if (output_all_bool)
					fprintf(outfile_all, "%s\tIntAPT\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%.6f\"; frac \"%.6f\";\n",chr.c_str(), exon_start, exon_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), exonNum, final_FPKM[i], frac_FPKM[i]);
				if (final_isoidx[i] == 1)
					fprintf(outfile, "%s\tIntAPT\texon\t%d\t%d\t1000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%d\"; FPKM \"%.6f\"; frac \"%.6f\";\n",chr.c_str(), exon_start, exon_end, strand.c_str() ,outgene.c_str(), outtrans.c_str(), exonNum, final_FPKM[i], frac_FPKM[i]);

				exonNum++;
			}
		}
	}
	if (output_all_bool)
		fclose(outfile_all);
	fclose(outfile);
}

void Info::getExonRegion(int isoidx, set<vector<int> >& exon_region)
{
	set<vector<int> > exon_temp;
	for (int i = 0; i < X_exon.rows(); i++)
	{
		if (X_exon(i,isoidx) > 0)
		{
			vector<int> s_exon_temp;
			s_exon_temp.push_back(exonbound[i][0]);
			s_exon_temp.push_back(exonbound[i][1]);
			exon_temp.insert(s_exon_temp);
		}
	}
	connect_region_recur(exon_temp, exon_region);


}

void connect_region_recur(set<vector<int> > region_in, set<vector<int> >& region_out)
{
	int quit = 0;
	set<vector<int> > temp = region_in;
	while (quit == 0)
	{	
		region_out.clear();
		connect_region(temp,region_out);
		if (region_out.size() == temp.size())
			quit = 1;
		else
			temp = region_out;
	}
}

void connect_region(set<vector<int> > region_in, set<vector<int> >& region_out)
{
	region_out.clear();
	int count = 1, jump = 0;
	vector<vector<int> > region_in_v;
	for (set<vector<int> >::iterator it = region_in.begin(); it != region_in.end(); it++)
		region_in_v.push_back(*it);
	vector<int> zeros;
	zeros.push_back(0);
	zeros.push_back(0);
	region_in_v.push_back(zeros);
	for (int i = 0; i < region_in_v.size()-1; i++)
	{
		int e_end = region_in_v[i][1];
		int e_start = 0;
		
		e_start = region_in_v[i+1][0];
		
		if (jump == 0)
		{
			if (e_start - e_end == 1)
			{
				int temp1 = region_in_v[i][0];
				int temp2 = region_in_v[i][1];
				
				if (temp1 > region_in_v[i+1][0])
					temp1 = region_in_v[i+1][0];
				if (temp2 < region_in_v[i+1][1])
					temp2 = region_in_v[i+1][1];

				vector<int> temp;
				temp.push_back(temp1);
				temp.push_back(temp2);
				region_out.insert(temp);
				region_in_v[i+1] = temp;
				jump = 1;
			}
			else
			{
				region_out.insert(region_in_v[i]);
			}
		}
		else
			jump = 0;
	}
}

void Info::calFPKM(VectorXd s_beta, VectorXd &FPKM, VectorXd &s_reads)
{
	double totalReads_mean = 0;
	for (int i = 0; i < TOTALREADSJ.size(); i++)
		totalReads_mean += (double) TOTALREADSJ[i];
	totalReads_mean /= TOTALREADSJ.size();
	s_reads = VectorXd::Zero(s_beta.size());
	VectorXd R_sum = VectorXd::Zero(R.rows());
	for (int i = 0; i < R.rows(); i++)
	{
		for (int j = 0; j < R.cols(); j++)
			R_sum(i) += R(i,j);
	}
	for (int i = 0; i < R_sum.size(); i++)
	{
		double sum_fpkm = 0;
		for (int j = 0; j < s_reads.size(); j++)
		{
			if (X(i,j) > 0)
			{
				sum_fpkm += s_beta(j);
			}
		}

		for (int j = 0; j < s_reads.size(); j++)
		{
			if (X(i,j) > 0)
			{
				s_reads(j) += s_beta(j) / sum_fpkm * R_sum(i);
			}
		}
	}
	FPKM = s_reads.array() / transL.array() * pow(10,9) / totalReads;
}

void groupiso(vector<vector<int> > &transloc, vector<int> &transgroup, int &max_group_idx)
{
    transgroup.clear();
    if (transloc.size() == 1)
    {
        transgroup.push_back(1);
        return;
    }
    sort(transloc.begin(),transloc.end(),sort2);
    sort(transloc.begin(),transloc.end(),sort1);
  

    int current_p1 = transloc[0][0];
    int current_p2 = transloc[0][1];
    int grp_idx = 1;
    transgroup.push_back(grp_idx);
    for (int i = 1; i < transloc.size(); i++)
    {
        if (transloc[i][0] <= current_p2)
        {
            if (transloc[i][1] > current_p2)
                current_p2 = transloc[i][1];
        }
        else
        {
            current_p1 = transloc[i][0];
            current_p2 = transloc[i][1];
            grp_idx++;
        }
        transgroup.push_back(grp_idx);
    }
    max_group_idx = transgroup[transgroup.size()-1];
    vector<int> tempgroup = transgroup;
    for (int i = 0; i < transloc.size(); i++)
    {
        for (int j = 0; j < transloc.size(); j++)
        {
            if (transloc[j][2] == i+1)
                transgroup[i] = tempgroup[j];
        }
    }

}
