#ifndef SAMPLING_H
#define SAMPLING_H

#include <iostream>
#include <boost/random.hpp>
#include <math.h>
#include <algorithm>
#include <vector>
#include <cstdlib>

using namespace std;

typedef boost::mt19937                     ENG;    // Mersenne Twister
typedef boost::normal_distribution<double> DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> GEN;    // Variate generator
 
class Sampler
{
public:
	Sampler(){}
	double normrnd(double mu, double sigma)
	{
		boost::normal_distribution<double> dist(mu,sigma);
		boost::variate_generator<ENG,boost::normal_distribution<double> > gen(eng,dist);
		return gen();
	}
	
	double gamrnd(double shape, double scale)
	{
		boost::gamma_distribution<double> dist(shape,scale);
		boost::variate_generator<ENG,boost::gamma_distribution<double> > gen(eng,dist);
		return gen();
	}
	
	double invgamrnd(double shape, double scale)
	{
		if (isnan(scale) or scale < pow(10,-6))
			scale = pow(10,-6);
		if (isinf(scale) or scale > pow(10,6))
			scale = pow(10,6);
		boost::gamma_distribution<double> dist(shape,1/scale);
		boost::variate_generator<ENG,boost::gamma_distribution<double> > gen(eng,dist);
		return 1/gen();
	}

	int binornd(int n, double p)
	{
		boost::binomial_distribution<int> dist(n,p);
		boost::variate_generator<ENG,boost::binomial_distribution<int> > gen(eng,dist);
		return gen();
	}	

	double trnormrnd(double x0, double mu, double sigma2, double a)
	{
		boost::uniform_real<double> uni_dist(0,1);
		boost::variate_generator<ENG,boost::uniform_real<double> > gen(eng,uni_dist);
		int Nsamples = 10;
		 //truncated upper
		double b = 100000; //truncated lower
		vector<double> x(Nsamples,0);
		double sigma = sqrt(sigma2);
		double az = (a-mu)/sigma;
		double bz = (b-mu)/sigma;
		double z = (x0-mu)/sigma;
		for (int i = 0; i < Nsamples; i++)
		{
			double logy = -0.5*pow(z,2) + log(gen());
			double cz = max(az, -sqrt(-2*logy));
			double dz = sqrt(-2*logy);
			z = cz + (dz-cz) * gen();
			x[i] = mu + sigma * z;
		}
		return x[Nsamples-1];
	}

	double uniform_01()
	{
		boost::uniform_real<double> uni_dist(0,1);
		boost::variate_generator<ENG,boost::uniform_real<double> > gen(eng,uni_dist);
		return gen();
	}

	ENG  eng;
};




#endif
