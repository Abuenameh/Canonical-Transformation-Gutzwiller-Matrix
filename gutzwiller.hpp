/* 
 * File:   gutzwiller.hpp
 * Author: Abuenameh
 *
 * Created on August 10, 2014, 10:45 PM
 */

#ifndef GUTZWILLER_HPP
#define	GUTZWILLER_HPP

#include <complex>
#include <vector>
#include <iostream>
#include <array>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef complex<double> doublecomplex;

#define L 50
#define nmax 11
#define dim (nmax+1)

typedef Matrix<doublecomplex, dim, dim> Hamiltonian;
typedef Matrix<doublecomplex, dim, 1> SiteWaveFunction;
typedef array<SiteWaveFunction, L> WaveFunction;
typedef Matrix<double, dim, 1> SiteWaveFunctionAbs;
typedef array<SiteWaveFunctionAbs, L> WaveFunctionAbs;
typedef array<double, L> Parameter;
typedef array<array<double, 6>, L> EnergyComponents;

template<class T>
complex<T> operator~(const complex<T> a) {
	return conj(a);
}

struct funcdata {
    bool canonical;
//	vector<double> U;
    Parameter U;
    Parameter J;
	double mu;
//	vector<double> J;
    double theta;
    double Emin;
    vector<double> xmin;
//	double* U;
//	double mu;
//	double* J;
};

struct phase_parameters {
    double theta;
    bool canonical;
    double eps;
};

struct device_parameters {
	double* U;
	double mu;
	double* J;
    double theta;
};

inline int mod(int i) {
	return (i + L) % L;
}

inline double g(int n, int m) {
	return sqrt(1.0*(n + 1) * m);
}

extern vector<double> nu;

inline double eps(vector<double> U, int i, int j, int n, int m) {
	return n * U[i] - (m - 1) * U[j] + nu[i] - nu[j];
}

inline double eps(Parameter U, int i, int j, int n, int m) {
	return n * U[i] - (m - 1) * U[j] /*+ nu[i] - nu[j]*/;
}

//double Encfunc(unsigned ndim, const double *x, double *grad, void *data);
//double Ecfunc(unsigned ndim, const double *x, double *grad, void *data);
//double Ecfuncp(ProblemSetup *setup, double *x);

void createHamiltonian(Hamiltonian& H, WaveFunction& f, int i, Parameter& J, Parameter& U, double mu, double theta);
double Ecfunc(WaveFunction& f, Parameter& U, Parameter& J, double mu, double theta, double& Emin, WaveFunction& fmin, EnergyComponents& Ei);



#endif	/* GUTZWILLER_HPP */

