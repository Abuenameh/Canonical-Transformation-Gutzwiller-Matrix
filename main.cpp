/*
 * File:   main.cpp
 * Author: Abuenameh
 *
 * Created on August 6, 2014, 11:21 PM
 */

#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <nlopt.h>
#include <complex>
#include <iostream>
#include <queue>
//#include <thread>
#include <nlopt.hpp>

#include <boost/array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/date_time.hpp>
#include <boost/random.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/progress.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>

#include "mathematica.hpp"
#include "gutzwiller.hpp"
#include "PGSL.h"

using namespace std;

//using boost::lexical_cast;
using namespace boost;
using namespace boost::random;
using namespace boost::filesystem;
using namespace boost::posix_time;

//typedef boost::array<double, L> Parameter;

//template<typename T> void printMath(ostream& out, string name, T& t) {
//    out << name << "=" << ::math(t) << ";" << endl;
//}
//
//template<typename T> void printMath(ostream& out, string name, int i, T& t) {
//    out << name << "[" << i << "]" << "=" << ::math(t) << ";" << endl;
//}

double Encfunc(unsigned ndim, const double *x, double *grad, void *data) {
    return 0;
}

double Ecfunc(unsigned ndim, const double *x, double *grad, void *data) {
    return 0;
}

double M = 1000;
double g13 = 2.5e9;
double g24 = 2.5e9;
double delta = 1.0e12;
double Delta = -2.0e10;
double alpha = 1.1e7;

double Ng = sqrt(M) * g13;

double JW(double W) {
    return alpha * (W * W) / (Ng * Ng + W * W);
}

double JWij(double Wi, double Wj) {
    return alpha * (Wi * Wj) / (sqrt(Ng * Ng + Wi * Wi) * sqrt(Ng * Ng + Wj * Wj));
}

Parameter JW(Parameter W) {
    Parameter v;
    for (int i = 0; i < L; i++) {
        v[i] = W[i] / sqrt(Ng * Ng + W[i] * W[i]);
    }
    Parameter J;
    for (int i = 0; i < L - 1; i++) {
        J[i] = alpha * v[i] * v[i + 1];
    }
    J[L - 1] = alpha * v[L - 1] * v[0];
    return J;
}

double UW(double W) {
    return -2 * (g24 * g24) / Delta * (Ng * Ng * W * W) / ((Ng * Ng + W * W) * (Ng * Ng + W * W));
}

Parameter UW(Parameter W) {
    Parameter U;
    for (int i = 0; i < L; i++) {
        U[i] = -2 * (g24 * g24) / Delta * (Ng * Ng * W[i] * W[i]) / ((Ng * Ng + W[i] * W[i]) * (Ng * Ng + W[i] * W[i]));
    }
    return U;
}

boost::mutex progress_mutex;
boost::mutex points_mutex;

struct Point {
    int i;
    int j;
    double x;
    double mu;
};

double Efunc(unsigned ndim, const double *x, double *grad, void *data) {
    funcdata* parms = static_cast<funcdata*> (data);
    if (parms->canonical) {
        return Ecfunc(ndim, x, grad, data);
    } else {
        return Encfunc(ndim, x, grad, data);
    }
}

void norm2s(unsigned m, double *result, unsigned ndim, const double* x,
        double* grad, void* data);

//double norm(const vector<double> x, vector<double>& norms) {
//    const doublecomplex * f[L];
//    for (int i = 0; i < L; i++) {
//        f[i] = reinterpret_cast<const doublecomplex*> (&x[2 * i * dim]);
//    }
//
//    norms.resize(L);
//
//    //    double norm = 1;
//    for (int i = 0; i < L; i++) {
//        double normi = 0;
//        for (int n = 0; n <= nmax; n++) {
//            normi += norm(f[i][n]);
//        }
//        //        norm *= normi;
//        norms[i] = sqrt(normi);
//    }
//    //    return norm;
//    return 0;
//}

int min(int a, int b) {
    return a < b ? a : b;
}

int max(int a, int b) {
    return a > b ? a : b;
}

struct fresults {
    multi_array<double, 2>& fmin;
    multi_array<std::array<double, L>, 2>& fn0;
    multi_array<std::array<double, L>, 2>& fmax;
    multi_array<WaveFunction, 2>& f0;
    multi_array<WaveFunction, 2>& fth;
    multi_array<WaveFunction, 2>& f2th;
};

struct results {
    multi_array<double, 2 >& E0res;
    multi_array<double, 2 >& Ethres;
    multi_array<double, 2 >& E2thres;
    multi_array<double, 2 >& fs;
};

bool allClose(WaveFunction& psi, WaveFunction& phi, double eps) {
    for (int i = 0; i < L; i++) {
        if (!psi[i].cwiseAbs().isApprox(phi[i].cwiseAbs(), eps)) {
            return false;
        }
    }
    return true;
}

doublecomplex chop(doublecomplex x) {
    if (abs(x) < 1e-20) {
        return 0;
    } else {
        return x;
    }
}

#define MAXITER 500
#define NRUNS 50

static inline void progressbar(unsigned int x, unsigned int n, /*unsigned int y, unsigned int m,*/ unsigned int w = 50) {
    float ratio = x / (float) n;
    int c = ratio * w;

    printf("%3d%% [", (int)(ratio*100));
    for (int x=0; x<c; x++)
       printf("=");
    for (int x=c; x<w; x++)
       printf(" ");
    printf("]");
    printf("\033[G");
    fflush(stdout);
}

static inline void progressbar(int i, unsigned int x, unsigned int n, /*unsigned int y, unsigned int m,*/ unsigned int w = 50) {
    float ratio = x / (float) n;
    int c = ratio * w;

    printf("\033[%dB", i+1);
    printf("%3d: [", i+1);
    for (int x=0; x<c; x++)
       printf("=");
    for (int x=c; x<w; x++)
       printf(" ");
    printf("]");
    printf("\033[%dA", i+1);
    printf("\033[G");
    fflush(stdout);
}

void phasepoints2(int thread, Parameter& xi, phase_parameters parms, queue<Point>& points, fresults& fres, results& results) {

    mt19937 rng;
    rng.seed(time(NULL));
    uniform_real_distribution<> uni(-1, 1);

    Parameter W, U, J;

    double scale = 1;
    
    int numpoints = points.size();

    for (;;) {
        {
            boost::mutex::scoped_lock lock(progress_mutex);
                progressbar(thread, 0, NRUNS);
            //                ++progress;
        }
        Point point;
        {
            boost::mutex::scoped_lock lock(points_mutex);
            if (points.empty()) {
                break;
            }
            point = points.front();
            points.pop();
        }

        for (int i = 0; i < L; i++) {
            W[i] = xi[i] * point.x;
        }

        for (int i = 0; i < L; i++) {
            U[i] = UW(W[i]) / UW(point.x) / scale;
            J[i] = JWij(W[i], W[mod(i + 1)]) / UW(point.x) / scale;
            //            U[i] = 1;
            //            J[i] = JW(point.x)/UW(point.x);

            ////		U[i] = 0.1 * sqrt(i + 1);
            ////		J[i] = 0.1 * min(i + 1, mod(i + 1) + 1)
            ////			+ 0.2 * max(i + 1, mod(i + 1) + 1);
            //            U[i] = 1+0.2*uni(rng);
            //                        J[i] = point.x;
        }

        std::ofstream out("/Users/Abuenameh/Documents/Temp/mindata5.txt");

        WaveFunction f, fnew;

        Hamiltonian H;

        bool failed;
        int iter;
        double theta;

        funcdata data;
        data.Emin = DBL_MAX;
        data.xmin = vector<double>(2 * L * dim);

        double nan = numeric_limits<double>::quiet_NaN();
        SiteWaveFunction fnan;
        fnan.fill(nan);

        double E0min = DBL_MAX;
        WaveFunction f0min;
        EnergyComponents Ecomps;

        double Ethmin = DBL_MAX;
        WaveFunction fthmin;

        double E2thmin = DBL_MAX;
        WaveFunction f2thmin;

        fres.fmin[point.i][point.j] = nan;
        fres.fn0[point.i][point.j].fill(nan);
        fres.fmax[point.i][point.j].fill(nan);
        fres.f0[point.i][point.j].fill(fnan);
        results.E0res[point.i][point.j] = nan;
        fres.fth[point.i][point.j].fill(fnan);
        results.Ethres[point.i][point.j] = nan;
        fres.f2th[point.i][point.j].fill(fnan);
        results.E2thres[point.i][point.j] = nan;
        results.fs[point.i][point.j] = nan;

        for (int j = 0; j < NRUNS; j++) {

            SiteWaveFunction psii;
            //                    rng.seed(0);

            for (int i = 0; i < L; i++) {
                SiteWaveFunction fi;
                for (int n = 0; n <= nmax; n++) {
                    //                    fi[n] = doublecomplex(uni(rng), uni(rng)); //1;
                    fi[n] = uni(rng);
                }
                fi.normalize();
                fnew[i] = fi;
            }

            failed = false;
            iter = 0;
            theta = 0;
            do {
                for (int i = 0; i < L; i++) {
                    f[i] = fnew[i];
                }

                for (int i = 0; i < L; i++) {
                    createHamiltonian(H, fnew, i, J, U, point.mu, theta);
                    SelfAdjointEigenSolver<Hamiltonian> es(H);
                    if (es.info() != Success) {
                        cout << "Diagonalization failed: " << es.info() << endl;
                        failed = true;
                        break;
                    }
                    fnew[i] = es.eigenvectors().col(0).unaryExpr(std::ptr_fun(chop));
                }
                iter++;
            } while (!allClose(f, fnew, parms.eps) && !failed && iter < MAXITER);

            if (iter == MAXITER) {
                continue;
            }

            for (int i = 0; i < L; i++) {
                f[i] = fnew[i];
            }

            double E0 = Ecfunc(f, U, J, point.mu, theta, E0min, f0min, Ecomps);

            //            out << "Ej5[" << j << "]=" << ::math(Ej) << ";" << endl;
            //            out << "fj5[" << j << "]=" << ::math(f) << ";" << endl;
            //            out << "Ecomps5[" << j << "]=" << ::math(Ecomps) << ";" << endl;

            for (int i = 0; i < L; i++) {
                SiteWaveFunction fi;
                for (int n = 0; n <= nmax; n++) {
                    //                    fi[n] = doublecomplex(uni(rng), uni(rng)); //1;
                    fi[n] = uni(rng);
                }
                fi.normalize();
                fnew[i] = fi;
            }

            failed = false;
            iter = 0;
            theta = parms.theta;
            do {
                for (int i = 0; i < L; i++) {
                    f[i] = fnew[i];
                }

                for (int i = 0; i < L; i++) {
                    createHamiltonian(H, fnew, i, J, U, point.mu, theta);
                    SelfAdjointEigenSolver<Hamiltonian> es(H);
                    if (es.info() != Success) {
                        cout << "Diagonalization failed: " << es.info() << endl;
                        failed = true;
                        break;
                    }
                    fnew[i] = es.eigenvectors().col(0).unaryExpr(std::ptr_fun(chop));
                }
                iter++;
            } while (!allClose(f, fnew, parms.eps) && !failed && iter < MAXITER);

            if (iter == MAXITER) {
                continue;
            }

            for (int i = 0; i < L; i++) {
                f[i] = fnew[i];
            }

            double Eth = Ecfunc(f, U, J, point.mu, theta, Ethmin, fthmin, Ecomps);

            for (int i = 0; i < L; i++) {
                SiteWaveFunction fi;
                for (int n = 0; n <= nmax; n++) {
                    //                    fi[n] = doublecomplex(uni(rng), uni(rng)); //1;
                    fi[n] = uni(rng);
                }
                fi.normalize();
                fnew[i] = fi;
            }

            failed = false;
            iter = 0;
            theta = 2 * parms.theta;
            do {
                for (int i = 0; i < L; i++) {
                    f[i] = fnew[i];
                }

                for (int i = 0; i < L; i++) {
                    createHamiltonian(H, fnew, i, J, U, point.mu, theta);
                    SelfAdjointEigenSolver<Hamiltonian> es(H);
                    if (es.info() != Success) {
                        cout << "Diagonalization failed: " << es.info() << endl;
                        failed = true;
                        break;
                    }
                    fnew[i] = es.eigenvectors().col(0).unaryExpr(std::ptr_fun(chop));
                }
                iter++;
            } while (!allClose(f, fnew, parms.eps) && !failed && iter < MAXITER);

            if (iter == MAXITER) {
                continue;
            }

            for (int i = 0; i < L; i++) {
                f[i] = fnew[i];
            }

            double E2th = Ecfunc(f, U, J, point.mu, theta, E2thmin, f2thmin, Ecomps);

            std::array<double, L> fmax;
            std::array<double, L> fn0;

            for (int i = 0; i < L; i++) {
                fmax[i] = f0min[i].cwiseAbs().maxCoeff();
                fn0[i] = f0min[i].cwiseAbs()[1];
            }

            fres.fmin[point.i][point.j] = *min_element(fn0.begin(), fn0.end());
            fres.fn0[point.i][point.j] = fn0;
            fres.fmax[point.i][point.j] = fmax;
            fres.f0[point.i][point.j] = f0min;
            results.E0res[point.i][point.j] = E0min;

            fres.fth[point.i][point.j] = fthmin;
            results.Ethres[point.i][point.j] = Ethmin;

            fres.f2th[point.i][point.j] = f2thmin;
            results.E2thres[point.i][point.j] = E2thmin;

            results.fs[point.i][point.j] = (E2th - 2 * Eth + E0) / (L * parms.theta * parms.theta);

            {
                boost::mutex::scoped_lock lock(progress_mutex);
                progressbar(thread, j, NRUNS);
                //                ++progress;
            }
        }
        //        cout << "Emin = " << ::math(Emin) << endl;
        //        cout << "fmin = " << ::math(fmin) << endl;

        {
            boost::mutex::scoped_lock progresslock(progress_mutex);
                progressbar(thread, NRUNS, NRUNS);
                
            boost::mutex::scoped_lock pointslock(points_mutex);
//                point++;
                progressbar(numpoints - points.size(), numpoints);
            //                ++progress;
        }
    }

}

#define NSAMP 1

//void phasepoints(Parameter& xi, phase_parameters pparms, queue<Point>& points, /*multi_array<vector<double>, 2 >& f0*/fresults& fres, /*multi_array<double, 2 >& E0res, multi_array<double, 2 >& Ethres, multi_array<double, 2 >& Eth2res, multi_array<double, 2 >& fs,*/results& results, progress_display& progress) {
//
//    mt19937 rng;
//    rng.seed(time(NULL));
//    uniform_real_distribution<> uni(-1, 1);
//
//    int ndim = 2 * L * dim;
//
//    vector<double> x(ndim);
//    doublecomplex * f[L];
//    for (int i = 0; i < L; i++) {
//        f[i] = reinterpret_cast<doublecomplex*> (&x[2 * i * dim]);
//    }
//
//    //    vector<double> U(L), J(L);
//    Parameter U, J;
//
//    vector<double> x0(ndim);
//    doublecomplex * f0[L];
//    for (int i = 0; i < L; i++) {
//        f0[i] = reinterpret_cast<doublecomplex*> (&x0[2 * i * dim]);
//    }
//
//    vector<double> xabs(ndim / 2);
//    double* fabs[L];
//    for (int i = 0; i < L; i++) {
//        fabs[i] = &xabs[i * dim];
//    }
//
//    vector<double> fn0(L);
//    vector<double> fmax(L);
//
//    vector<double> norms(L);
//
//    for (int i = 0; i < L; i++) {
//        //        U[i] = 1 + 0.2 * uni(rng);
//        //        U[i] = 1;
//    }
//
//    funcdata data;
//    data.canonical = pparms.canonical;
//    //        parms.theta = theta;
//
//    double theta = pparms.theta;
//
//    double scale = 1;
//
//    for (;;) {
//        Point point;
//        {
//            boost::mutex::scoped_lock lock(points_mutex);
//            if (points.empty()) {
//                break;
//            }
//            point = points.front();
//            points.pop();
//        }
//        //        cout << "Got queued" << endl;
//
//        //
//        //    vector<double> U(L), J(L);
//        double W[L];
//        for (int i = 0; i < L; i++) {
//            W[i] = xi[i] * point.x;
//        }
//        for (int i = 0; i < L; i++) {
//            //            double Wi = xi[i] * point.x;
//            //            double Wj = xi[i] * point.x;
//
//            U[i] = UW(W[i]) / UW(point.x) / scale;
//            J[i] = JWij(W[i], W[mod(i + 1)]) / UW(point.x) / scale;
//            //            U[i] = 1;
//            //            J[i] = JW(point.x)/UW(point.x);
//
//            ////		U[i] = 0.1 * sqrt(i + 1);
//            ////		J[i] = 0.1 * min(i + 1, mod(i + 1) + 1)
//            ////			+ 0.2 * max(i + 1, mod(i + 1) + 1);
//            //            U[i] = 1+0.2*uni(rng);
//            //                        J[i] = point.x;
//        }
//        //        {
//        //            boost::mutex::scoped_lock lock(points_mutex);
//        //            cout << ::math(U) << endl;
//        //        }
//
//        //    
//        //	parameters parms;
//        data.J = J;
//        data.U = U;
//        data.mu = point.mu / scale;
//        data.Emin = DBL_MAX;
//        data.xmin = vector<double>(ndim);
//
//        //
//        ////    Efuncth(ndim, &x[0], NULL, &parms);
//        ////    return 0;
//        //
//        //        cout << "Setting up optimizer" << endl;
//        nlopt_srand_time();
//        nlopt::opt opt(nlopt::GD_MLSL_LDS, ndim);
//        //                nlopt::opt opt(nlopt::GN_CRS2_LM, ndim);
//        //                nlopt::opt opt(nlopt::GN_CRS2_LM, ndim);
//        //        nlopt::opt opt(nlopt::LD_MMA, ndim);
//        //        nlopt::opt opt(nlopt::LD_LBFGS, ndim);
//        opt.set_lower_bounds(-1);
//        opt.set_upper_bounds(1);
//        opt.set_maxtime(15 * 60);
//        //        opt.set_ftol_rel(1e-5);
//        //        opt.set_xtol_rel(1e-5);
//        //        opt.set_population(1000000);
//        //        opt.set_population(ndim+1);
//        //        opt.set_population(10);
//        //        opt.set_vector_storage(200);
//        //        opt.set_maxtime(10);
//        //                opt.set_xtol_rel(1e-14);
//        //                opt.set_xtol_abs(1e-14);
//        //                opt.set_ftol_rel(1e-14);
//        //                opt.set_ftol_abs(1e-14);
//        //        opt.set_xtol_abs()
//
//        nlopt::opt lbopt(nlopt::LD_LBFGS, ndim);
//        lbopt.set_lower_bounds(-1);
//        lbopt.set_upper_bounds(1);
//        lbopt.set_ftol_rel(1e-6);
//        lbopt.set_xtol_rel(1e-6);
//        //        lbopt.set_vector_storage(10);
//        opt.set_local_optimizer(lbopt);
//
//        nlopt::opt lopt(nlopt::LD_LBFGS, ndim);
//        lopt.set_lower_bounds(-1);
//        lopt.set_upper_bounds(1);
//        //        lopt.set_ftol_rel(1e-2);
//        //        lopt.set_xtol_rel(1e-2);
//        //        opt.set_local_optimizer(lopt);
//
//        opt.set_min_objective(Efunc, &data);
//        lopt.set_min_objective(Efunc, &data);
//        //            cout << "Optimizer set up. Doing optimization" << endl;
//
//        //        int NFC = 240;
//        //        int NSDC = 80;
//        //        double threshold = -1000;
//        //        ProblemSetup* setup = ProblemSetup_create(ndim, NFC, NSDC, threshold);
//        //        setup->userData = &parms;
//        //        setup->costFunction = Ecfuncp;
//        //        setup->randomSeed = time(NULL);
//        //        for(int i = 0; i < ndim; i++) {
//        //            setup->axes[i] = PAxis_create(-1,1);
//        //            setup->axes[i]->axisPrecision = 1e-7;
//        //        }
//        //        
//        //        double minvalue = findMinimum(setup);
//        //        cout << "minvalue = " << minvalue << endl;
//        //        exit(0);
//
//        int res = 0;
//
//        for (int i = 0; i < L; i++) {
//            for (int n = 0; n <= nmax; n++) {
//                //                                f[i][n] = 1 / sqrt(dim);
//                f[i][n] = uni(rng);
//                //                f[i][n] = n == 1 ? 1 : 1e-2;//0;
//                //                f[i][n] = (n == 1 || n == 2) ? 1/sqrt(2) : 0;
//            }
//        }
//
//        vector<double> xmin(ndim);
//        copy(x.begin(), x.end(), xmin.begin());
//
//        data.theta = 0;
//        double E0 = 0;
//        double E0min = DBL_MAX;
//        //        for (int k = 0; k < NSAMP; k++) {
//        //            for (int i = 0; i < L; i++) {
//        //                for (int n = 0; n <= nmax; n++) {
//        //                    f[i][n] = doublecomplex(uni(rng), uni(rng));
//        //                }
//        //            }
//        //            try {
//        //                res = opt.optimize(x, E0);
//        //            } catch (std::exception& e) {
//        //                printf("nlopt failed!: E0: %d, %d\n", point.i, point.j);
//        //                cout << e.what() << endl;
//        //                E0 = numeric_limits<double>::quiet_NaN();
//        //                res = -10;
//        //            }
//        //            if (E0 < E0min) {
//        //                E0min = E0;
//        //                copy(x.begin(), x.end(), xmin.begin());
//        //            }
//        //        }
//        //        cout << "E0 = " << E0 << endl;
//        //        cout << "E0min = " << E0min << endl;
//        //        exit(0);
//
//        try {
//            //            res = lopt.optimize(xmin, E0);
//            res = opt.optimize(x, E0);
//        } catch (std::exception& e) {
//            printf("nlopt failed!: E0 refine: %d, %d\n", point.i, point.j);
//            cout << e.what() << endl;
//            E0 = numeric_limits<double>::quiet_NaN();
//            res = -10;
//        }
//        copy(xmin.begin(), xmin.end(), x.begin());
//        cout << "E0 = " << data.Emin << endl;
//        double E01 = data.Emin;
//
//        try {
//            res = lopt.optimize(data.xmin, E0);
//            //            res = lopt.optimize(x, E0);
//        } catch (std::exception& e) {
//            printf("nlopt failed!: E0 refine: %d, %d\n", point.i, point.j);
//            cout << e.what() << endl;
//            E0 = numeric_limits<double>::quiet_NaN();
//            res = -10;
//        }
//        copy(xmin.begin(), xmin.end(), x.begin());
//        //        cout << "E0 = " << E0 << endl;
//        cout << "Emin = " << E01 << "\tE0 = " << E0 << endl;
//        exit(0);
//
//        norm(x, norms);
//        for (int i = 0; i < L; i++) {
//            for (int n = 0; n <= nmax; n++) {
//                x0[2 * (i * dim + n)] = x[2 * (i * dim + n)] / norms[i];
//                x0[2 * (i * dim + n) + 1] = x[2 * (i * dim + n) + 1] / norms[i];
//            }
//            transform(f0[i], f0[i] + dim, fabs[i], std::ptr_fun<const doublecomplex&, double>(abs));
//            fmax[i] = *max_element(fabs[i], fabs[i] + dim);
//            fn0[i] = fabs[i][1];
//        }
//
//        results.res0[point.i][point.j] = res;
//        fres.fmin[point.i][point.j] = *min_element(fn0.begin(), fn0.end());
//        fres.fn0[point.i][point.j] = fn0;
//        fres.fmax[point.i][point.j] = fmax;
//        fres.f0[point.i][point.j] = x0;
//        results.E0res[point.i][point.j] = E0;
//
//        //        opt.set_min_objective(Ethfunc, &parms);
//
//        for (int i = 0; i < L; i++) {
//            for (int n = 0; n <= nmax; n++) {
//                //                f[i][n] = 1 / sqrt(dim);
//                //                //                f[i][n] = doublecomplex(1/sqrt(dim),1/sqrt(dim));
//                //                f[i][n] = uni(rng);
//                //                f[i][n] = n == 1 ? 1 : 1e-2;//0;
//                //                f[i][n] = (n == 1 || n == 2) ? 1/sqrt(2) : 0;
//            }
//        }
//        data.theta = theta;
//        double Eth = 0;
//        double Ethmin = DBL_MAX;
//        for (int k = 0; k < NSAMP; k++) {
//            for (int i = 0; i < L; i++) {
//                for (int n = 0; n <= nmax; n++) {
//                    f[i][n] = doublecomplex(uni(rng), uni(rng));
//                }
//            }
//            try {
//                res = opt.optimize(x, Eth);
//            } catch (std::exception& e) {
//                printf("nlopt failed!: Eth: %d, %d\n", point.i, point.j);
//                cout << e.what() << endl;
//                Eth = numeric_limits<double>::quiet_NaN();
//                res = -10;
//            }
//            if (Eth < Ethmin) {
//                Ethmin = Eth;
//                copy(x.begin(), x.end(), xmin.begin());
//            }
//        }
//        try {
//            res = lopt.optimize(xmin, Eth);
//        } catch (std::exception& e) {
//            printf("nlopt failed!: Eth refine: %d, %d\n", point.i, point.j);
//            cout << e.what() << endl;
//            Eth = numeric_limits<double>::quiet_NaN();
//            res = -10;
//        }
//        copy(xmin.begin(), xmin.end(), x.begin());
//
//        //        try {
//        //            res = lopt.optimize(x, Eth);
//        //        } catch (std::exception& e) {
//        //            printf("nlopt failed!: Eth: %d, %d\n", point.i, point.j);
//        //            cout << e.what() << endl;
//        //            Eth = numeric_limits<double>::quiet_NaN();
//        //            res = -10;
//        //        }
//
//        //        Encfunc(ndim, x.data(), grad.data(), &parms);
//        //        transform(grad.begin(), grad.end(), grad0.begin(), std::ptr_fun<double,double>(std::abs));
//        //        cout << "th: " << *min_element(grad0.begin(), grad0.end()) << " - " << *max_element(grad0.begin(), grad0.end()) << endl;
//
//        norm(x, norms);
//        for (int i = 0; i < L; i++) {
//            for (int n = 0; n <= nmax; n++) {
//                x0[2 * (i * dim + n)] = x[2 * (i * dim + n)] / norms[i];
//                x0[2 * (i * dim + n) + 1] = x[2 * (i * dim + n) + 1] / norms[i];
//            }
//        }
//
//        results.resth[point.i][point.j] = res;
//        fres.fth[point.i][point.j] = x0;
//        results.Ethres[point.i][point.j] = Eth;
//
//        for (int i = 0; i < L; i++) {
//            for (int n = 0; n <= nmax; n++) {
//                //                f[i][n] = 1 / sqrt(dim);
//                //                f[i][n] = uni(rng);
//                //                f[i][n] = n == 1 ? 1 : 1e-2;//0;
//                //                f[i][n] = (n == 1 || n == 2) ? 1/sqrt(2) : 0;
//            }
//        }
//        data.theta = 2 * theta;
//        double Eth2 = 0;
//        double Eth2min = DBL_MAX;
//        for (int k = 0; k < NSAMP; k++) {
//            for (int i = 0; i < L; i++) {
//                for (int n = 0; n <= nmax; n++) {
//                    f[i][n] = doublecomplex(uni(rng), uni(rng));
//                }
//            }
//            try {
//                res = opt.optimize(x, Eth2);
//            } catch (std::exception& e) {
//                printf("nlopt failed!: Eth2: %d, %d\n", point.i, point.j);
//                cout << e.what() << endl;
//                Eth2 = numeric_limits<double>::quiet_NaN();
//                res = -10;
//            }
//            if (Eth2 < Eth2min) {
//                Eth2min = Eth2;
//                copy(x.begin(), x.end(), xmin.begin());
//            }
//        }
//        try {
//            res = lopt.optimize(xmin, Eth2);
//            //            printf("Twisted energy 2: %0.10g\n", Eth2);
//        } catch (std::exception& e) {
//            printf("nlopt failed!: Eth2 refine: %d, %d\n", point.i, point.j);
//            cout << e.what() << endl;
//            Eth2 = numeric_limits<double>::quiet_NaN();
//            res = -10;
//        }
//        copy(xmin.begin(), xmin.end(), x.begin());
//
//        //        try {
//        //            res = lopt.optimize(x, Eth2);
//        //            //            printf("Twisted energy 2: %0.10g\n", Eth2);
//        //        } catch (std::exception& e) {
//        //            printf("nlopt failed!: Eth2: %d, %d\n", point.i, point.j);
//        //            cout << e.what() << endl;
//        //            Eth2 = numeric_limits<double>::quiet_NaN();
//        //            res = -10;
//        //        }
//
//        //        Encfunc(ndim, x.data(), grad.data(), &parms);
//        //        transform(grad.begin(), grad.end(), grad0.begin(), std::ptr_fun<double,double>(std::abs));
//        //        cout << "th2: " << *min_element(grad0.begin(), grad0.end()) << " - " << *max_element(grad0.begin(), grad0.end()) << endl;
//
//        norm(x, norms);
//        for (int i = 0; i < L; i++) {
//            for (int n = 0; n <= nmax; n++) {
//                x0[2 * (i * dim + n)] = x[2 * (i * dim + n)] / norms[i];
//                x0[2 * (i * dim + n) + 1] = x[2 * (i * dim + n) + 1] / norms[i];
//            }
//        }
//
//        fres.fth2[point.i][point.j] = x0;
//        results.Eth2res[point.i][point.j] = Eth2;
//
//        results.resth2[point.i][point.j] = res;
//        results.fs[point.i][point.j] = (Eth2 - 2 * Eth + E0) / (L * theta * theta);
//        //        cout << "fs = " << (Eth2-2*Eth+E0)/(0.01*0.01) << endl;
//
//        //    
//        //        cout << "Eth - E0 = " << Eth-E0 << endl << endl;
//
//        {
//            boost::mutex::scoped_lock lock(progress_mutex);
//            ++progress;
//        }
//    }
//
//}

vector<double> nu;

/*
 *
 */
int main(int argc, char** argv) {

    //        Hamiltonian H;
    //        SiteWaveFunction psii;
    //        WaveFunction psi;
    //        WaveFunction f;
    //        
    //////        psii.fill(1);
    //////        psii.normalize();
    //////        psi.fill(psii);
    //        double Ww = 280000000000;
    //        Parameter J,U;
    ////        J.fill(0.1);
    ////        U.fill(1);
    //        Parameter W;
    //        Matrix<double, L, 1> xi;
    //        xi << 0.958510999218561,1.24859240406659,1.11016224464402,1.21627868060023,0.750057190540247;
    //        for (int i = 0; i < L; i++) {
    //            W[i] = xi[i] * Ww;
    //        }
    //
    //        for (int i = 0; i < L; i++) {
    //            U[i] = UW(W[i]) / UW(Ww);
    //            J[i] = JWij(W[i], W[mod(i + 1)]) / UW(Ww);
    //        }
    //
    ////        J.fill(JW(Ww)/UW(Ww));
    ////        U.fill(UW(Ww)/UW(Ww));
    //        double theta1 = 0;
    //        cout << "Matrix: " << theta1 << endl;
    //        double Etot = 0;
    //        for(int i = 0; i < L; i++) {
    //        createHamiltonian(H, f, i, J, U, 0, theta1);
    //        doublecomplex E3 = (f[i].adjoint()*H*f[i]).value();
    //        Etot += E3.real();
    //        cout << ::math(E3) << endl;
    //        }
    //        cout << "Etot = " << ::math(Etot) << endl;
    //////        cout << ::math(H) << endl;
    //        return 0;

    //    cout << ::math(cos(0.1)) << endl;
    //    cout << ::math(exp(doublecomplex(0,1)*0.1).real()) << endl;
    //    cout << ::math(exp(-doublecomplex(0,1)*0.1).real()) << endl;

    //                mt19937 rng2;
    //                uniform_real_distribution<> uni2(-1, 1);
    //            
    //                rng2.seed(time(NULL));
    //                nu = vector<double>(L, 0);
    //                for(int i = 0; i < L; i++) {
    //                    nu[i] = 0;//0.5*uni2(rng2);
    //                }
    //            
    //                vector<double> f(2*L*dim,0);
    //                for(int i = 0; i <2*L*dim; i++) {
    //                    f[i] = (i%2==0) ? 1 : 0;//uni2(rng2);
    //                }
    //                    vector<double> g(2*L*dim,0);
    //                    	funcdata parms;
    //                    	parms.J = vector<double>(L,0.1);
    //                    	parms.U = vector<double>(L,1);
    //                    	parms.mu = 0.5;
    //                        parms.theta = 0.2;
    //                    double E1 = Ecfunc(2*L*dim,f.data(),g.data(),&parms);
    //                    cout << "Function: " << parms.theta << endl;
    //                    cout << ::math(E1) << endl;
    //                    return 0;
    //            int id = 2;
    //            for(int id = 2; id < 2*L*dim; id++) {
    //            double df = 1e-7;
    //            f[id] += df;
    //                parms.theta = 0.1;
    //            double E2 = Ecfunc(2*L*dim,f.data(),g.data(),&parms);
    ////            cout << ::math(E2) << endl;
    //            cout << ::math(g[id]) << "\t";//endl;
    //            cout << ::math((E2-E1)/df) << "\t";//endl;
    //            cout << ::math((E2-E1)/df-g[id]) << endl;
    //            f[id] -= df;
    //            }
    //            
    //            return 0;

    mt19937 rng;
    uniform_real_distribution<> uni(-1, 1);

    int seed = lexical_cast<int>(argv[1]);
    int nseed = lexical_cast<int>(argv[2]);

    double xmin = lexical_cast<double>(argv[3]);
    double xmax = lexical_cast<double>(argv[4]);
    int nx = lexical_cast<int>(argv[5]);

    deque<double> x(nx);
    if (nx == 1) {
        x[0] = xmin;
    } else {
        double dx = (xmax - xmin) / (nx - 1);
        for (int ix = 0; ix < nx; ix++) {
            x[ix] = xmin + ix * dx;
        }
    }

    double mumin = lexical_cast<double>(argv[6]);
    double mumax = lexical_cast<double>(argv[7]);
    int nmu = lexical_cast<int>(argv[8]);

    deque<double> mu(nmu);
    if (nmu == 1) {
        mu[0] = mumin;
    } else {
        double dmu = (mumax - mumin) / (nmu - 1);
        for (int imu = 0; imu < nmu; imu++) {
            mu[imu] = mumin + imu * dmu;
        }
    }

    double D = lexical_cast<double>(argv[9]);
    double theta = lexical_cast<double>(argv[10]);

    int numthreads = lexical_cast<int>(argv[11]);

    int resi = lexical_cast<int>(argv[12]);

    //    bool canonical = lexical_cast<bool>(argv[13]);

    double eps = lexical_cast<double>(argv[13]);

#ifdef AMAZON
    //    path resdir("/home/ubuntu/Dropbox/Amazon EC2/Simulation Results/Gutzwiller Phase Diagram");
    path resdir("/home/ubuntu/Dropbox/Amazon EC2/Simulation Results/Canonical Transformation Gutzwiller");
#else
    //    path resdir("/Users/Abuenameh/Dropbox/Amazon EC2/Simulation Results/Gutzwiller Phase Diagram");
    path resdir("/Users/Abuenameh/Documents/Simulation Results/Canonical Transformation Gutzwiller");
#endif
    if (!exists(resdir)) {
        cerr << "Results directory " << resdir << " does not exist!" << endl;
        exit(1);
    }
    for (int iseed = 0; iseed < nseed; iseed++, seed++) {
        ptime begin = microsec_clock::local_time();


        ostringstream oss;
        oss << "res." << resi << ".txt";
        path resfile = resdir / oss.str();
        while (exists(resfile)) {
            resi++;
            oss.str("");
            oss << "res." << resi << ".txt";
            resfile = resdir / oss.str();
        }
        if (seed < 0) {
            resi = seed;
            oss.str("");
            oss << "res." << resi << ".txt";
            resfile = resdir / oss.str();
        }

        Parameter xi;
        xi.fill(1);
        //        xi.assign(1);
        rng.seed(seed);
        if (seed > -1) {
            for (int j = 0; j < L; j++) {
                xi[j] = (1 + D * uni(rng));
            }
        }

        //        rng.seed(seed);
        nu = vector<double>(L, 0);
        //        for (int i = 0; i < L; i++) {
        //            nu[i] = 0.25 * uni(rng);
        //        }

        int Lres = L;
        int nmaxres = nmax;

        boost::filesystem::ofstream os(resfile);
        //        printMath(os, "canonical", resi, canonical);
        printMath(os, "Lres", resi, Lres);
        printMath(os, "nmaxres", resi, nmaxres);
        printMath(os, "seed", resi, seed);
        printMath(os, "theta", resi, theta);
        printMath(os, "Delta", resi, D);
        printMath(os, "xres", resi, x);
        printMath(os, "mures", resi, mu);
        printMath(os, "xires", resi, xi);
        os << flush;

        cout << "Res: " << resi << endl;

        multi_array<double, 2 > fsres(extents[nx][nmu]);
        multi_array<double, 2 > fminres(extents[nx][nmu]);
        multi_array<std::array<double, L>, 2> fn0res(extents[nx][nmu]);
        multi_array<std::array<double, L>, 2> fmaxres(extents[nx][nmu]);
        multi_array<WaveFunction, 2> f0res(extents[nx][nmu]);
        multi_array<WaveFunction, 2> fthres(extents[nx][nmu]);
        multi_array<WaveFunction, 2> f2thres(extents[nx][nmu]);
        multi_array<double, 2> E0res(extents[nx][nmu]);
        multi_array<double, 2> Ethres(extents[nx][nmu]);
        multi_array<double, 2> E2thres(extents[nx][nmu]);

        fresults fres = {fminres, fn0res, fmaxres, f0res, fthres, f2thres};
        results results = {E0res, Ethres, E2thres, fsres};

//        progress_display progress(nx * nmu);

        //        cout << "Queueing" << endl;
        queue<Point> points;
        for (int imu = 0; imu < nmu; imu++) {
            queue<Point> rowpoints;
            for (int ix = 0; ix < nx; ix++) {
                Point point;
                point.i = ix;
                point.j = imu;
                point.x = x[ix];
                point.mu = mu[imu];
                points.push(point);
            }
        }
                progressbar(0, points.size());

        phase_parameters parms;
        parms.theta = theta;
        //        parms.canonical = canonical;
        parms.eps = eps;

        //        cout << "Dispatching" << endl;
        thread_group threads;
        //        vector<thread> threads;
        for (int i = 0; i < numthreads; i++) {
            //                        threads.emplace_back(phasepoints, std::ref(xi), theta, std::ref(points), std::ref(f0res), std::ref(E0res), std::ref(Ethres), std::ref(fsres), std::ref(progress));
            threads.create_thread(bind(&phasepoints2, i, boost::ref(xi), parms, boost::ref(points), boost::ref(fres), boost::ref(results)));
        }
        //        for (thread& t : threads) {
        //            t.join();
        //        }
        threads.join_all();
    printf("\033[%dB", numthreads+2);


        printMath(os, "fsres", resi, fsres);
        printMath(os, "fn0", resi, fn0res);
        printMath(os, "fmin", resi, fminres);
        printMath(os, "fmax", resi, fmaxres);
        printMath(os, "f0res", resi, f0res);
        printMath(os, "fthres", resi, fthres);
        printMath(os, "f2thres", resi, f2thres);
        printMath(os, "E0res", resi, E0res);
        printMath(os, "Ethres", resi, Ethres);
        printMath(os, "E2thres", resi, E2thres);

        ptime end = microsec_clock::local_time();
        time_period period(begin, end);
        cout << endl << period.length() << endl << endl;

        os << "runtime[" << resi << "]=\"" << period.length() << "\";" << endl;
    }

    //    time_t start = time(NULL);
    //
    //    int ndim = 2 * L * dim;

    //    vector<double> x(ndim);
    //
    //    vector<double> U(L), J(L);
    //	for (int i = 0; i < L; i++) {
    ////		U[i] = 0.1 * sqrt(i + 1);
    ////		J[i] = 0.1 * min(i + 1, mod(i + 1) + 1)
    ////			+ 0.2 * max(i + 1, mod(i + 1) + 1);
    //        U[i] = 1;
    //        J[i] = 0.2;
    //	}
    //
    //	doublecomplex * f[L];
    //	for (int i = 0; i < L; i++) {
    //		f[i] = reinterpret_cast<doublecomplex*>(&x[2 * i * dim]);
    //	}
    //
    //	for (int i = 0; i < L; i++) {
    //		for (int n = 0; n <= nmax; n++) {
    //            f[i][n] = 1/sqrt(dim);
    //		}
    //	}
    //    
    //	parameters parms;
    //	parms.J = J;
    //	parms.U = U;
    //	parms.mu = 0.5;
    //    parms.theta = 0.1;
    //
    ////    Efuncth(ndim, &x[0], NULL, &parms);
    ////    return 0;
    //
    //    nlopt::opt opt(/*nlopt::GN_ISRES*/nlopt::LN_COBYLA/*nlopt::LD_SLSQP*/, ndim);
    ////    nlopt::opt opt(nlopt::AUGLAG/*nlopt::GN_ISRES*//*nlopt::LN_COBYLA*//*nlopt::LD_SLSQP*/, ndim);
    ////    nlopt::opt local_opt(nlopt::LN_SBPLX, ndim);
    ////    opt.set_local_optimizer(local_opt);
    //    opt.set_lower_bounds(-1);
    //    opt.set_upper_bounds(1);
    //    vector<double> ctol(L, 1e-8);
    //    opt.add_equality_mconstraint(norm2s, NULL, ctol);
    //    opt.set_min_objective(Efunc, &parms);
    //    opt.set_xtol_rel(1e-8);
    //
    //	int res = 0;
    //    
    //    double E0 = 0;
    //    try {
    //        res = opt.optimize(x, E0);
    //        printf("Found minimum: %0.10g\n", E0);
    //    }
    //    catch(exception& e) {
    //        printf("nlopt failed! %d\n", res);
    //        cout << e.what() << endl;
    //    }
    //    
    //    opt.set_min_objective(Efuncth, &parms);
    //    
    //    double Eth = 0;
    //    try {
    //        res = opt.optimize(x, Eth);
    //        printf("Found minimum: %0.10g\n", Eth);
    //    }
    //    catch(exception& e) {
    //        printf("nlopt failed! %d\n", res);
    //        cout << e.what() << endl;
    //    }
    //    
    //    cout << "Eth - E0 = " << Eth-E0 << endl << endl;
    //
    //    for(int i = 0; i < 1; i++) {
    //        for(int n = 0; n <= nmax; n++) {
    //        cout << norm(f[i][n]) << endl;
    //        }
    //        cout << endl;
    //    }

    //    vector<double> norm2is;
    //	cout << norm2(x, norm2is) << endl;
    //	cout << E0 / norm2(x, norm2is) << endl;

    //    nlopt_set_xtol_rel(opt, 1e-8);
    //
    //    double x[2] = {1.234, 5.678}; /* some initial guess */
    //    double minf; /* the minimum objective value, upon return */
    //
    //    int res = 0;
    //    if ((res = nlopt_optimize(opt, x, &minf)) < 0) {
    //        printf("nlopt failed! %d\n",res);
    //    } else {
    //        printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
    //    }
    //
    //nlopt_destroy(opt);


    //    time_t end = time(NULL);
    //
    //    printf("Runtime: %ld", end - start);

    return 0;
}

