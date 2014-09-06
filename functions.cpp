#include "gutzwiller.hpp"

#include <iomanip>

//double Encfunc(unsigned ndim, WaveFunction& f, Parameter& U, Parameter& J, double mu, double theta, EnergyComponents& Ei) {
////    funcdata* parms = static_cast<funcdata*> (data);
////    vector<double> U = parms->U;
////    Parameter U = parms->U;
////    Parameter J = parms->J;
////    double mu = parms->mu;
//////    vector<double> J = parms->J;
////    double theta = parms->theta;
//
//    doublecomplex expth = exp(doublecomplex(0, 1) * theta);
//    doublecomplex expmth = exp(-doublecomplex(0, 1) * theta);
//
//    doublecomplex Ec = 0;
//
//    vector<doublecomplex> E0s(L, 0), E1j1s(L, 0), E1j2s(L, 0);
//
//
//    for (int i = 0; i < L; i++) {
//
//        int j1 = mod(i - 1);
//        int j2 = mod(i + 1);
//
//        doublecomplex E0 = 0;
//        doublecomplex E1j1 = 0;
//        doublecomplex E1j2 = 0;
//
//        for (int n = 0; n <= nmax; n++) {
//
//            E0 += (0.5 * U[i] * n * (n - 1) - mu * n + nu[i] * n) * ~f[i][n] * f[i][n];
//
//            if (n < nmax) {
//                for (int m = 1; m <= nmax; m++) {
//                    E1j1 += -J[j1] * expth * g(n, m) * ~f[i][n + 1] * ~f[j1][m - 1]
//                            * f[i][n] * f[j1][m];
//                    E1j2 += -J[i] * expmth * g(n, m) * ~f[i][n + 1] * ~f[j2][m - 1]
//                            * f[i][n] * f[j2][m];
//                }
//            }
//
//        }
//
//        Ec += E0;
//        Ec += E1j1;
//        Ec += E1j2;
//
//        E0s[i] += E0;
//        E1j1s[i] += E1j1;
//        E1j2s[i] += E1j2;
//
//    }
//
//    //    cout << Ec.real() << endl;
//    return Ec.real();
//}

double Ecfunc(WaveFunction& f, Parameter& U, Parameter& J, double mu, double theta, double& Emin, WaveFunction& fmin, EnergyComponents& Ecomps) {
    //    funcdata* fdata = static_cast<funcdata*> (data);
    //    Parameter U = fdata->U;
    //    Parameter J = fdata->J;
    ////    vector<double> U = fdata->U;
    //    double mu = fdata->mu;
    ////    vector<double> J = fdata->J;
    //    double theta = fdata->theta;
    
//    EnergyComponents Ecomps;

    doublecomplex expth = exp(doublecomplex(0, 1) * theta);
    doublecomplex expmth = exp(-doublecomplex(0, 1) * theta);
    doublecomplex exp2th = exp(doublecomplex(0, 1)*2.0 * theta);
    doublecomplex expm2th = exp(-doublecomplex(0, 1)*2.0 * theta);

    doublecomplex Ec = 0;

    //    vector<doublecomplex> E0s(L, 0), E1j1s(L, 0), E1j2s(L, 0),
    //            E2j1s(L, 0), E2j2s(L, 0), E3j1s(L, 0), E3j2s(L, 0),
    //            E4j1j2s(L, 0), E4j1k1s(L, 0), E4j2k2s(L, 0),
    //            E5j1j2s(L, 0), E5j1k1s(L, 0), E5j2k2s(L, 0);




    for (int i = 0; i < L; i++) {

        int k1 = mod(i - 2);
        int j1 = mod(i - 1);
        int j2 = mod(i + 1);
        int k2 = mod(i + 2);

        doublecomplex E0 = 0;
        doublecomplex E1j1 = 0;
        doublecomplex E1j2 = 0;
        doublecomplex E2j1 = 0;
        doublecomplex E2j2 = 0;
        doublecomplex E3j1 = 0;
        doublecomplex E3j2 = 0;
        doublecomplex E4j1j2 = 0;
        doublecomplex E4j1k1 = 0;
        doublecomplex E4j2k2 = 0;
        doublecomplex E5j1j2 = 0;
        doublecomplex E5j1k1 = 0;
        doublecomplex E5j2k2 = 0;

        for (int n = 0; n <= nmax; n++) {
            E0 += (0.5 * U[i] * n * (n - 1) - mu * n + nu[i] * n) * ~f[i][n] * f[i][n];

            if (n < nmax) {
                E1j1 += -J[j1] * expth * g(n, n + 1) * ~f[i][n + 1] * ~f[j1][n]
                        * f[i][n] * f[j1][n + 1];
                E1j2 += -J[i] * expmth * g(n, n + 1) * ~f[i][n + 1] * ~f[j2][n] * f[i][n]
                        * f[j2][n + 1];

                if (n > 0) {
                    E2j1 += 0.5 * J[j1] * J[j1] * exp2th * g(n, n) * g(n - 1, n + 1)
                            * ~f[i][n + 1] * ~f[j1][n - 1] * f[i][n - 1] * f[j1][n + 1]
                            * (1 / eps(U, i, j1, n, n) - 1 / eps(U, i, j1, n - 1, n + 1));
                    E2j2 += 0.5 * J[i] * J[i] * expm2th * g(n, n) * g(n - 1, n + 1)
                            * ~f[i][n + 1] * ~f[j2][n - 1] * f[i][n - 1] * f[j2][n + 1]
                            * (1 / eps(U, i, j2, n, n) - 1 / eps(U, i, j2, n - 1, n + 1));
                }

                for (int m = 1; m <= nmax; m++) {
                    if (n != m - 1) {
                        E3j1 += 0.5 * (J[j1] * J[j1] / eps(U, i, j1, n, m)) * g(n, m)
                                * g(m - 1, n + 1)
                                * (~f[i][n + 1] * ~f[j1][m - 1] * f[i][n + 1] * f[j1][m - 1]
                                - ~f[i][n] * ~f[j1][m] * f[i][n] * f[j1][m]);
                        E3j2 += 0.5 * (J[i] * J[i] / eps(U, i, j2, n, m)) * g(n, m)
                                * g(m - 1, n + 1)
                                * (~f[i][n + 1] * ~f[j2][m - 1] * f[i][n + 1] * f[j2][m - 1]
                                - ~f[i][n] * ~f[j2][m] * f[i][n] * f[j2][m]);
                    }
                }

                if (n > 0) {
                    E4j1j2 += 0.5 * (J[j1] * J[i] / eps(U, i, j1, n, n)) * g(n, n)
                            * g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j1][n - 1] * ~f[j2][n]
                            * f[i][n - 1] * f[j1][n] * f[j2][n + 1];
                    E4j1j2 += 0.5 * (J[i] * J[j1] / eps(U, i, j2, n, n)) * g(n, n)
                            * g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j2][n - 1] * ~f[j1][n]
                            * f[i][n - 1] * f[j2][n] * f[j1][n + 1];
                    E4j1k1 += 0.5 * (J[j1] * J[k1] / eps(U, i, j1, n, n)) * g(n, n)
                            * g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j1][n - 1] * ~f[k1][n]
                            * f[i][n] * f[j1][n + 1] * f[k1][n - 1];
                    E4j2k2 += 0.5 * (J[i] * J[j2] / eps(U, i, j2, n, n)) * g(n, n)
                            * g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j2][n - 1] * ~f[k2][n]
                            * f[i][n] * f[j2][n + 1] * f[k2][n - 1];
                    E4j1j2 -= 0.5 * (J[j1] * J[i] / eps(U, i, j1, n - 1, n + 1))
                            * g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j1][n]
                            * ~f[j2][n - 1] * f[i][n - 1] * f[j1][n + 1] * f[j2][n];
                    E4j1j2 -= 0.5 * (J[i] * J[j1] / eps(U, i, j2, n - 1, n + 1))
                            * g(n, n) * g(n - 1, n + 1) * ~f[i][n + 1] * ~f[j2][n]
                            * ~f[j1][n - 1] * f[i][n - 1] * f[j2][n + 1] * f[j1][n];
                    E4j1k1 -= 0.5 * (J[j1] * J[k1] / eps(U, i, j1, n - 1, n + 1))
                            * g(n, n) * g(n - 1, n + 1) * ~f[i][n] * ~f[j1][n - 1]
                            * ~f[k1][n + 1] * f[i][n - 1] * f[j1][n + 1] * f[k1][n];
                    E4j2k2 -= 0.5 * (J[i] * J[j2] / eps(U, i, j2, n - 1, n + 1))
                            * g(n, n) * g(n - 1, n + 1) * ~f[i][n] * ~f[j2][n - 1]
                            * ~f[k2][n + 1] * f[i][n - 1] * f[j2][n + 1] * f[k2][n];
                }

                for (int m = 1; m <= nmax; m++) {
                    if (n != m - 1 && n < nmax) {
                        E5j1j2 += 0.5 * (J[j1] * J[i] * exp2th / eps(U, i, j1, n, m))
                                * g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1] * ~f[j1][m - 1]
                                * ~f[j2][m] * f[i][n + 1] * f[j1][m] * f[j2][m - 1];
                        E5j1j2 += 0.5 * (J[i] * J[j1] * expm2th / eps(U, i, j2, n, m))
                                * g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1] * ~f[j2][m - 1]
                                * ~f[j1][m] * f[i][n + 1] * f[j2][m] * f[j1][m - 1];
                        E5j1k1 += 0.5 * (J[j1] * J[k1] * exp2th / eps(U, i, j1, n, m))
                                * g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1] * ~f[j1][m - 1]
                                * ~f[k1][n] * f[i][n] * f[j1][m - 1] * f[k1][n + 1];
                        E5j2k2 += 0.5 * (J[i] * J[j2] * expm2th / eps(U, i, j2, n, m))
                                * g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1] * ~f[j2][m - 1]
                                * ~f[k2][n] * f[i][n] * f[j2][m - 1] * f[k2][n + 1];
                        E5j1j2 -= 0.5 * (J[j1] * J[i] * exp2th / eps(U, i, j1, n, m))
                                * g(n, m) * g(m - 1, n + 1) * ~f[i][n] * ~f[j1][m - 1]
                                * ~f[j2][m] * f[i][n] * f[j1][m] * f[j2][m - 1];
                        E5j1j2 -= 0.5 * (J[i] * J[j1] * expm2th / eps(U, i, j2, n, m))
                                * g(n, m) * g(m - 1, n + 1) * ~f[i][n] * ~f[j2][m - 1]
                                * ~f[j1][m] * f[i][n] * f[j2][m] * f[j1][m - 1];
                        E5j1k1 -= 0.5 * (J[j1] * J[k1] * exp2th / eps(U, i, j1, n, m))
                                * g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1] * ~f[j1][m]
                                * ~f[k1][n] * f[i][n] * f[j1][m] * f[k1][n + 1];
                        E5j2k2 -= 0.5 * (J[i] * J[j2] * expm2th / eps(U, i, j2, n, m))
                                * g(n, m) * g(m - 1, n + 1) * ~f[i][n + 1] * ~f[j2][m]
                                * ~f[k2][n] * f[i][n] * f[j2][m] * f[k2][n + 1];
                    }
                }
            }

        }

        Ec += E0;

        Ec += E1j1;
        Ec += E1j2;

        Ec += E2j1;
        Ec += E2j2;

        Ec += E3j1;
        Ec += E3j2;

        Ec += E4j1j2;
        Ec += E4j1k1;
        Ec += E4j2k2;

        Ec += E5j1j2;
        Ec += E5j1k1;
        Ec += E5j2k2;

        Ecomps[i][0] = E0.real();
        Ecomps[i][1] = (E1j1 + E1j2).real();
        Ecomps[i][2] = (E2j1 + E2j2).real();
        Ecomps[i][3] = (E3j1 + E3j2).real();
        Ecomps[i][4] = (E4j1j2 + E4j1k1 + E4j2k2).real();
        Ecomps[i][5] = (E5j1j2 + E5j1k1 + E5j2k2).real();
    }

    if (Ec.real() < Emin) {
        Emin = Ec.real();
        fmin = f;
//        Ecompso = Ecomps;
//        copy(x, x + ndim, fdata->xmin.begin());
    }

    //    static double Emin = 1e10;
    //    Emin = min(Emin, Ec.real());
    ////    cout << setprecision(10) << Emin << endl;
    //    cout << setprecision(10) << Ec.real() << "\t" << Emin << endl;

    //    cout << setprecision(10) << Ec.real() << endl;
    return Ec.real();

}

