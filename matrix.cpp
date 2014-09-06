#include "gutzwiller.hpp"

void createHamiltonian(Hamiltonian& Hout, WaveFunction& f, int i, Parameter& J, Parameter& U, double mu, double theta) {
    doublecomplex expth = exp(doublecomplex(0, 1) * theta);
    doublecomplex expmth = exp(-doublecomplex(0, 1) * theta);
    doublecomplex exp2th = exp(doublecomplex(0, 1)*2.0 * theta);
    doublecomplex expm2th = exp(-doublecomplex(0, 1)*2.0 * theta);

    Hamiltonian H = Hamiltonian::Zero();

    int k1 = mod(i - 2);
    int j1 = mod(i - 1);
    int j2 = mod(i + 1);
    int k2 = mod(i + 2);

    for (int n = 0; n <= nmax; n++) {
        H(n, n) += 0.5 * U[i] * n * (n - 1) - mu * n/* + nu[i] * n*/;

        if (n < nmax) {
            H(n + 1, n) += -J[j1] * expth * g(n, n + 1) * ~f[j1][n] * f[j1][n + 1];
            H(n, n + 1) += -J[j1] * expmth * g(n, n + 1) * ~f[j1][n + 1] * f[j1][n];
            H(n + 1, n) += -J[i] * expmth * g(n, n + 1) * ~f[j2][n] * f[j2][n + 1];
            H(n, n + 1) += -J[i] * expth * g(n, n + 1) * ~f[j2][n + 1] * f[j2][n];

            if (n > 0) {
                H(n + 1, n - 1) += 0.5 * J[j1] * J[j1] * exp2th * g(n, n) * g(n - 1, n + 1)
                        * ~f[j1][n - 1] * f[j1][n + 1]
                        * (1 / eps(U, i, j1, n, n) - 1 / eps(U, i, j1, n - 1, n + 1));
                H(n + 1, n - 1) += 0.5 * J[i] * J[i] * expm2th * g(n, n) * g(n - 1, n + 1)
                        * ~f[j2][n - 1] * f[j2][n + 1]
                        * (1 / eps(U, i, j2, n, n) - 1 / eps(U, i, j2, n - 1, n + 1));
                H(n - 1, n + 1) += 0.5 * J[j1] * J[j1] * expm2th * g(n, n) * g(n - 1, n + 1)
                        * ~f[j1][n + 1] * f[j1][n - 1]
                        * (1 / eps(U, i, j1, n, n) - 1 / eps(U, i, j1, n - 1, n + 1));
                H(n - 1, n + 1) += 0.5 * J[i] * J[i] * exp2th * g(n, n) * g(n - 1, n + 1)
                        * ~f[j2][n + 1] * f[j2][n - 1]
                        * (1 / eps(U, i, j2, n, n) - 1 / eps(U, i, j2, n - 1, n + 1));
            }

            for (int m = 1; m <= nmax; m++) {
                if (n != m - 1) {
                    H(n + 1, n + 1) += 2*0.5 * (J[j1] * J[j1] / eps(U, i, j1, n, m)) * g(n, m)
                            * g(m - 1, n + 1)
                            * ~f[j1][m - 1] * f[j1][m - 1];
                    H(n, n) -= 2*0.5 * (J[j1] * J[j1] / eps(U, i, j1, n, m)) * g(n, m)
                            * g(m - 1, n + 1)
                            * ~f[j1][m] * f[j1][m];
                    H(n + 1, n + 1) += 2*0.5 * (J[i] * J[i] / eps(U, i, j2, n, m)) * g(n, m)
                            * g(m - 1, n + 1)
                            * ~f[j2][m - 1] * f[j2][m - 1];
                    H(n, n) -= 2*0.5 * (J[i] * J[i] / eps(U, i, j2, n, m)) * g(n, m)
                            * g(m - 1, n + 1)
                            * ~f[j2][m] * f[j2][m];
                }
            }

            if (n > 0) {
                H(n + 1, n - 1) += 0.5 * (J[j1] * J[i] / eps(U, i, j1, n, n)) * g(n, n)
                        * g(n - 1, n + 1) * ~f[j1][n - 1] * ~f[j2][n]
                        * f[j1][n] * f[j2][n + 1];
                H(n - 1, n + 1) -= 0.5 * (J[j1] * J[i] / eps(U, j1, i, n - 1, n + 1))
                        * g(n, n) * g(n - 1, n + 1) * ~f[j1][n]
                        * ~f[j2][n + 1] * f[j1][n - 1] * f[j2][n];

                H(n + 1, n - 1) += 0.5 * (J[i] * J[j1] / eps(U, i, j2, n, n)) * g(n, n)
                        * g(n - 1, n + 1) * ~f[j2][n - 1] * ~f[j1][n]
                        * f[j2][n] * f[j1][n + 1];
                H(n - 1, n + 1) -= 0.5 * (J[i] * J[j1] / eps(U, j2, i, n - 1, n + 1))
                        * g(n, n) * g(n - 1, n + 1) * ~f[j2][n]
                        * ~f[j1][n + 1] * f[j2][n - 1] * f[j1][n];

                H(n + 1, n) += 0.5 * (J[j1] * J[k1] / eps(U, i, j1, n, n)) * g(n, n)
                        * g(n - 1, n + 1) * ~f[j1][n - 1] * ~f[k1][n]
                        * f[j1][n + 1] * f[k1][n - 1];
                H(n, n + 1) -= 0.5 * (J[j1] * J[k1] / eps(U, j1, i, n - 1, n + 1))
                        * g(n, n) * g(n - 1, n + 1) * ~f[j1][n + 1]
                        * ~f[k1][n - 1] * f[j1][n - 1] * f[k1][n];

                H(n + 1, n) += 0.5 * (J[i] * J[j2] / eps(U, i, j2, n, n)) * g(n, n)
                        * g(n - 1, n + 1) * ~f[j2][n - 1] * ~f[k2][n]
                        * f[j2][n + 1] * f[k2][n - 1];
                H(n, n + 1) -= 0.5 * (J[i] * J[j2] / eps(U, j2, i, n - 1, n + 1))
                        * g(n, n) * g(n - 1, n + 1) * ~f[j2][n + 1]
                        * ~f[k2][n - 1] * f[j2][n - 1] * f[k2][n];

                H(n + 1, n - 1) -= 0.5 * (J[j1] * J[i] / eps(U, i, j1, n - 1, n + 1))
                        * g(n, n) * g(n - 1, n + 1) * ~f[j1][n]
                        * ~f[j2][n - 1] * f[j1][n + 1] * f[j2][n];
                H(n - 1, n + 1) += 0.5 * (J[j1] * J[i] / eps(U, j1, i, n, n)) * g(n, n)
                        * g(n - 1, n + 1) * ~f[j1][n + 1] * ~f[j2][n]
                        * f[j1][n] * f[j2][n - 1];

                H(n + 1, n - 1) -= 0.5 * (J[i] * J[j1] / eps(U, i, j2, n - 1, n + 1))
                        * g(n, n) * g(n - 1, n + 1) * ~f[j2][n]
                        * ~f[j1][n - 1] * f[j2][n + 1] * f[j1][n];
                H(n - 1, n + 1) += 0.5 * (J[i] * J[j1] / eps(U, j2, i, n, n)) * g(n, n)
                        * g(n - 1, n + 1) * ~f[j2][n + 1] * ~f[j1][n]
                        * f[j2][n] * f[j1][n - 1];

                H(n, n - 1) -= 0.5 * (J[j1] * J[k1] / eps(U, i, j1, n - 1, n + 1))
                        * g(n, n) * g(n - 1, n + 1) * ~f[j1][n - 1]
                        * ~f[k1][n + 1] * f[j1][n + 1] * f[k1][n];
                H(n - 1, n) += 0.5 * (J[j1] * J[k1] / eps(U, j1, i, n, n)) * g(n, n)
                        * g(n - 1, n + 1) * ~f[j1][n + 1] * ~f[k1][n]
                        * f[j1][n - 1] * f[k1][n + 1];

                H(n, n - 1) -= 0.5 * (J[i] * J[j2] / eps(U, i, j2, n - 1, n + 1))
                        * g(n, n) * g(n - 1, n + 1) * ~f[j2][n - 1]
                        * ~f[k2][n + 1] * f[j2][n + 1] * f[k2][n];
                H(n - 1, n) += 0.5 * (J[i] * J[j2] / eps(U, j2, i, n, n)) * g(n, n)
                        * g(n - 1, n + 1) * ~f[j2][n + 1] * ~f[k2][n]
                        * f[j2][n - 1] * f[k2][n + 1];

            }

            for (int m = 1; m <= nmax; m++) {
                if (n != m - 1 && n < nmax) {

                    H(n + 1, n + 1) += 0.5 * (J[j1] * J[i] * exp2th / eps(U, i, j1, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j1][m - 1]
                            * ~f[j2][m] * f[j1][m] * f[j2][m - 1];
                    H(n + 1, n + 1) += 0.5 * (J[i] * J[j1] * expm2th / eps(U, i, j2, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j2][m - 1]
                            * ~f[j1][m] * f[j2][m] * f[j1][m - 1];

                    H(n + 1, n + 1) += 0.5 * (J[j1] * J[i] * exp2th / eps(U, i, j2, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j1][m - 1]
                            * ~f[j2][m] * f[j1][m] * f[j2][m - 1];
                    H(n + 1, n + 1) += 0.5 * (J[i] * J[j1] * expm2th / eps(U, i, j1, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j2][m - 1]
                            * ~f[j1][m] * f[j2][m] * f[j1][m - 1];


                    H(n + 1, n) += 0.5 * (J[j1] * J[k1] * exp2th / eps(U, i, j1, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j1][m - 1]
                            * ~f[k1][n] * f[j1][m - 1] * f[k1][n + 1];
                    H(n, n + 1) += 0.5 * (J[k1] * J[j1] * expm2th / eps(U, i, j1, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[k1][n + 1] * ~f[j1][m - 1]
                            * f[k1][n] * f[j1][m - 1];
                    H(n + 1, n) += 0.5 * (J[i] * J[j2] * expm2th / eps(U, i, j2, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j2][m - 1]
                            * ~f[k2][n] * f[j2][m - 1] * f[k2][n + 1];
                    H(n, n + 1) += 0.5 * (J[j2] * J[i] * exp2th / eps(U, i, j2, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[k2][n + 1] * ~f[j2][m - 1]
                            * f[k2][n] * f[j2][m - 1];

                    H(n, n) -= 0.5 * (J[j1] * J[i] * exp2th / eps(U, i, j1, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j1][m - 1]
                            * ~f[j2][m] * f[j1][m] * f[j2][m - 1];
                    H(n, n) -= 0.5 * (J[i] * J[j1] * expm2th / eps(U, i, j2, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j2][m - 1]
                            * ~f[j1][m] * f[j2][m] * f[j1][m - 1];

                    H(n, n) -= 0.5 * (J[j1] * J[i] * exp2th / eps(U, i, j2, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j1][m - 1]
                            * ~f[j2][m] * f[j1][m] * f[j2][m - 1];
                    H(n, n) -= 0.5 * (J[i] * J[j1] * expm2th / eps(U, i, j1, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j2][m - 1]
                            * ~f[j1][m] * f[j2][m] * f[j1][m - 1];

                    H(n + 1, n) -= 0.5 * (J[j1] * J[k1] * exp2th / eps(U, i, j1, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j1][m]
                            * ~f[k1][n] * f[j1][m] * f[k1][n + 1];
                    H(n, n + 1) -= 0.5 * (J[k1] * J[j1] * expm2th / eps(U, i, j1, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[k1][n + 1] * ~f[j1][m]
                            * f[k1][n] * f[j1][m];
                    H(n + 1, n) -= 0.5 * (J[i] * J[j2] * expm2th / eps(U, i, j2, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[j2][m]
                            * ~f[k2][n] * f[j2][m] * f[k2][n + 1];
                    H(n, n + 1) -= 0.5 * (J[j2] * J[i] * exp2th / eps(U, i, j2, n, m))
                            * g(n, m) * g(m - 1, n + 1) * ~f[k2][n + 1] * ~f[j2][m]
                            * f[k2][n] * f[j2][m];
                }
            }
        }
    }

    //    Hout = 0.5 * (H + H.adjoint());
    Hout = H;
}
