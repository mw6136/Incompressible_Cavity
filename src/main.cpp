#include <iomanip>

#include <iostream>

#include "grid.hpp"

#include "omp.h"

void
ProblemGenerator(SimpleMesh * pmesh, double p0, double T0) {

    double u0 = 0.0;
    double v0 = 0.0;

    double rho0 = p0 / (pmesh -> R * T0);
    double e0 = 1 / (pmesh -> gamma - 1) * pmesh -> R * T0;

    for (int i = 0; i < pmesh -> m_nx1 + 2 * pmesh -> m_ng; ++i) {
        for (int j = 0; j < pmesh -> m_nx2 + 2 * pmesh -> m_ng; ++j) {

            pmesh -> un(i, j) = u0;
            pmesh -> vn(i, j) = v0;
            pmesh -> p(i, j) = p0;
        }
    }
}

int
main(int argc, char * argv[]) {

    int nx1 = 64;
    int nx2 = 64;
    int ng = 1;

    double gamma = 1.4;
    double R = 287;

    double CFL = 0.5;
    int printfreq = 100;

    double Re = 100;
    double T0 = 300; // K
    double Ma = 0.8; // Mach number
    double P0 = 101330;

    double L = 1.0 / (Ma / 0.025);
    double x1max = L;
    double x1min = 0.0;
    double x2max = L;
    double x2min = 0.0;

    double a = std::sqrt(gamma * R * T0);
    double Uw = a * Ma;
    double nu = Uw * L / Re;
    double omega = 1 / (L * L) * 2 * nu;

    double rho = P0 / (R * T0);

    double omegat_max = 10.0;
    double savedt = 1.0; // in omega ts

    double tol = 1e-5;

    omp_set_num_threads(16);

    std::string time_int;
    time_int = "euler";

    std::string savepref = "test.output.";

    // --------------------------------------------------------------------------------------------//

    SimpleMesh * pmesh;
    pmesh = new SimpleMesh(nx1, nx2, x1max, x1min, x2max, x2min, ng, gamma, R, "adiabatic");

    ProblemGenerator(pmesh, 101330, 300);

    // pmesh->saveOutput(69, 0.0, savepref);

    double t = 0;
    int saveiter = 0;
    double dt = CFL * pmesh -> m_dx / (Uw);
    int iter = 0;
    int savefreq = (int)(savedt / (omega * dt));

    while (omega * t <= omegat_max) {

        if (omega * t > omegat_max) {
            break;
        }

        if (iter % printfreq == 0) {
            std::cout << "Iteration : " << iter << ", Time (omega * t): " << omega * t <<
                ", dt = " << dt << std::endl;
        }

        pmesh -> new_enforceBCs(t, Uw, omega, nu);

        pmesh -> fillHs(nu);

        pmesh -> pressureCorrect(iter, rho, tol, t, Uw, omega, nu);

        pmesh -> new_enforceBCs(t, Uw, omega, nu);

        if (time_int == "euler") {

            double pipoh, pimoh, pjpoh, pjmoh;

            for (int i = 1; i <= pmesh -> m_nx1; ++i) {
                for (int j = 1; j <= pmesh -> m_nx2; ++j) {

                    pipoh = (pmesh -> p(i + 1, j) + pmesh -> p(i, j)) / 2;
                    pimoh = (pmesh -> p(i, j) + pmesh -> p(i - 1, j)) / 2;
                    pjpoh = (pmesh -> p(i, j + 1) + pmesh -> p(i, j)) / 2;
                    pjmoh = (pmesh -> p(i, j) + pmesh -> p(i, j - 1)) / 2;

                    pmesh -> un(i, j) += dt * pmesh -> Hx(i, j) - dt / rho * ((pipoh - pimoh) / pmesh -> m_dx);
                    pmesh -> vn(i, j) += dt * pmesh -> Hy(i, j) - dt / rho * ((pjpoh - pjmoh) / pmesh -> m_dy);

                    // if (i == 1 && j == 50) {
                    //     std::cout << "dpdx Term (i,j) = (1,50): " << (pipoh - pimoh)/pmesh->m_dy <<
                    //     std::endl;

                    // } else if (i == 50 && j == 1) {
                    //     std::cout << "dpdy Term (i,j) = (50,1): " <<  (pjpoh - pjmoh)/pmesh->m_dy
                    //     << std::endl;
                    // }
                }
            }
        }

        pmesh -> new_enforceBCs(t, Uw, omega, nu);

        // pmesh->swapArrays();

        if (iter % savefreq == 0) {
            std::cout << "Iteration : " << iter << ", Time (omega * t): " << omega * t <<
                ", dt = " << dt << ", Data saved" << std::endl;
            pmesh -> saveOutput(saveiter, t, savepref);
            saveiter += 1;
        }

        t += dt;
        iter += 1;
    }

    return 0;
}