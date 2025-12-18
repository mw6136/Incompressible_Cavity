#include <cmath>

#include <cstdlib>

#include <fstream>

#include <sstream>

#include <string>

#include <vector>

#include "H5Cpp.h"

#ifndef GRID_HPP_
#define GRID_HPP_

std::string
ZeroPadNumber(int num) {
    std::ostringstream ss;
    ss << std::setw(4) << std::setfill('0') << num;
    std::string result = ss.str();
    if (result.length() > 7) {
        result.erase(0, result.length() - 7);
    }
    return result;
}

std::string
getFname(int iter, std::string prefix) {

    std::string fname;
    std::stringstream ss;
    std::stringstream ssfname;
    std::string newnum = ZeroPadNumber(iter);

    ss << prefix << newnum << ".h5";

    ss >> fname;

    return fname;
}

class Array2D {
    int x1s_, x2s_;

    public:
        std::vector < double > m_data;

    Array2D(size_t x1, size_t x2, double init = 0.0): x1s_(x1), x2s_(x2), m_data(x1s_ * x2s_, init) {}
    double &
        operator()(size_t x, size_t y) {
            return m_data[y + x1s_ * (x)];
        }
    double
    operator()(size_t x, size_t y) const {
        return m_data[y + x1s_ * (x)];
    }
};

class SimpleMesh {

    public: std::string m_wall_cond;

    // Number of grid points in x1, x2 direction and number of ghost zones
    int m_nx1,
    m_nx2,
    m_ng;
    double m_x1max,
    m_x1min,
    m_x2max,
    m_x2min,
    m_dx,
    m_dy;

    // Fluid quantities to track
    double gamma,
    R;

    // Vector for Un, Unp1, x direction flux (F), y direction flux (F)
    Array2D un,
    vn,
    unp1,
    vnp1,
    Hx,
    Hy,
    p,
    pnew,
    dpdx,
    dpdy;

    // Arrays for other important intermediate quantities

    // Arrays to store the derivative information for viscosity, temp diffusion
    // at the locations of where the fluxes are

    SimpleMesh(int nx1, int nx2, double x1max, double x1min, double x2max, double x2min, int ng,
        double i_gamma, double i_R, std::string wall_cond): un(nx1 + 2 * ng, nx2 + 2 * ng),
    vn(nx1 + 2 * ng, nx2 + 2 * ng),
    unp1(nx1 + 2 * ng, nx2 + 2 * ng),
    vnp1(nx1 + 2 * ng, nx2 + 2 * ng),
    Hx(nx1 + 2 * ng, nx2 + 2 * ng),
    Hy(nx1 + 2 * ng, nx2 + 2 * ng),
    p(nx1 + 2 * ng, nx2 + 2 * ng),
    dpdx(nx1 + 2 * ng, nx2 + 2 * ng),
    dpdy(nx1 + 2 * ng, nx2 + 2 * ng),
    pnew(nx1 + 2 * ng, nx2 + 2 * ng) {

        m_ng = ng;
        m_nx1 = nx1;
        m_nx2 = nx2;
        m_wall_cond = wall_cond;

        m_x1max = x1max;
        m_x1min = x1min;
        m_x2max = x2max;
        m_x2min = x2min;

        m_dx = (x1max - x1min) / nx1;
        m_dy = (x2max - x2min) / nx2;

        gamma = i_gamma;
        R = i_R;
    }

    void saveOutput(int iter, double t, std::string fpref);
    // void SingleTimestep();
    void swapArrays();
    void fillHs(double nu);
    void pressureCorrect(int simiter, double rho, double tol, double t, double Uw, double omega,
        double nu);
    // void enforceBCs(double t, double Uw, double omega);

    // void new_SingleTimestep();
    void new_enforceBCs(double t, double Uw, double omega, double nu);
    // void fillDerivativeArrays();
    // void tempFillFluxesInefficient();
};

void
SimpleMesh::pressureCorrect(int simiter, double rho, double tol, double t, double Uw, double omega,
    double nu) {

    double diff = 1e3;
    double maxdiff, iterdiff, newp;
    int numiters = 0;
    int tot_iters = 40000;
    int i, j;

    if (simiter == 0) {
        tol = 1e-5;
    } else {
        tol = tol;
    }

    while (diff > tol) {
        // std::cout << "Pressure at top wall : " << p(50,10) << std::endl;

        maxdiff = 0.0;

        // For everywhere not on a boundary
        for (i = 2; i <= m_nx1 - 1; ++i) {
            for (j = 2; j <= m_nx2 - 1; ++j) {

                // if (i == 2 && j == 2) {
                //     std::cout << "newp = " << 0.25 * (p(i+2,j) + p(i-2,j) + p(i,j+2) + p(i,j-2)) <<
                //     std::endl;
                // }
                newp
                    = 0.25 *
                    (p(i + 2, j) + p(i - 2, j) + p(i, j + 2) + p(i, j - 2) -
                        4. * rho * m_dx * 0.5 *
                        ((Hx(i + 1, j) - Hx(i - 1, j)) + (Hy(i, j + 1) - Hy(i, j - 1))));

                iterdiff = std::abs((newp - p(i, j)));

                p(i, j) = newp;

                if (iterdiff > maxdiff) {
                    maxdiff = iterdiff;
                }
            }
        }

        // Left wall
        i = 1;
        for (j = 2; j <= m_nx2 - 1; ++j) {
            newp = 0.3333 *
                (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) + p(i, j + 2) + p(i, j - 2) -
                    4. * rho * m_dx * 0.5 *
                    (Hx(i + 1, j) + Hx(i, j) + Hy(i, j + 1) - Hy(i, j - 1)));

            iterdiff = std::abs((newp - p(i, j)));

            p(i, j) = newp;

            if (iterdiff > maxdiff) {
                maxdiff = iterdiff;
            }
        }

        // Right wall
        i = m_nx1;
        for (j = 2; j <= m_nx2 - 1; ++j) {
            newp = 0.3333 *
                (p(i - 1, j) - p(i + 1, j) + p(i - 2, j) + p(i, j + 2) + p(i, j - 2) -
                    4. * rho * m_dx * 0.5 *
                    (-(Hx(i, j) + Hx(i - 1, j)) + Hy(i, j + 1) - Hy(i, j - 1)));

            iterdiff = std::abs((newp - p(i, j)));

            p(i, j) = newp;

            if (iterdiff > maxdiff) {
                maxdiff = iterdiff;
            }
        }

        // Bottom wall
        j = 1;
        for (i = 2; i <= m_nx1 - 1; ++i) {
            newp = 0.3333 *
                (p(i + 2, j) + p(i - 2, j) + p(i, j + 2) + p(i, j + 1) - p(i, j - 1) -
                    4. * rho * m_dx * 0.5 *
                    (Hx(i + 1, j) - Hx(i - 1, j) + Hy(i, j + 1) + Hy(i, j)));

            iterdiff = std::abs((newp - p(i, j)));

            p(i, j) = newp;

            if (iterdiff > maxdiff) {
                maxdiff = iterdiff;
            }
        }

        // Top wall
        j = m_nx2;
        for (i = 2; i <= m_nx1 - 1; ++i) {
            newp = 0.3333 *
                (p(i + 2, j) + p(i - 2, j) - p(i, j + 1) + p(i, j - 1) + p(i, j - 2) -
                    4. * rho * m_dx * 0.5 *
                    (Hx(i + 1, j) - Hx(i - 1, j) - (Hy(i, j) + Hy(i, j - 1))));

            iterdiff = std::abs((newp - p(i, j)));

            p(i, j) = newp;

            if (iterdiff > maxdiff) {
                maxdiff = iterdiff;
            }
        }

        // Top left corner
        i = 1;
        j = m_nx2;

        newp
            = 0.5 *
            (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) - p(i, j + 1) + p(i, j - 1) +
                p(i, j - 2) -
                4. * rho * m_dx * 0.5 * (Hx(i + 1, j) + Hx(i, j) - (Hy(i, j) + Hy(i, j - 1))));
        iterdiff = std::abs((newp - p(i, j)));
        if (iterdiff > maxdiff) {
            maxdiff = iterdiff;
        }
        p(i, j) = newp;

        // Top right corner
        i = m_nx1;
        j = m_nx2;

        newp = 0.5 *
            (p(i - 1, j) - p(i + 1, j) + p(i - 2, j) - p(i, j + 1) + p(i, j - 1) +
                p(i, j - 2) +
                4. * rho * m_dx * 0.5 * (Hx(i, j) + Hx(i - 1, j) + Hy(i, j) + Hy(i, j - 1)));
        iterdiff = std::abs((newp - p(i, j)));
        if (iterdiff > maxdiff) {
            maxdiff = iterdiff;
        }
        p(i, j) = newp;

        // Bottom left corner
        i = 1;
        j = 1;

        newp = 0.5 *
            (p(i + 2, j) + p(i + 1, j) - p(i - 1, j) + p(i, j + 2) + p(i, j + 1) -
                p(i, j - 1) -
                4. * rho * m_dx * 0.5 * (Hx(i + 1, j) + Hx(i, j) + Hy(i, j + 1) + Hy(i, j)));
        iterdiff = std::abs((newp - p(i, j)));
        if (iterdiff > maxdiff) {
            maxdiff = iterdiff;
        }
        p(i, j) = newp;

        // Bottom right corner
        i = m_nx1;
        j = 1;

        newp
            = 0.5 *
            (p(i, j + 2) + p(i, j + 1) - p(i, j - 1) - p(i + 1, j) + p(i - 1, j) +
                p(i - 2, j) -
                4. * rho * m_dx * 0.5 * (Hy(i, j + 1) + Hy(i, j) - (Hx(i, j) + Hx(i - 1, j))));
        iterdiff = std::abs((newp - p(i, j)));
        if (iterdiff > maxdiff) {
            maxdiff = iterdiff;
        }
        p(i, j) = newp;

        diff = maxdiff;
        SimpleMesh::new_enforceBCs(t, Uw, omega, nu);

        if (numiters % 1000 == 0) {
            std::cout << "Pressure correction error : " << diff << std::endl;
        }
        numiters += 1;
    }

    // std::cout << "Finish pressure correction after " << numiters << " iterations, with error = " <<
    // diff << std::endl;
}

void
SimpleMesh::new_enforceBCs(double t, double Uw, double omega, double nu) {

    double T_inner, e_inner;
    int nx1t = m_nx1 + 2 * m_ng;
    int nx2t = m_nx2 + 2 * m_ng;

    // Left wall
    for (int j = 1; j <= m_nx2; ++j) {

        un(0, j) = -un(1, j);
        vn(0, j) = -vn(1, j);
        p(0, j) = p(1, j);
        pnew(0, j) = pnew(1, j);

        Hx(0, j) = -Hx(1, j);
        // Hy(0,j) = -Hy(1,j);
    }

    // Right wall
    for (int j = 1; j <= m_nx2; ++j) {

        un(nx1t - 1, j) = -un(nx1t - 2, j);
        vn(nx1t - 1, j) = -vn(nx1t - 2, j);
        p(nx1t - 1, j) = p(nx1t - 2, j);
        pnew(nx1t - 1, j) = pnew(nx1t - 2, j);

        Hx(nx1t - 1, j) = -Hx(nx1t - 2, j);
        // Hy(nx1t-1,j) = -Hy(nx1t-2,j);
    }

    // Bottom wall
    for (int i = 0; i <= m_nx1 + 1; ++i) {

        un(i, 0) = -un(i, 1);
        vn(i, 0) = -vn(i, 1);
        p(i, 0) = p(i, 1);
        pnew(i, 0) = pnew(i, 1);

        // Hx(i,0) = -Hx(i,1);
        Hy(i, 0) = -Hy(i, 1);
    }

    // Top wall
    for (int i = 0; i <= m_nx1 + 1; ++i) {

        un(i, nx2t - 1) = 2.0 * Uw * std::sin(omega * t) - un(i, nx2t - 2);
        // un(i,nx2t-1) = 2 * Uw - un(i,nx2t-2);
        vn(i, nx2t - 1) = -vn(i, nx2t - 2);
        p(i, nx2t - 1) = p(i, nx2t - 2);
        pnew(i, nx2t - 1) = pnew(i, nx2t - 2);

        // Hx(i,nx2t-1) = -Hx(i,nx2t-2);
        Hy(i, nx2t - 1) = -Hy(i, nx2t - 2);
    }
}

void
SimpleMesh::fillHs(double nu) {
    /*
    Fills arrays with their values to be used to calculate thermal
    diffusion and viscous forces
    */

    double uipoh, uimoh, ujpoh, ujmoh;
    double vipoh, vimoh, vjpoh, vjmoh;

    for (int i = 1; i <= m_nx1; ++i) {
        for (int j = 1; j <= m_nx2; ++j) {

            uipoh = (un(i + 1, j) + un(i, j)) / 2;
            uimoh = (un(i, j) + un(i - 1, j)) / 2;
            ujpoh = (un(i, j + 1) + un(i, j)) / 2;
            ujmoh = (un(i, j) + un(i, j - 1)) / 2;

            vipoh = (vn(i + 1, j) + vn(i, j)) / 2;
            vimoh = (vn(i, j) + vn(i - 1, j)) / 2;
            vjpoh = (vn(i, j + 1) + vn(i, j)) / 2;
            vjmoh = (vn(i, j) + vn(i, j - 1)) / 2;

            Hx(i, j) = -((uipoh * uipoh - uimoh * uimoh) / m_dx) -
                ((ujpoh * vjpoh - ujmoh * vjmoh) / m_dy) +
                nu *
                ((un(i + 1, j) - 2 * un(i, j) + un(i - 1, j)) / (m_dx * m_dx) +
                    (un(i, j + 1) - 2 * un(i, j) + un(i, j - 1)) / (m_dy * m_dy));

            Hy(i, j) = -((vjpoh * vjpoh - vjmoh * vjmoh) / m_dy) -
                ((uipoh * vipoh - uimoh * vimoh) / m_dx) +
                nu *
                ((vn(i + 1, j) - 2 * vn(i, j) + vn(i - 1, j)) / (m_dx * m_dx) +
                    (vn(i, j + 1) - 2 * vn(i, j) + vn(i, j - 1)) / (m_dy * m_dy));
        }
    }
}

void
SimpleMesh::saveOutput(int iter, double t, std::string fpref) {

    std::string fname = getFname(iter, fpref);
    double pipoh, pimoh, pjpoh, pjmoh;
    for (int i = 1; i <= m_nx1; ++i) {
        for (int j = 1; j <= m_nx2; ++j) {

            pipoh = (p(i + 1, j) + p(i, j)) / 2;
            pimoh = (p(i, j) + p(i - 1, j)) / 2;
            pjpoh = (p(i, j + 1) + p(i, j)) / 2;
            pjmoh = (p(i, j) + p(i, j - 1)) / 2;

            dpdx(i, j) = (pipoh - pimoh) / m_dx;
            dpdy(i, j) = (pjpoh - pjmoh) / m_dy;
        }
    }

    int nvar = 4;
    int nx1 = m_nx1 + 2 * m_ng;
    int nx2 = m_nx2 + 2 * m_ng;

    H5::H5File file(fname.c_str(), H5F_ACC_TRUNC);

    hsize_t dims[1];
    dims[0] = nx2 * nx1;

    H5::DataSpace dataspace(1, dims);

    H5::DataSet XVEL = file.createDataSet("VelX1", H5::PredType::NATIVE_DOUBLE, dataspace);
    XVEL.write(un.m_data.data(), H5::PredType::NATIVE_DOUBLE);

    H5::DataSet YVEL = file.createDataSet("VelX2", H5::PredType::NATIVE_DOUBLE, dataspace);
    YVEL.write(vn.m_data.data(), H5::PredType::NATIVE_DOUBLE);

    H5::DataSet PRES = file.createDataSet("Press", H5::PredType::NATIVE_DOUBLE, dataspace);
    PRES.write(p.m_data.data(), H5::PredType::NATIVE_DOUBLE);

    H5::DataSet DPDX = file.createDataSet("dpdx", H5::PredType::NATIVE_DOUBLE, dataspace);
    DPDX.write(dpdx.m_data.data(), H5::PredType::NATIVE_DOUBLE);

    H5::DataSet DPDY = file.createDataSet("dpdy", H5::PredType::NATIVE_DOUBLE, dataspace);
    DPDY.write(dpdy.m_data.data(), H5::PredType::NATIVE_DOUBLE);
}

void
SimpleMesh::swapArrays() {

    for (int i = 1; i <= m_nx1; ++i) {
        #pragma omp simd
        for (int j = 1; j <= m_nx2; ++j) {
            un(i, j) = unp1(i, j);
            vn(i, j) = vnp1(i, j);
        }
    }
}

// for (int i=1; i<=m_nx1; ++i) {
//     for (int j=1; j<=m_nx2; ++j) {

//         // Left Wall
//         if (i == 1 && j > 1 && j < m_nx2) {
//             newp = 1/4 * (p(i+2,j) + p(i-1,j) + p(i,j+2) + p(i,j-2) - 2 * rho * (m_dx) * (2 *
//             (Hx(i+1,j) - Hx(i,j)) + (Hy(i,j+1) - Hy(i,j-1))));
//         }

//         // Right wall
//         else if (i == m_nx1 && j > 1 && j < m_nx2) {
//             newp = 1/4 * (p(i+1,j) + p(i-2,j) + p(i,j+2) + p(i,j-2) - 2 * rho * (m_dx) * (2 *
//             (Hx(i,j) - Hx(i-1,j)) + (Hy(i,j+1) - Hy(i,j-1))));
//         }

//         // Bottom wall
//         else if (j == 1 && i > 1 && i < m_nx1) {
//             newp = 1/4 * (p(i+2,j) + p(i-2,j) + p(i,j+2) + p(i,j-1) - 2 * rho * (m_dx) *
//             ((Hx(i+1,j) - Hx(i-1,j)) + 2 * (Hy(i,j+1) - Hy(i,j))));
//         }

//         // Top wall
//         else if (j == m_nx2 && i > 1 && i < m_nx1) {
//             newp = 1/4 * (p(i+2,j) + p(i-2,j) + p(i,j+1) + p(i,j-2) - 2 * rho * (m_dx) *
//             ((Hx(i+1,j) - Hx(i-1,j)) + 2 * (Hy(i,j) - Hy(i,j-1))));
//         }

//         // Bottom left corner
//         else if (i == 1 && j == 1) {
//             newp = 1/4 * (p(i+2,j) + p(i-1,j) + p(i,j+2) + p(i,j-1) - 4 * rho * (m_dx) *
//             ((Hx(i+1,j) - Hx(i,j)) + (Hy(i,j+1) - Hy(i,j))));
//         }

//         // Bottom right corner
//         else if (i == 1 && j == m_nx2) {
//             newp = 1/4 * (p(i+2,j) + p(i-1,j) + p(i,j+1) + p(i,j-2) - 4 * rho * (m_dx) *
//             ((Hx(i+1,j) - Hx(i,j)) + (Hy(i,j) - Hy(i,j-1))));
//         }

//         // Top right corner
//         else if (i == m_nx1 && j == m_nx2) {
//             newp = 1/4 * (p(i+1,j) + p(i-2,j) + p(i,j+1) + p(i,j-2) - 4 * rho * (m_dx) *
//             ((Hx(i,j) - Hx(i-1,j)) + (Hy(i,j) - Hy(i,j-1))));
//         }

//         // Top left corner
//         else if (i == m_nx1 && j == 1) {
//             newp = 1/4 * (p(i+1,j) + p(i-2,j) + p(i,j+2) + p(i,j-1) - 4 * rho * (m_dx) *
//             ((Hx(i,j) - Hx(i-1,j)) + (Hy(i,j+1) - Hy(i,j))));
//         }

//         // Everything else
//         else {
//             newp = 1/4 * (p(i+2,j) + p(i-2,j) + p(i,j+2) + p(i,j-2) - 2 * rho * (m_dx) *
//             ((Hx(i+1,j) - Hx(i-1,j)) + (Hy(i,j+1) - Hy(i,j-1))));
//         }

//         iterdiff = std::abs((newp - p(i,j)));

//         p(i,j) = newp;

//         if (iterdiff > maxdiff) {
//             maxdiff = iterdiff;
//         }
//     }
// }

#endif