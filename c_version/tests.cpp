#include <iostream>

#include "constants.hpp"
#include "utils.hpp"
#include "satellite.hpp"

void test_theta_GMST1982()
{
    double jd_ut1, fraction_ut1, theta, theta2, theta_dot;
    double dt = 0.1;

    jd_ut1 = 2460125.;
    fraction_ut1 = 0;
    theta_GMST1982(jd_ut1, fraction_ut1, &theta, &theta_dot);
    theta_GMST1982(jd_ut1, fraction_ut1 + dt / 86400., &theta2, &theta_dot);

    std::cout << theta_dot << std::endl;
    std::cout << (theta2 - theta) / dt << std::endl;
}

void test_teme_orbit_mat()
{
    const double lat = 80 * M_PI / 180;
    const double inc = 80 * M_PI / 180;
    const double alt_km = 600;
    const double sma = Req + alt_km * 1e3;
    const double raan = 0.;
    const double argp = 0.0;
    const double dt = 0.1;
    Matrix3d R1, R2, dR;
    VectorXd obs(6);
    obs = llavpa_to_itrf(0, lat, 0, 0, 0, 0);
    Satellite sat;

    sat = Satellite(
        sma,
        inc,
        0,
        0,
        raan * M_PI / 180);

    R1 = sat.getTEMEOrbitRotationMatrix(0, false).toRotationMatrix();
    R2 = sat.getTEMEOrbitRotationMatrix(dt, false).toRotationMatrix();
    dR = sat.getTEMEOrbitRotationMatrix(0, true).toRotationMatrix();

    std::cout << "dR" << std::endl;
    std::cout << dR << std::endl;
    std::cout << "(R2 - R1) / dt" << std::endl;
    std::cout << (R2 - R1) / dt << std::endl;
}

void test_teme_mat()
{
    const double dt = 0.1;
    Matrix3d R1, R2, dR;

    R1 = teme_transition_matrix(0, false, false).toRotationMatrix();
    R2 = teme_transition_matrix(dt, false, false).toRotationMatrix();
    dR = teme_transition_matrix(0, false, true).toRotationMatrix();

    std::cout << "dR" << std::endl;
    std::cout << dR << std::endl;
    std::cout << "(R2 - R1) / dt" << std::endl;
    std::cout << (R2 - R1) / dt << std::endl;
}
