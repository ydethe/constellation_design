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

    R1 = sat.getTEMEOrbitRotationMatrix(0, false);
    R2 = sat.getTEMEOrbitRotationMatrix(dt, false);
    dR = sat.getTEMEOrbitRotationMatrix(0, true);

    std::cout << "dR" << std::endl;
    std::cout << dR << std::endl;
    std::cout << "(R2 - R1) / dt" << std::endl;
    std::cout << (R2 - R1) / dt << std::endl;
}

void test_teme_mat()
{
    const double dt = 0.1;
    Matrix3d R1, R2, dR;

    R1 = teme_transition_matrix(0, false, false);
    R2 = teme_transition_matrix(dt, false, false);
    dR = teme_transition_matrix(0, false, true);

    std::cout << "\ndR" << std::endl;
    std::cout << dR << std::endl;
    std::cout << "\n(R2 - R1) / dt" << std::endl;
    std::cout << (R2 - R1) * (1 / dt) << std::endl;
}

double _culm_func(unsigned n, const double *x, double *grad, void *my_func_data);

typedef struct
{
    double s;
    Satellite *sat;
    Vector3d *M0;
    Vector3d *pos_rx;
    bool minimize;
} _find_event_data;

void test_func()
{
    const double s = 0.998429;
    const double dt = 0.1;
    double x[1], grad[1];
    double Jm, J0, Jp;
    Vector3d M0(-0.375765, 0.0975621, 0.921565);
    Vector3d M1(0.174783, 0, 0.984607);
    Satellite sat(7e6, 3 * M_PI / 8, 0.5, 1.0, 1.5);
    double max_err, err;

    _find_event_data data = {s, &sat, &M0, &M1, true};

    max_err = 0;
    for (double t = 0; t < 1000; t += 100)
    {
        x[0] = 1000 - dt;
        Jm = _culm_func(1, x, NULL, &data);

        x[0] = 1000;
        J0 = _culm_func(1, x, grad, &data);

        x[0] = 1000 + dt;
        Jp = _culm_func(1, x, NULL, &data);

        err = (Jp - Jm) / (2 * dt) - grad[0];

        if (err > max_err)
            max_err = err;
    }

    std::cout << max_err << std::endl;
}
