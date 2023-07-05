#include <cmath>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <nlopt.h>

#include "constants.hpp"

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::seq;
using Eigen::Vector3d;
using Eigen::VectorXd;

typedef struct
{
    double ecc, M;
} _anomaly_mean_to_ecc_data;

double _anomaly_mean_to_ecc_fun(unsigned n, const double *x, double *grad, void *data)
{
    _anomaly_mean_to_ecc_data *d = (_anomaly_mean_to_ecc_data *)data;
    double E = x[0];
    double J = E - d->ecc * sin(E) - d->M;
    if (grad)
    {
        grad[0] = J * (1 - d->ecc * cos(E));
    }
    return pow(J, 2) / 2;
}

double anomaly_mean_to_ecc(double ecc, double M)
{
    if (ecc == 0)
        return M;

    double lb[2] = {-M_PI, M_PI};
    _anomaly_mean_to_ecc_data data = {ecc, M};
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LD_MMA, 1); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_min_objective(opt, _anomaly_mean_to_ecc_fun, &data);

    nlopt_set_xtol_rel(opt, 1e-4);

    double x[1] = {0}; /* `*`some` `initial` `guess`*` */
    double minf;       /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */
    if (nlopt_optimize(opt, x, &minf) < 0)
    {
        printf("nlopt failed!\n");
    }

    nlopt_destroy(opt);

    return x[0];
}

double anomaly_ecc_to_true(double ecc, double E)
{
    double tv2 = sqrt((1 + ecc) / (1 - ecc)) * tan(E / 2);
    return 2 * atan(tv2);
}

double anomaly_mean_to_true(double ecc, double M)
{
    double E = anomaly_mean_to_ecc(ecc, M);
    double v = anomaly_ecc_to_true(ecc, E);
    return v;
}

double anomaly_true_to_ecc(double ecc, double v)
{
    double tE2 = sqrt((1 - ecc) / (1 + ecc)) * tan(v / 2);
    return 2 * atan(tE2);
}

void time_to_jd_fraction(double t_epoch, double *jd, double *fraction)
{
    /*

    Args:
        t_epoch: Time since 31/12/1949 00:00 UT (s)

    */
    double epoch = t_epoch / 86400;

    *fraction = fmod(epoch, 1.0);
    double whole = epoch - *fraction;
    *jd = whole + 2433281.5;
}

void theta_GMST1982(double jd_ut1, double fraction_ut1, double *theta, double *theta_dot)
{
    double t = (jd_ut1 - T0 + fraction_ut1) / 36525.0;
    double g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t;
    double dg = 8640184.812866 + (0.093104 * 2.0 + (-6.2e-6 * 3.0) * t) * t;
    *theta = fmod(fmod(jd_ut1, 1.0) + fraction_ut1 + fmod(g / DAY_S, 1.0), 1.) * M_PI*2;
    *theta_dot = (1.0 + dg / (DAY_S * 36525.0)) * M_PI*2;
}

VectorXd orbital_to_teme(double sma, double ecc, double argp, double inc, double mano, double node)
{
    // https://en.wikipedia.org/wiki/True_anomaly#From_the_mean_anomaly
    double p = sma * (1 - pow(ecc, 2));
    double tano = anomaly_mean_to_true(ecc, mano);
    double r = p / (1 + ecc * cos(tano));
    double n = sqrt(mu / pow(sma, 3));

    double x = r * (cos(node) * cos(argp + tano) - sin(node) * cos(inc) * sin(argp + tano));
    double y = r * (sin(node) * cos(argp + tano) + cos(node) * cos(inc) * sin(argp + tano));
    double z = r * sin(inc) * sin(argp + tano);

    double E = anomaly_true_to_ecc(ecc, tano);

    double dE = n / (1 - cos(E) * ecc);
    double rr = sqrt((1 + ecc) / (1 - ecc));
    double dtano = dE * rr * pow(cos(tano / 2) / cos(E / 2), 2);
    double vr = ecc * p * sin(tano) * dtano / pow(ecc * cos(tano) + 1, 2);

    double vx = r * (-cos(node) * dtano * sin(tano + argp) - cos(inc) * sin(node) * dtano * cos(tano + argp)) + vr * (cos(node) * cos(tano + argp) - cos(inc) * sin(node) * sin(tano + argp));
    double vy = r * (cos(inc) * cos(node) * dtano * cos(tano + argp) - sin(node) * dtano * sin(tano + argp)) + vr * (cos(inc) * cos(node) * sin(tano + argp) + sin(node) * cos(tano + argp));
    double vz = sin(inc) * vr * sin(tano + argp) + sin(inc) * r * dtano * cos(tano + argp);

    VectorXd pv(6);

    pv << x, y, z, vx, vy, vz;

    return pv;
}

VectorXd teme_to_itrf(double t_epoch, VectorXd pv_teme)
{
    double jd, fraction;
    double theta, theta_dot;

    time_to_jd_fraction(t_epoch, &jd, &fraction);

    theta_GMST1982(jd, fraction, &theta, &theta_dot);

    Vector3d uz(0.0, 0.0, 1.0);
    Vector3d angular_velocity = uz * (-theta_dot / 86400.0);

    Matrix3d R;
    R = AngleAxisd(-theta, uz);

    Vector3d rTEME = pv_teme(seq(0, 2));
    Vector3d vTEME = pv_teme(seq(3, 5));

    Vector3d rITRF = R * rTEME;
    Vector3d vITRF = R * vTEME + angular_velocity.cross(rITRF);

    VectorXd pv(6);
    pv(seq(0, 2)) = rITRF;
    pv(seq(3, 5)) = vITRF;

    return pv;
}

VectorXd itrf_to_teme(double t_epoch, VectorXd pv_itrf)
{
    double jd, fraction;
    double theta, theta_dot;

    time_to_jd_fraction(t_epoch, &jd, &fraction);

    theta_GMST1982(jd, fraction, &theta, &theta_dot);

    Vector3d uz(0.0, 0.0, 1.0);
    Vector3d angular_velocity = uz * (-theta_dot / 86400.0);

    Matrix3d R;
    R = AngleAxisd(theta, uz);

    Vector3d rITRF = pv_itrf(seq(0, 2));
    Vector3d vITRF = pv_itrf(seq(3, 5));

    Vector3d rTEME = R * rITRF;
    Vector3d vTEME = R * vITRF - R * angular_velocity.cross(rITRF);

    VectorXd pv(6);
    pv(seq(0, 2)) = rTEME;
    pv(seq(3, 5)) = vTEME;

    return pv;
}

AngleAxisd teme_transition_matrix(double t_epoch, bool reciprocal)
{
    double jd, fraction;
    double theta, theta_dot;
    AngleAxisd R;

    time_to_jd_fraction(t_epoch, &jd, &fraction);

    theta_GMST1982(jd, fraction, &theta, &theta_dot);

    Vector3d uz(0.0, 0.0, 1.0);

    if (reciprocal)
        R = AngleAxisd(-theta, uz);
    else
        R = AngleAxisd(theta, uz);

    return R;
}
