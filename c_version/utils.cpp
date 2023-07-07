#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/Polynomials>
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
    *theta = fmod(fmod(jd_ut1, 1.0) + fraction_ut1 + fmod(g / DAY_S, 1.0), 1.) * M_PI * 2;
    *theta_dot = (1.0 + dg / (DAY_S * 36525.0)) * M_PI * 2 / 86400.;
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

AngleAxisd teme_transition_matrix(double t_epoch, bool reciprocal, bool derivative)
{
    double jd, fraction;
    double theta, theta_dot;
    AngleAxisd R;
    Matrix3d id3 = Matrix3d::Identity(3, 3);
    Matrix3d t_hat;
    Matrix3d t_t;

    time_to_jd_fraction(t_epoch, &jd, &fraction);

    theta_GMST1982(jd, fraction, &theta, &theta_dot);

    Vector3d uz(0.0, 0.0, 1.0);
    t_hat << 0, -uz(2), uz(1),
        uz(2), 0, -uz(0),
        -uz(1), uz(0), 0;
    t_t << uz(0) * uz(0), uz(1) * uz(0), uz(2) * uz(0),
        uz(0) * uz(1), uz(1) * uz(1), uz(2) * uz(1),
        uz(0) * uz(2), uz(1) * uz(2), uz(2) * uz(2);

    if (reciprocal)
        theta *= -1;

    if (derivative)
        R = theta_dot * (-sin(theta) * id3 + cos(theta) * t_hat + sin(theta) * t_t);
    else
        // R = AngleAxisd(theta, uz);
        R = cos(theta)*id3 + sin(theta)*t_hat + (1-cos(theta))*t_t;

    return R;
}

double compute_elevation_mask(double px_dbm, double nf_db, double eta, double fcarrier_mhz, double cn0_lim_dbhz, double alpha, double t0_k, double alt_km)
{
    // Parametres du probleme
    // ======================
    const double wl = clum / (fcarrier_mhz * 1e6);

    // Calculs des constantes K et Q
    // =============================
    const double sma = Req + alt_km * 1e3;
    const double K = (-36 + px_dbm - nf_db + 10 * log10((pow(eta * wl, 2) * pow(z1_J1 * sma, 4)) / (pow(4 * M_PI, 2) * kb * t0_k * pow(Req, 4))));
    const double Q = pow(10, (K - cn0_lim_dbhz) / 20);
    const double sma2 = pow(sma, 2);
    const double alpha2 = pow(alpha, 2);
    const double req2 = pow(Req, 2);
    const double q2 = pow(Q, 2);

    // Calcul de l'élévation
    // =====================
    Eigen::VectorXd coeff(7);
    coeff << alpha,
        2,
        -(
            (2 * alpha2 - 1) * sma2 + (1 - 2 * alpha2) * req2 + 2 * alpha * Q * Req + q2) /
            (alpha * sma2 - alpha * req2),
        -(4 * alpha * sma2 - 4 * alpha * req2 + 2 * Q * Req) / (alpha * sma2 - alpha * req2),
        ((alpha2 - 2) * sma2 + (2 - alpha2) * req2 + 2 * alpha * Q * Req) / (alpha * sma2 - alpha * req2),
        (2 * alpha * sma2 - 2 * alpha * req2 + 2 * Q * Req) / (alpha * sma2 - alpha * req2),
        1 / alpha;

    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    solver.compute(coeff);
    const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots();

    double re, im, re0;
    re0 = -1.;
    for (int i = 0; i < r.size(); i++)
    {
        re = std::real(r(i));
        im = std::imag(r(i));
        if (fabs(im) < 1e-6 && re > 0 && re < 1)
            re0 = re;
    }

    return asin(re0);
}

Vector3d geodetic_to_itrf(double lon, double lat, double h)
{
    /*Compute the Geocentric (Cartesian) Coordinates X, Y, Z
    given the Geodetic Coordinates lat, lon + Ellipsoid Height h

    Args:
        lon: Longitude (rad)
        lat: Latitude (rad)
        h: Altitude (m)

    Returns:
        A array of x, y, z coordinates (m)

    Examples:
        >>> x,y,z = geodetic_to_itrf(0,0,0)

    */
    Vector3d pos;
    double N = Req / sqrt(1 - (1 - pow(1 - 1 / rf, 2)) * pow(sin(lat), 2));
    double X = (N + h) * cos(lat) * cos(lon);
    double Y = (N + h) * cos(lat) * sin(lon);
    double Z = (pow(1 - 1 / rf, 2) * N + h) * sin(lat);

    pos << X, Y, Z;

    return pos;
}

Matrix3d build_local_matrix(Vector3d pos)
{
    /*Builds a ENV frame at a given position

    Args:
        pos: Position (m) of a point in ITRF

    Returns:
        Matrix whose columns are:

        * Local East vector
        * Local North vector
        * Local Vertical vector, orthogonal to xvec

    */
    // Local ENV for the observer
    Vector3d vert = pos;
    Vector3d uz(0, 0, 1);
    vert.normalize();

    Vector3d east = uz.cross(pos);
    east.normalize();

    Vector3d north = vert.cross(east);

    Matrix3d env;

    env(seq(0, 2), 0) = east;
    env(seq(0, 2), 1) = north;
    env(seq(0, 2), 2) = vert;

    return env;
}

VectorXd llavpa_to_itrf(double lon, double lat, double alt, double vel, double vs, double va)
{
    /*Converts a LLAVPA into ITRF
    A LLAVPA is an array with
    longitude, latitude, altitude (WGS84) and velocity, slope of velocity, azimut of velocity
    Velocity is the speed of sat in ITRF frame

    Args:
        llavpa: lon,lat,alt,vel,vs,va

    Returns:
        An ITRF position and velocity

    */
    VectorXd pv(6);
    pv(seq(0, 2)) = geodetic_to_itrf(lon, lat, alt);

    // Local ENV for the satellite
    Matrix3d M = build_local_matrix(pv(seq(0, 2)));

    Vector3d local_v;
    double ve = cos(vs) * sin(va) * vel;
    double vn = cos(vs) * cos(va) * vel;
    double vv = sin(vs) * vel;
    local_v << ve, vn, vv;

    pv(seq(3, 5)) = M * local_v;

    return pv;
}
