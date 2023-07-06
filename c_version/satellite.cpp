#include <nlopt.h>
#include <iostream>
#include <iomanip>

#include "satellite.hpp"
#include "constants.hpp"
#include "utils.hpp"

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::seq;
using Eigen::Vector3d;
using Eigen::VectorXd;

Satellite::Satellite()
{
    m_sma = 0.;
    m_inc = 0.;
    m_argp = 0.;
    m_mano = 0.;
    m_node = 0.;
}

Satellite::Satellite(
    double sma,
    double inc,
    double argp,
    double mano,
    double node)
{
    m_sma = sma;
    m_inc = inc;
    m_argp = argp;
    m_mano = mano;
    m_node = node;

    VectorXd pv_teme = orbital_to_teme(sma, 0, argp, inc, mano, node);

    Vector3d pos = pv_teme(seq(0, 2));
    Vector3d vel = pv_teme(seq(3, 5));

    double a = pos.norm();
    m_sat_puls = sqrt(mu / pow(a, 3));

    Vector3d n = pos.cross(vel);
    n.normalize();

    vel = n.cross(pos) * m_sat_puls;

    m_r1.resize(6);
    m_r2.resize(6);
    m_initial_itrf.resize(6);

    m_r1(seq(0, 2)) = pos;
    m_r1(seq(3, 5)) = vel;

    m_r2(seq(0, 2)) = n.cross(pos);
    m_r2(seq(3, 5)) = n.cross(vel);

    m_initial_itrf = teme_to_itrf(t0_epoch, m_r1);
}

VectorXd Satellite::getGeocentricITRFPositionAt(double td)
{
    double t_epoch = t0_epoch + td;
    double th = m_sat_puls * td;

    VectorXd newpv_teme = m_r1 * cos(th) + m_r2 * sin(th);
    VectorXd newpv = teme_to_itrf(t_epoch, newpv_teme);

    return newpv;
}

AngleAxisd Satellite::getTEMEOrbitRotationMatrix(double t)
{
    /*Compute the transition matrix R from ITRF to TEME, i.e.
    for X in ITRF, R @ X is the vector in TEME

    Args:
        t: Time from Satellite's tsync (s)

    Returns:
        See `blocksim.utils.rotation_matrix`
    */
    VectorXd pv0 = this->getGeocentricITRFPositionAt(0);
    VectorXd pv_teme = itrf_to_teme(t0_epoch, pv0);
    Vector3d pos = pv_teme(seq(0, 2));
    Vector3d vel = pv_teme(seq(3, 5));
    double angle = m_sat_puls * t;
    Vector3d axe = pos.cross(vel);
    axe.normalize();

    AngleAxisd R(angle, axe);

    return R;
}

double Satellite::getOrbitalPeriod()
{
    return 2 * M_PI / m_sat_puls;
}

typedef struct
{
    double s;
    Satellite *sat;
    Vector3d *M0;
    Vector3d *pos_rx;
    bool minimize;
} _find_event_data;

double _culm_func(unsigned n, const double *x, double *grad, void *my_func_data)
{
    double t = x[0];
    _find_event_data *d = (_find_event_data *)my_func_data;
    double t_epoch = t0_epoch + t;
    AngleAxisd R1 = d->sat->getTEMEOrbitRotationMatrix(t);
    AngleAxisd R2 = teme_transition_matrix(t_epoch, true);

    Vector3d M0 = *d->M0;
    Vector3d prx = *d->pos_rx;

    double v = prx.dot(R2 * R1 * M0);

    double J = d->s - v;

    if (d->minimize)
        return J;
    else
        return J * J;
}

void Satellite::find_events(VectorXd obs, double t0, double elevation, event_type *events)
{
    // ===========================================================
    // Computing the problem's parameters
    // ===========================================================
    const double Torb = this->getOrbitalPeriod();
    VectorXd pv0 = this->getGeocentricITRFPositionAt(0);
    const double r = pv0(seq(0, 2)).norm();
    int opt_status;

    VectorXd pv_teme = itrf_to_teme(t0_epoch, pv0);
    Vector3d M0 = pv_teme(seq(0, 2));
    M0.normalize();
    Vector3d M1 = obs(seq(0, 2));
    M1.normalize();

    const double d = -Req * sin(elevation) + sqrt(r * r - Req * Req * pow(cos(elevation), 2));
    const double beta = acos((d * d + r * r - Req * Req) / (2 * r * d));
    const double alpha = M_PI / 2 - (elevation + beta);
    const double Tup_max = Torb * alpha / M_PI;
    const double s = cos(alpha);
    _find_event_data data = {s, this, &M0, &M1, true};

    // ===========================================================
    // Setting up te optimizer
    // ===========================================================
    double lb[2];
    double x[1];
    double minf;
    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_COBYLA, 1); /* algorithm and dimensionality */
    nlopt_set_min_objective(opt, _culm_func, &data);
    nlopt_set_xtol_rel(opt, 1e-4);

    // ===========================================================
    // Calling the optimizer to find the culmination
    // ===========================================================
    lb[0] = t0;
    lb[1] = t0 + 1.2 * Torb;
    x[0] = t0 + Torb / 2;
    data.minimize = true;
    nlopt_set_lower_bounds(opt, lb);
    opt_status = nlopt_optimize(opt, x, &minf);
    if (opt_status < 0)
    {
        std::cerr << nlopt_get_errmsg(opt) << std::endl;
        return;
    }

    const double alpha_max = acos(s - minf);
    const double d_max = sqrt(Req * Req + r * r - 2 * Req * r * cos(alpha_max));
    const double elev_max = -asin((d_max * d_max + Req * Req - r * r) / (2 * Req * d_max));

    events->t_culmination = x[0];
    events->e_culmination = elev_max;
    events->tup_max = Tup_max;
    events->t_rise = -1;
    events->t_set = -1;

    if (elev_max < elevation)
        return;

    // ===========================================================
    // Searching rise time
    // ===========================================================
    data.minimize = false;
    lb[0] = events->t_culmination - Tup_max;
    lb[1] = events->t_culmination;
    x[0] = events->t_culmination - Tup_max / 2;
    nlopt_set_lower_bounds(opt, lb);
    opt_status = nlopt_optimize(opt, x, &minf);
    if (opt_status < 0)
    {
        std::cerr << nlopt_get_errmsg(opt) << std::endl;
        return;
    }

    events->t_rise = x[0];

    // ===========================================================
    // Searching set time
    // ===========================================================
    data.minimize = false;
    lb[0] = events->t_culmination;
    lb[1] = events->t_culmination + Tup_max;
    x[0] = events->t_culmination + Tup_max / 2;
    nlopt_set_lower_bounds(opt, lb);
    opt_status = nlopt_optimize(opt, x, &minf);
    if (opt_status < 0)
    {
        std::cerr << nlopt_get_errmsg(opt) << std::endl;
        return;
    }

    events->t_set = x[0];

    nlopt_destroy(opt);

    return;
}

// int main(int argc, char **argv)
// {
//     int status;
//     event_type events;
//     VectorXd obs(6);
//     obs = llavpa_to_itrf(0, M_PI/4, 0, 0, 0, 0);
//     Satellite sat(7e6, M_PI / 3, 0.5, 1.0, 1.5);
//     status = sat.find_events(obs, 0, 0.15, &events);
//     std::cout << status << "," << events.t_rise << "," << events.t_culmination << "," << events.t_set << std::endl;

//     double elev_mask = compute_elevation_mask(36, 3, 0.7, 2000., 46.46595211, 0.0305, 290, 1000);
//     std::cout << elev_mask << std::endl;
// }
