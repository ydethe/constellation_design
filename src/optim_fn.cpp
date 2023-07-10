#include "constants.hpp"
#include "utils.hpp"
#include "optim_fn.h"

double
quadratic(double t, void *params)
{
    _find_event_data *d = (_find_event_data *)params;
    double t_epoch = t0_epoch + t;
    Matrix3d R1 = d->sat->getTEMEOrbitRotationMatrix(t, false);
    Matrix3d R2 = teme_transition_matrix(t_epoch, true, false);

    Vector3d M0 = *d->M0;
    Vector3d prx = *d->pos_rx;

    double v = prx.dot(R2 * R1 * M0);
    double J = d->s - v;

    return J;
}

double
quadratic_deriv(double t, void *params)
{
    _find_event_data *d = (_find_event_data *)params;
    double t_epoch = t0_epoch + t;
    Matrix3d R1 = d->sat->getTEMEOrbitRotationMatrix(t, false);
    Matrix3d R2 = teme_transition_matrix(t_epoch, true, false);
    Matrix3d dR1 = d->sat->getTEMEOrbitRotationMatrix(t, true);
    Matrix3d dR2 = teme_transition_matrix(t_epoch, true, true);

    Vector3d M0 = *d->M0;
    Vector3d prx = *d->pos_rx;

    double v = prx.dot(R2 * R1 * M0);
    double dJ = -prx.dot(dR2 * R1 * M0 + R2 * dR1 * M0);

    return dJ;
}

void quadratic_fdf(double t, void *params,
                   double *y, double *dy)
{
    *y = quadratic(t, params);
    *dy = quadratic_deriv(t, params);
}

double _culm_func(unsigned n, const double *x, double *grad, void *params)
{
    double J = quadratic(x[0], params);
    double dJ = quadratic_deriv(x[0], params);

    if (grad)
        grad[0]=dJ;

    return J;

}
