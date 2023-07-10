#ifndef OPTIM_FN_H_
#define OPTIM_FN_H_

#include "satellite.hpp"

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::seq;
using Eigen::Vector3d;
using Eigen::VectorXd;

typedef struct
{
    double s;
    Satellite *sat;
    Vector3d *M0;
    Vector3d *pos_rx;
} _find_event_data;

double quadratic (double x, void *params);
double quadratic_deriv (double x, void *params);
void quadratic_fdf (double x, void *params, double *y, double *dy);
double _culm_func(unsigned n, const double *x, double *grad, void *params);

#endif
