#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>

using Eigen::MatrixXd;
using Eigen::VectorXd;

double anomaly_mean_to_ecc(double ecc, double M);

double anomaly_ecc_to_true(double ecc, double E);

double anomaly_mean_to_true(double ecc, double M);

double anomaly_true_to_ecc(double ecc, double v);

VectorXd orbital_to_teme(double sma, double ecc, double argp, double inc, double mano, double node);

VectorXd teme_to_itrf(double t_epoch, VectorXd pv_teme);

VectorXd itrf_to_teme(double t_epoch, VectorXd pv_itrf);

AngleAxisd teme_transition_matrix(double t_epoch, bool reciprocal);

#endif
