#ifndef UTILS_H_
#define UTILS_H_

#include <vector>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using Eigen::MatrixXd;
using Eigen::VectorXd;


/**
 * Argsort(currently support ascending sort)
 * @tparam T array element type
 * @param array input array
 * @return indices w.r.t sorted array
 */
template<typename T>
std::vector<size_t> argsort(const std::vector<T> &array) {
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&array](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return array[left] < array[right];
              });

    return indices;
}

double anomaly_mean_to_ecc(double ecc, double M);

double anomaly_ecc_to_true(double ecc, double E);

double anomaly_mean_to_true(double ecc, double M);

double anomaly_true_to_ecc(double ecc, double v);

VectorXd orbital_to_teme(double sma, double ecc, double argp, double inc, double mano, double node);

VectorXd teme_to_itrf(double t_epoch, VectorXd pv_teme);

VectorXd itrf_to_teme(double t_epoch, VectorXd pv_itrf);

AngleAxisd teme_transition_matrix(double t_epoch, bool reciprocal, bool derivative);

double compute_elevation_mask(double px_dbm, double nf_db, double eta, double fcarrier_mhz, double cn0_lin_dbhz, double alpha, double t0_k, double alt_km);

Matrix3d build_local_matrix(Vector3d pos);

Vector3d geodetic_to_itrf(double lon,double lat,double h);

VectorXd llavpa_to_itrf(double lon, double lat, double alt, double vel, double vs, double va);

#endif
