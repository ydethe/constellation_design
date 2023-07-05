#ifndef SATELLITE_H_
#define SATELLITE_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>

using Eigen::AngleAxisd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::seq;
using Eigen::Vector3d;
using Eigen::VectorXd;



typedef struct {
    double t_rise;
    double t_culmination;
    double t_set;
    double e_culmination;
} event_type;

class Satellite
{
private:
    double m_sma;
    double m_inc;
    double m_argp;
    double m_mano;
    double m_node;
    double m_sat_puls;
    VectorXd m_r1;
    VectorXd m_r2;
    VectorXd m_initial_itrf;

public:
    Satellite(
        double sma,
        double inc,
        double argp,
        double mano,
        double node);
    VectorXd getGeocentricITRFPositionAt(double td);
    AngleAxisd getTEMEOrbitRotationMatrix(double t);
    double getOrbitalPeriod();
    int _find_events(VectorXd obs, double t0, double elevation, event_type *events);

};

#endif
