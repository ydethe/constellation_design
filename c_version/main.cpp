#include <iostream>

#include "constants.hpp"
#include "satellite.hpp"
#include "utils.hpp"
#include "simulate.hpp"

void test_theta_GMST1982();
void test_teme_orbit_mat();
void test_teme_mat();
void test_func();

int main(int argc, char **argv)
{
    const double lat = 80 * M_PI / 180;
    const double inc = 80 * M_PI / 180;
    const double nsat = 1;
    const double npla = 1;
    const double pha = 193;
    const double alt_km = 607.4;
    double t_blind;
    int nsat_max;
    Satellite sat;
    event_type events;
    VectorXd obs(6);
    obs = llavpa_to_itrf(0, lat, 0, 0, 0, 0);
    Vector3d M0(1,0,0);
    Vector3d M1(0.174783,0,0.984607);
    double J1,J2,s=0.998494;
    double dt = 0.1;
    double x[1]  ,grad[1] ;

    sat = Satellite(
                alt_km*1000+Req,
                inc,
                0,
                0,0);
    // sat.find_events(obs,0,1,&events);

    _find_event_data data = {s, &sat, &M0, &M1, true};

    x[0] =0;
    J1=_culm_func(1, x,grad, &data);

    x[0] =dt;
    J2=_culm_func(1, x,NULL, &data);

    // test_theta_GMST1982();
    // test_teme_orbit_mat();
    // test_teme_mat();
    // test_func();
    simulate(lat, inc, nsat, npla, pha, alt_km, &t_blind, &nsat_max);
    std::cout << t_blind << "," << nsat_max << std::endl;

}
