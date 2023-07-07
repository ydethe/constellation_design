#include <iostream>

#include "satellite.hpp"
#include "utils.hpp"
#include "simulate.hpp"

void test_theta_GMST1982();
void test_teme_orbit_mat();
void test_teme_mat();

int main(int argc, char **argv)
{
    const double lat = 80 * M_PI / 180;
    const double inc = 80 * M_PI / 180;
    const double nsat = 2;
    const double npla = 2;
    const double pha = 193;
    const double alt_km = 607.4;
    double t_blind;
    int nsat_max;

    test_theta_GMST1982();
    test_teme_mat();
    // simulate(lat, inc, nsat, npla, pha, alt_km, &t_blind, &nsat_max);

    // std::cout << t_blind << "," << nsat_max << std::endl;
}
