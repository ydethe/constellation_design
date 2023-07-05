#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const double t0_epoch=2319105600.0;

// Boltzmann constant (J/K)
const double kb = 1.380649e-23;

// Speed of light in vacuum (m/s)
const double clum = 299792458.0;

// Equatorial Earth radius in WGS84 (m)
const double Req = 6378137.0;

// Polar Earth radius in WGS84 (m)
const double Rpo = 6356752.3;

// Inverse of flattening in WGS84 (-), =Req/(Req-Rpo)
const double rf = 298.257223563;

// Earth sideral revolution pulsation (rad/s)
const double omega = 7.2921151467064e-5;

// Gravitationnal parameter (m^3/s^2)
const double mu = 3.986004418e14;

const double DAY_S = 86400.0;

const double T0 = 2451545.0;

const double z1_J1 = 1.616339347;

#endif
