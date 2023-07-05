#include <math.h>

#include "constants.hpp"
#include "satellite.hpp"
#include "utils.hpp"

void simulate(double lat, double inc, int nsat, int npla, int pha, double alt_km)
{
    int s = nsat / npla;
    int status;
    double elev_mask;
    Satellite sat;
    event_type events;
    VectorXd obs(6), pv_sat(6);
    obs << 4517590.87884893, 0., 4487348.40886592, 0., 0., 0.;
    const double sma = Req + alt_km * 1e3;
    const double px_dbm = 36;
    const double nf_db = 3;
    const double eta = 0.7;
    const double fcarrier_mhz = 2000;
    const double cn0_lim_dbhz = 46.46595211;
    const double alpha = 0.0305;
    const double t0_k = 290;
    const double t_sim = 5 * 86400;
    const double firstraan = 0.;
    const double argp = 0.0;
    double raan, meanAnomaly;

    for (int idxP = 0; idxP < npla; idxP++)
    {
        raan = firstraan * 180 / M_PI + idxP * 360.0 / npla;
        for (int idxS = 0; idxS < s; idxS++)
        {
            meanAnomaly = idxP * pha * 360.0 / nsat + idxS * 360.0 / s;
            elev_mask = compute_elevation_mask(px_dbm, nf_db, eta, fcarrier_mhz, cn0_lim_dbhz, alpha, t0_k, alt_km);

            sat = Satellite(
                sma,
                inc,
                argp,
                meanAnomaly * M_PI / 180,
                raan * M_PI / 180);

            status = sat._find_events(obs, 0, elev_mask, &events);
        }
    }
}
