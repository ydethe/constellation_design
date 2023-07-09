#include <math.h>
#include <iostream>

#include "constants.hpp"
#include "satellite.hpp"
#include "utils.hpp"

#define RISE_EVENT 0
#define SET_EVENT 1

void analyse_timeline(const int nb_sat_ini, const std::vector<int> event_types, const std::vector<double> event_dates, const double total_sim_time, double *t_blind, int *nsat_max)
{
    int nsat, sorted_index;
    int ns;
    double last_blind_date;
    std::vector<size_t> indices;
    indices = argsort(event_dates);

    ns = event_dates.size();
    nsat = nb_sat_ini;
    *t_blind = 0;
    *nsat_max = 0;
    last_blind_date = 0;
    for (int i = 0; i < ns; i++)
    {
        sorted_index = indices[i];

        if (event_types[sorted_index] == RISE_EVENT)
        {
            if (nsat == 0)
                *t_blind += event_dates[sorted_index] - last_blind_date;
            nsat += 1;
            std::cout << event_dates[sorted_index] << " rise, nsat=" << nsat << ", t_blind=" << *t_blind << std::endl;
        }
        else if (event_types[sorted_index] == SET_EVENT)
        {
            nsat -= 1;
            if (nsat == 0)
                last_blind_date = event_dates[sorted_index];
            else if (nsat < 0) {
                std::cerr << "nsat:" << nsat << std::endl;
                throw;
            }
            std::cout << event_dates[sorted_index] << " set, nsat=" << nsat << ", t_blind=" << *t_blind << std::endl;
        }
        else {
            std::cerr << "unknown event:" << event_types[sorted_index]  << std::endl;
            throw;
        }

        if (*nsat_max < nsat)
            *nsat_max = nsat;
    }

    if (*nsat_max == 0)
        *t_blind = total_sim_time;
}

void simulate(double lat, double inc, int nsat, int npla, int pha, double alt_km, double *t_blind, int *nsat_max)
{
    std::vector<int> event_types;
    std::vector<double> event_dates;
    int s = nsat / npla;
    int status;
    int satellite_id, nb_sat;
    double elev_mask;
    Satellite sat;
    event_type events;
    VectorXd obs(6), pv_sat(6);
    obs = llavpa_to_itrf(0, lat, 0, 0, 0, 0);

    const double sma = Req + alt_km * 1e3;
    const double px_dbm = 36;
    const double nf_db = 3;
    const double eta = 0.7;
    const double fcarrier_mhz = 2000;
    const double cn0_lim_dbhz = 46.46595211;
    const double alpha = 0.0305;
    const double t0_k = 290;
    // const double t_sim = 5 * 86400;
    const double t_sim = 20000;
    const double firstraan = 0.;
    const double argp = 0.0;
    double raan, meanAnomaly, t_start;

    nb_sat = 0;
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

            t_start = 0;
            while (true)
            {
                sat.find_events(obs, t_start, elev_mask, &events);
                t_start = events.t_culmination + events.tup_max;
                std::cout << events.t_culmination << ", " << events.t_rise << ", " << events.is_initially_visible << std::endl;
                if (events.is_initially_visible)
                    nb_sat++;
                if (events.t_culmination > t_sim)
                    break;
                if (events.t_rise > 0)
                {
                    event_types.push_back(RISE_EVENT);
                    event_dates.push_back(events.t_rise);
                    event_types.push_back(SET_EVENT);
                    event_dates.push_back(events.t_set);
                }
            }
        }
    }

    std::cout << event_dates.size() << std::endl;
    analyse_timeline(nb_sat, event_types, event_dates, t_sim, t_blind, nsat_max);
}
