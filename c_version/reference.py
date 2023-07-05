from datetime import datetime, timezone

from numpy import pi
import numpy as np
from blocksim.satellite.Satellite import CircleSatellite


obs = np.array([4517590.87884893, 0.0, 4487348.40886592, 0.0, 0.0, 0.0])
t0 = datetime(2023, 6, 27, 12, 0, 0, tzinfo=timezone.utc)
sat = CircleSatellite.fromOrbitalElements("sat", t0, 7e6, pi / 3, 0.5, 1.0, 1.5)
Torb = sat.orbit_period.total_seconds()
print(Torb)
events = sat.find_events(obs, 0, 1.2 * Torb, 0.15)
print(events)
