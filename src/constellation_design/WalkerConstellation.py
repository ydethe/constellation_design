from datetime import datetime, timezone
from typing import List

from numpy import pi
from blocksim.satellite import CircleSatellite


# Adapted from https://gitlab.isae-supaero.fr/jsatorb-dev/jsatorb-common/-/blob/master/src/constellation/WalkerConstellation.py


class WalkerConstellation:
    """Class that generates a constellation according to the Walker Delta Pattern t:p:f

    Args:
        name: Name of the constellation
        sma : Semi-major axis (m)
        inc: Inclination of orbital planes (deg)
        firstraan: RAAN of the first orbital plane (deg)
        t: Number of satellites
        p: Number of equally spaced planes
        f: Relative spacing between satellites in adjacent planes

    """

    def __init__(self, name: str, sma: float, inc: float, firstraan: float, t: int, p: int, f: int):
        # Name
        self.name = str(name)

        # SMA or radius
        self.sma = float(sma)

        # Inclination
        self.inc = float(inc)

        # RAAN of first orbit
        self.firstraan = float(firstraan)

        # Number of satellites
        self.t = int(t)

        # Number of equally spaced planes
        self.p = int(p)

        # Relative spacing between satellites in adjacent planes
        self.f = int(f)

    def generate(self, tsync: datetime) -> List[CircleSatellite]:
        satellites = []

        # Number of satellites per plane
        s = self.t / self.p
        if s.is_integer():
            s = int(s)
        else:
            raise (ValueError("Number of satellites per plane (t/p) should be integer"))

        for idxP in range(self.p):
            raan = self.firstraan + idxP * 360.0 / self.p
            for idxS in range(s):
                meanAnomaly = idxP * self.f * 360.0 / self.t + idxS * 360.0 / s

                nameSat = self.name + "_" + str(s * idxP + idxS)

                satCur = CircleSatellite.fromOrbitalElements(
                    name=nameSat,
                    tsync=tsync,
                    a=self.sma,
                    inc=self.inc * pi / 180,
                    argp=0.0,
                    mano=meanAnomaly * pi / 180,
                    node=raan * pi / 180,
                )

                satellites.append(satCur)

        return satellites


if __name__ == "__main__":
    argConst = {
        "name": "myConst",
        "sma": 7000000.0,
        "inc": 15.0,
        "firstraan": 10.0,
        "t": 9,
        "p": 3,
        "f": 1,
    }
    walkerConst = WalkerConstellation(**argConst)
    satellites = walkerConst.generate(tsync=datetime.now().replace(tzinfo=timezone.utc))

    sat = satellites[0]
    pv = sat.getGeocentricITRFPositionAt(0)
    print(pv)
