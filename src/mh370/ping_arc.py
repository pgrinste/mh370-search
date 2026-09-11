"""Geometry of the 'seventh arc' from the final ATOC handshake (2014-03-09 00:19:54 UTC).

The BTO measurement fixes the *slant range* between the aircraft and the
Inmarsat-3 F1 geostationary satellite (sub-satellite point at 0 deg, 64.5 deg E).
We model the constant-range locus exactly as the intersection of a sphere
centred on the satellite with radius R (slant range) and the WGS84 ellipsoid,
rather than as a ground circle about the sub-satellite point.

R is fitted so that the locus passes through the documented anchor point
(34.13 S, 93.95 E) - the seventh-arc longitude at 34.13 S for an aircraft
altitude of 20,000 ft (public analysis; see data/curated/data_dictionary.md).

Note on the old ground-circle model: with both models calibrated on the same
anchor they agree to within ~10 km over the southern band (the difference is
purely ellipsoidal), so this change removes a systematic bias by construction
rather than shifting the result. See tests/test_ping_arc.py.

Validation: JACC stated the arc's total extent reaches from latitude 20 S to
39 S (MarineLink, 20 Jun 2014). See tests/test_ping_arc.py.
"""

import math
from dataclasses import dataclass

from .geodesy import ecef_to_geodetic, geodetic_to_ecef, haversine_nm

SATELLITE_SUBPOINT = (0.0, 64.5)   # Inmarsat-3 F1 geostationary slot
GEO_RADIUS_KM = 42164.0            # geostationary orbit radius from Earth centre
ARC_ANCHOR = (-34.13, 93.95)       # documented point on the seventh arc


def satellite_ecef():
    """Inmarsat-3 F1 position in ECEF km (geostationary slot, equatorial)."""
    l = math.radians(SATELLITE_SUBPOINT[1])
    return (
        GEO_RADIUS_KM * math.cos(l),
        GEO_RADIUS_KM * math.sin(l),
        0.0,
    )


@dataclass(frozen=True)
class Arc:
    sat_ecef: tuple
    slant_range_km: float

    def lat_at_lon(self, lon_deg):
        """Latitude (deg) of the southern-branch locus at a given longitude.

        Solves |P(lat, lon) - S| = R for lat by bisection on [-80, -10].
        Returns None where the locus does not cross that longitude in the south.
        """
        s = self.sat_ecef
        r = self.slant_range_km

        def resid(lat):
            p = geodetic_to_ecef(lat, lon_deg, 0.0)
            return math.dist(p, s) - r

        lo, hi = -80.0, -10.0
        if resid(lo) * resid(hi) > 0:
            return None
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if resid(lo) * resid(mid) <= 0:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    def surface_points(self, lon_min=60.0, lon_max=112.0, step_deg=0.1):
        """Sample the southern branch as a polyline sorted by longitude."""
        pts = []
        n = int(round((lon_max - lon_min) / step_deg)) + 1
        for i in range(n):
            lon = lon_min + i * step_deg
            lat = self.lat_at_lon(lon)
            if lat is not None:
                pts.append((lat, lon))
        return pts

    def distance_to_arc_nm(self, lat, lon):
        """Approximate surface distance from a point to the locus (dense sample)."""
        best = float("inf")
        for p_lat, p_lon in self.surface_points(step_deg=0.25):
            d = haversine_nm(lat, lon, p_lat, p_lon)
            if d < best:
                best = d
        return best


def build_seventh_arc():
    """Fit the constant-slant-range locus through the documented anchor point."""
    s = satellite_ecef()
    anchor_ecef = geodetic_to_ecef(ARC_ANCHOR[0], ARC_ANCHOR[1], 0.0)
    slant_range_km = math.dist(s, anchor_ecef)
    return Arc(sat_ecef=s, slant_range_km=slant_range_km)


def arc_segment(arc, lon_min=60.0, lon_max=112.0):
    """Sample points along the southern branch inside [lon_min, lon_max]."""
    return arc.surface_points(lon_min=lon_min, lon_max=lon_max, step_deg=0.1)
