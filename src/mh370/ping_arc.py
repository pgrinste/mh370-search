"""Geometry of the 'seventh arc' from the 2014-03-09 08:19:54 UTC ATOC ping.

The BTO measurement fixes the slant range between the aircraft and the
Inmarsat-3 F1 geostationary satellite (sub-satellite point at 0 deg, 64.5 deg E).
As a first-order approximation we model the constant-range locus as a circle on
the Earth's surface centred on that sub-satellite point, with radius fitted to
pass through the documented anchor point (34.13 S, 93.95 E) - the seventh-arc
longitude at 34.13 S for an aircraft altitude of 20,000 ft.

Validation: JACC stated the arc's total extent reaches from latitude 20 S to
39 S (MarineLink, 20 Jun 2014). See tests/test_ping_arc.py.
"""

from dataclasses import dataclass

from .geodesy import destination_point, haversine_nm

SATELLITE_SUBPOINT = (0.0, 64.5)  # Inmarsat-3 F1 geostationary slot
ARC_ANCHOR = (-34.13, 93.95)      # documented point on the seventh arc


@dataclass(frozen=True)
class Arc:
    center_lat: float
    center_lon: float
    radius_nm: float

    def point_at_bearing(self, bearing_deg):
        return destination_point(
            self.center_lat, self.center_lon, bearing_deg, self.radius_nm
        )

    def distance_to_arc_nm(self, lat, lon):
        """Approximate surface distance from a point to the arc circle."""
        d_center = haversine_nm(lat, lon, self.center_lat, self.center_lon)
        return abs(d_center - self.radius_nm)


def build_seventh_arc():
    radius = haversine_nm(
        SATELLITE_SUBPOINT[0], SATELLITE_SUBPOINT[1], ARC_ANCHOR[0], ARC_ANCHOR[1]
    )
    return Arc(SATELLITE_SUBPOINT[0], SATELLITE_SUBPOINT[1], radius)


def arc_segment(arc, lat_min=-40.0, lat_max=-18.0, n_points=360):
    """Sample points along the part of the circle inside [lat_min, lat_max] (south band)."""
    pts = []
    for i in range(n_points * 4):
        bearing = i * (360.0 / (n_points * 4))
        lat, lon = arc.point_at_bearing(bearing)
        if lat_min <= lat <= lat_max:
            pts.append((lat, lon))
    # de-duplicate by rounding and sort by longitude for a clean polyline
    seen = {}
    for lat, lon in pts:
        key = (round(lat, 3), round(lon, 3))
        seen[key] = (lat, lon)
    return sorted(seen.values(), key=lambda p: p[1])
