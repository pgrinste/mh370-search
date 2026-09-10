"""Spherical geometry helpers (WGS84, nautical miles)."""

import math

EARTH_RADIUS_NM = 3440.065


def _to_rad(x):
    return math.radians(x)


def haversine_nm(lat1, lon1, lat2, lon2):
    """Great-circle distance in nautical miles."""
    p1, l1, p2, l2 = map(_to_rad, (lat1, lon1, lat2, lon2))
    d_p = (p2 - p1) / 2.0
    d_l = (l2 - l1) / 2.0
    a = math.sin(d_p) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(d_l) ** 2
    return 2.0 * EARTH_RADIUS_NM * math.asin(min(1.0, math.sqrt(a)))


def initial_bearing_deg(lat1, lon1, lat2, lon2):
    """Initial great-circle bearing from point 1 to point 2, degrees true."""
    p1, l1, p2, l2 = map(_to_rad, (lat1, lon1, lat2, lon2))
    d_l = l2 - l1
    x = math.sin(d_l) * math.cos(p2)
    y = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(d_l)
    return (math.degrees(math.atan2(x, y)) + 360.0) % 360.0


def destination_point(lat, lon, bearing_deg, dist_nm):
    """Point reached by travelling dist_nm from (lat, lon) along a great circle at initial bearing."""
    p1 = _to_rad(lat)
    l1 = _to_rad(lon)
    th = _to_rad(bearing_deg)
    dr = dist_nm / EARTH_RADIUS_NM
    p2 = math.asin(
        math.sin(p1) * math.cos(dr) + math.cos(p1) * math.sin(dr) * math.cos(th)
    )
    l2 = l1 + math.atan2(
        math.sin(th) * math.sin(dr) * math.cos(p1),
        math.cos(dr) - math.sin(p1) * math.sin(p2),
    )
    return math.degrees(p2), (math.degrees(l2) + 540.0) % 360.0 - 180.0
