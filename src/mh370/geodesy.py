"""Spherical and WGS84 ellipsoidal geometry helpers (nautical miles / km)."""

import math

EARTH_RADIUS_NM = 3440.065

# WGS84 ellipsoid
WGS84_A_KM = 6378.137          # semi-major axis, km
WGS84_F = 1.0 / 298.257223563  # flattening
_WGS84_E2 = WGS84_F * (2.0 - WGS84_F)


def geodetic_to_ecef(lat_deg, lon_deg, h_km=0.0):
    """Geodetic (lat, lon, height above ellipsoid) -> ECEF km."""
    p = math.radians(lat_deg)
    l = math.radians(lon_deg)
    n = WGS84_A_KM / math.sqrt(1.0 - _WGS84_E2 * math.sin(p) ** 2)
    return (
        (n + h_km) * math.cos(p) * math.cos(l),
        (n + h_km) * math.cos(p) * math.sin(l),
        (n * (1.0 - _WGS84_E2) + h_km) * math.sin(p),
    )


def ecef_to_geodetic(x, y, z):
    """ECEF km -> (lat_deg, lon_deg). Iterative solution, good to ~1e-6 deg."""
    lat = math.atan2(z, math.hypot(x, y) * (1.0 - _WGS84_E2))
    for _ in range(12):
        n = WGS84_A_KM / math.sqrt(1.0 - _WGS84_E2 * math.sin(lat) ** 2)
        h = math.hypot(x, y) / math.cos(lat) - n
        lat = math.atan2(z + _WGS84_E2 * n * math.sin(lat), math.hypot(x, y))
    return math.degrees(lat), math.degrees(math.atan2(y, x))


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


def surface_tangent_ecef(lat_deg, lon_deg, bearing_deg):
    """Unit vector in ECEF of the local horizontal direction at true bearing (deg)."""
    p = math.radians(lat_deg)
    l = math.radians(lon_deg)
    b = math.radians(bearing_deg)
    east = (-math.sin(l), math.cos(l), 0.0)
    north = (
        -math.sin(p) * math.cos(l),
        -math.sin(p) * math.sin(l),
        math.cos(p),
    )
    return tuple(
        east[i] * math.sin(b) + north[i] * math.cos(b) for i in range(3)
    )
