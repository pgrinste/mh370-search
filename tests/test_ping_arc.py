import math

from mh370.geodesy import destination_point, geodetic_to_ecef, haversine_nm
from mh370.ping_arc import ARC_ANCHOR, build_seventh_arc


def _scan(arc):
    return arc.surface_points(lon_min=60.0, lon_max=112.0, step_deg=0.05)


def test_arc_passes_through_documented_anchor():
    arc = build_seventh_arc()
    pts = _scan(arc)
    d = min(haversine_nm(la, lo, ARC_ANCHOR[0], ARC_ANCHOR[1]) for la, lo in pts)
    assert d < 2.0


def test_slant_range_is_constant_along_locus():
    """Defining property: every locus point sits at the same slant range from the satellite."""
    arc = build_seventh_arc()
    ranges = []
    for la, lo in _scan(arc)[::20]:
        p = geodetic_to_ecef(la, lo, 0.0)
        ranges.append(math.dist(p, arc.sat_ecef))
    spread_km = max(ranges) - min(ranges)
    assert spread_km < 1.0


def test_arc_extent_matches_jacc_statement():
    """JACC (Jun 2014): the seventh arc reaches from latitude 20S to 39S."""
    arc = build_seventh_arc()
    pts = _scan(arc)

    def crossing_lons(lat_target, lon_band):
        return [lo for la, lo in pts if abs(la - lat_target) < 0.3 and lon_band[0] <= lo <= lon_band[1]]

    north_crossings = crossing_lons(-20.0, (90.0, 115.0))
    south_crossings = crossing_lons(-39.0, (75.0, 95.0))
    assert north_crossings, "arc should cross lat -20 in the eastern Indian Ocean"
    assert south_crossings, "arc should cross lat -39 in the southern Indian Ocean"


def test_ground_circle_agreement():
    """With both models calibrated on the same anchor they agree to ~15 km over the southern band.

    The difference is purely ellipsoidal (a constant-slant-range locus is an
    exact ground circle only for a spherical Earth), so this pins down how much
    of the old v1 result was approximation error.
    """
    arc = build_seventh_arc()
    radius_nm = haversine_nm(0.0, 64.5, ARC_ANCHOR[0], ARC_ANCHOR[1])

    def old_lon_at_lat(lat_target):
        for i in range(7200):
            la, lo = destination_point(0.0, 64.5, i * 0.05, radius_nm)
            if abs(la - lat_target) < 0.01 and lo > 60.0:
                return lo
        return None

    for lat in (-32.0, -34.13, -36.0, -38.0):
        old_lon = old_lon_at_lat(lat)
        new_lons = [lo for la, lo in _scan(arc) if abs(la - lat) < 0.05]
        assert old_lon is not None and new_lons
        d_km = haversine_nm(lat, old_lon, lat, new_lons[0]) * 1.852
        assert d_km < 15.0
