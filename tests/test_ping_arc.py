from mh370.geodesy import haversine_nm
from mh370.ping_arc import ARC_ANCHOR, build_seventh_arc


def _scan(arc, step_deg):
    return [arc.point_at_bearing(i * step_deg) for i in range(int(360.0 / step_deg))]


def test_arc_passes_through_documented_anchor():
    arc = build_seventh_arc()
    pts = _scan(arc, 0.05)
    d = min(haversine_nm(la, lo, ARC_ANCHOR[0], ARC_ANCHOR[1]) for la, lo in pts)
    assert d < 2.0


def test_arc_extent_matches_jacc_statement():
    """JACC (Jun 2014): the seventh arc reaches from latitude 20S to 39S."""
    arc = build_seventh_arc()
    pts = _scan(arc, 0.1)

    def crossing_lons(lat_target, lon_band):
        return [lo for la, lo in pts if abs(la - lat_target) < 0.3 and lon_band[0] <= lo <= lon_band[1]]

    north_crossings = crossing_lons(-20.0, (90.0, 115.0))
    south_crossings = crossing_lons(-39.0, (75.0, 95.0))
    assert north_crossings, "arc should cross lat -20 in the eastern Indian Ocean"
    assert south_crossings, "arc should cross lat -39 in the southern Indian Ocean"
