from mh370.drift import advect
from mh370.geodesy import haversine_nm


def test_flaperon_drift_reaches_near_reunion():
    # from the southern seventh-arc region, 476 days of drift (Mar 2014 -> Jul 2015)
    end_lat, end_lon = advect(-35.0, 93.0, 476)
    d = haversine_nm(end_lat, end_lon, -21.08, 55.0)
    assert d < 700.0


def test_drift_moves_northwest_in_open_ocean():
    lat0, lon0 = advect(-30.0, 90.0, 30)
    assert lat0 > -30.0      # northward component
    assert lon0 < 90.0       # westward component
