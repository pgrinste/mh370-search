from mh370.drift import advect, min_distance_to
from mh370.geodesy import haversine_nm


def test_flaperon_calibration_reaches_near_reunion():
    # Calibration target: from the southern seventh-arc region, 507 days of
    # drift (Mar 2014 -> Jul 2015) must reach near the Reunion flaperon.
    end_lat, end_lon = advect(-35.0, 93.0, 507)
    d = haversine_nm(end_lat, end_lon, -21.08, 55.0)
    assert d < 700.0


def test_mozambique_heldout_validation():
    # Held out: the current field was NOT fitted to this find (Mar 2017).
    # Beach finds are not point-in-time targets - check closest approach over
    # the full ~1100-day window instead.
    d, day = min_distance_to(-35.0, 93.0, 1102, -26.0, 32.6)
    assert d < 1000.0
    # and the approach happens well before the find date (beaching lead time)
    assert day < 1000


def test_drift_moves_northwest_in_open_ocean():
    lat0, lon0 = advect(-30.0, 90.0, 30)
    assert lat0 > -30.0      # northward component
    assert lon0 < 90.0       # westward component
