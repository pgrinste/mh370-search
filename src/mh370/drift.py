"""First-order ocean drift of floating debris.

A two-region constant current field approximating the Agulhas extension along
Africa's east coast plus monsoon-driven surface drift in the central Indian
Ocean (see data/curated/data_dictionary.md).

Calibration design (to avoid circularity): the field is fitted so that mean
drift from the southern seventh-arc region reaches near the Reunion flaperon
find (Jul 2015, ~507 days of advection). The Mozambique belly panel (Mar
2017, ~1100 days) and all other finds are held out as validation: the field's
trajectory passes within a few hundred nm of the Maputo find location partway
through the window (~day 740), consistent with the panel beaching there before
its March 2017 discovery. See tests/test_drift.py.

A real model would use reanalysis currents (e.g. CMEMS GLORYS12V1) with
per-debris-type wind leeway; that is the documented v2 upgrade path.

Units: knots, degrees of latitude/longitude per step handled in nm.
"""

from .geodesy import destination_point


def current_at(lat, lon):
    """Return (u_ks eastward, v_ks northward) at a location."""
    if lon <= 60.0 and lat < -15.0:
        # approaching Africa / Agulhas extension: faster northwest drift
        return (-0.28, 0.12)
    # open central Indian Ocean: slow west-northwest surface drift
    return (-0.13, 0.045)


def advect(lat, lon, days, step_h=6.0):
    """Advect a floating object for `days` under the regional current field."""
    steps = int(days * 24.0 / step_h)
    for _ in range(steps):
        u, v = current_at(lat, lon)
        dist_nm = (u**2 + v**2) ** 0.5 * step_h
        if dist_nm < 1e-6:
            break
        bearing = (90.0 - _bearing_from_uv(u, v)) % 360.0
        lat, lon = destination_point(lat, lon, bearing, dist_nm)
    return lat, lon


def min_distance_to(lat, lon, days, target_lat, target_lon, step_h=6.0):
    """Minimum distance (nm) from the drift trajectory to a target over `days`.

    Returns (min_distance_nm, day_of_closest_approach). Beach finds are not
    point-in-time targets - debris can sit on a shore for months before being
    reported - so the closest approach anywhere in the window is the fair
    comparison.
    """
    from .geodesy import haversine_nm

    best_d, best_day = float("inf"), 0
    for day in range(1, int(days) + 1):
        lat, lon = advect(lat, lon, 1.0, step_h=step_h)
        d = haversine_nm(lat, lon, target_lat, target_lon)
        if d < best_d:
            best_d, best_day = d, day
    return best_d, best_day


def _bearing_from_uv(u_ks, v_ks):
    """Math angle (deg, CCW from east) of a velocity vector."""
    import math

    return math.degrees(math.atan2(v_ks, u_ks))
