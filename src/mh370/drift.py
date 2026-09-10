"""First-order ocean drift of floating debris.

A constant regional current field, chosen so that mean drift from the southern
seventh-arc region reaches both confirmed find locations (Reunion flaperon,
Jul 2015; Mozambique belly panel, Mar 2017) within a few hundred nm over their
respective advection periods. The real flow is dominated by the Agulhas
extension along Africa's east coast plus monsoon-driven surface drift; this
two-region constant field is a documented v1 approximation (see
data/curated/data_dictionary.md).

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


def _bearing_from_uv(u_ks, v_ks):
    """Math angle (deg, CCW from east) of a velocity vector."""
    import math

    return math.degrees(math.atan2(v_ks, u_ks))
