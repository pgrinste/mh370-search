"""Grid scoring: combine path likelihood with debris-drift agreement."""

import math

from .drift import advect
from .geodesy import haversine_nm

STATUS_WEIGHT = {"confirmed": 1.0, "probable": 0.5, "unconfirmed": 0.2}


def _gauss(x, sigma):
    return math.exp(-0.5 * (x / sigma) ** 2)


def path_density(candidates, weights, lat, lon, kernel_nm=60.0):
    """Weighted kernel density of candidate impact points at a grid cell.

    `weights` are per-candidate likelihoods (e.g. BFO agreement); candidates
    with zero weight contribute nothing.
    """
    total = 0.0
    wsum = 0.0
    for c, w in zip(candidates, weights):
        if w <= 0.0:
            continue
        d = haversine_nm(lat, lon, c.impact_lat, c.impact_lon)
        total += w * _gauss(d, kernel_nm)
        wsum += w
    return total / wsum if wsum > 0 else 0.0


def drift_score(debris_items, lat, lon):
    """Weighted agreement between advected debris and actual find locations."""
    scores = []
    for item in debris_items:
        w = STATUS_WEIGHT.get(item["status"], 0.2)
        end_lat, end_lon = advect(lat, lon, item["drift_days"])
        d = haversine_nm(end_lat, end_lon, item["find_lat"], item["find_lon"])
        scores.append((w, _gauss(d, sigma=450.0)))
    if not scores:
        return 1.0
    num = sum(w * s for w, s in scores)
    den = sum(w for w, _ in scores)
    return num / den


def score_grid(candidates, weights, debris_items, lats, lons):
    """Return (grid of total scores, max score)."""
    grid = []
    best = 0.0
    for lat in lats:
        row = []
        for lon in lons:
            p = path_density(candidates, weights, lat, lon)
            d = drift_score(debris_items, lat, lon)
            v = p * d
            row.append(v)
            best = max(best, v)
        grid.append(row)
    return grid, best
