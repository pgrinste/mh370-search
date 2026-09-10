"""Monte Carlo candidate flight paths from the turn region to the seventh arc.

Simplified model (see data/curated/data_dictionary.md for sources):
- The aircraft turned back over the South China Sea / Andaman Sea region and
  flew south-southwest toward the southern Indian Ocean.
- It crossed the seventh arc at or just before fuel exhaustion, i.e. around
  the time of the seventh ping (2014-03-09 08:19:54 UTC). JACC considered the
  aircraft to have exhausted its fuel on reaching the arc and descending.
- A straight track from the turn region is far shorter than the distance the
  aircraft could fly in the ~13 h between last military contact and the
  seventh ping, so public analyses invoke a meandering ('zig-zag') track. We
  approximate that by allowing each candidate to cruise at a plausible speed
  and then loiter/meander near its arc crossing until the ping time.

Each candidate is accepted only if:
- required ground speed is within [375, 500] kt (the range used in Inmarsat's
  own track modelling), and
- the implied meander/loiter duration is within [1 h, 8 h].
"""

import random
from dataclasses import dataclass

from .geodesy import haversine_nm
from .ping_arc import Arc

SEVENTH_PING_SECONDS = (9 * 24 + 8) * 3600 + 19 * 60 + 54 - (8 * 24) * 3600
# seconds since 2014-03-08T00:00Z for the ping on 2014-03-09T08:19:54Z

START_BOX = {
    "lat": (-3.0, 5.0),     # turn region: Andaman Sea / S China Sea (deg)
    "lon": (97.0, 102.0),   # deg E
}
TURN_START_RANGE_H = (17.5, 20.0)   # plausible start of southbound track, hours UTC on Mar 8
SPEED_RANGE_KT = (375.0, 500.0)     # Inmarsat track-modelling range
MEANDER_RANGE_H = (1.0, 8.0)        # implied loiter/zig-zag duration


@dataclass(frozen=True)
class Candidate:
    start_lat: float
    start_lon: float
    impact_lat: float
    impact_lon: float
    speed_kt: float
    meander_h: float


def _arc_crossing_at_latitude(arc, lat):
    """Point on the arc circle at a given latitude (eastern branch)."""
    # scan bearings; keep points near the requested latitude in the east band
    best = None
    for i in range(720):
        bearing = i * 0.5
        p_lat, p_lon = arc.point_at_bearing(bearing)
        if abs(p_lat - lat) < 0.05 and p_lon > 60.0:
            best = (p_lat, p_lon)
            break
    return best


def simulate_candidates(arc, n=2000, seed=7):
    rng = random.Random(seed)
    out = []
    for _ in range(n * 8):
        if len(out) >= n:
            break
        start_lat = rng.uniform(*START_BOX["lat"])
        start_lon = rng.uniform(*START_BOX["lon"])
        # target latitude on the arc, weighted toward the ATSB high-probability band
        r = rng.random()
        if r < 0.7:
            target_lat = rng.uniform(-37.0, -32.5)   # ATSB FPR band
        else:
            target_lat = rng.uniform(-40.0, -20.0)   # full JACC extent
        crossing = _arc_crossing_at_latitude(arc, round(target_lat, 1))
        if crossing is None:
            continue
        impact_lat, impact_lon = crossing

        dist_nm = haversine_nm(start_lat, start_lon, impact_lat, impact_lon)
        speed = rng.uniform(*SPEED_RANGE_KT)
        cruise_h = dist_nm / speed
        turn_start_h = rng.uniform(*TURN_START_RANGE_H)
        meander_h = (SEVENTH_PING_SECONDS / 3600.0) - turn_start_h - cruise_h
        if not (MEANDER_RANGE_H[0] <= meander_h <= MEANDER_RANGE_H[1]):
            continue
        out.append(
            Candidate(start_lat, start_lon, impact_lat, impact_lon, speed, meander_h)
        )
    return out
