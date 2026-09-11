"""Monte Carlo candidate flight paths from the turn region to the seventh arc.

Simplified model (see data/curated/data_dictionary.md for sources):
- The aircraft turned back over the South China Sea / Andaman Sea region and
  flew south-southwest toward the southern Indian Ocean.
- It crossed the seventh arc at or just before fuel exhaustion, i.e. around
  the time of the final handshake (2014-03-09 00:19:54 UTC = 08:19:54 MYT).
  JACC considered the aircraft to have exhausted its fuel on reaching the arc
  and descending. Total flight time from departure is ~15.6 h, consistent with
  a heavily fuelled B777-200ER at cruise burn rates.
- A straight track from the turn region is far shorter than the distance the
  aircraft could fly in the ~6 h between last military radar contact (~18:30 UTC
  Mar 8) and the final handshake, so public analyses invoke a meandering
  ('zig-zag') track. We approximate that by allowing each candidate to cruise at
  a plausible speed and then loiter/meander near its arc crossing until the
  handshake time.

Each candidate is accepted only if:
- required ground speed is within [375, 500] kt (the range used in Inmarsat's
  own track modelling), and
- the implied meander/loiter duration is within [1 h, 8 h].
"""

import random
from dataclasses import dataclass

from .geodesy import destination_point, haversine_nm, initial_bearing_deg
from .ping_arc import Arc

# Seconds since 2014-03-08T00:00Z for the final handshake on 2014-03-09T00:19:54Z.
# (Public reporting often quotes "08:19", which is Malaysian time, UTC+8.)
SEVENTH_PING_SECONDS = 24 * 3600 + 19 * 60 + 54

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

    @property
    def cruise_hours(self):
        dist_nm = haversine_nm(
            self.start_lat, self.start_lon, self.impact_lat, self.impact_lon
        )
        return dist_nm / self.speed_kt

    @property
    def final_bearing_deg(self):
        """Initial great-circle bearing of the start -> impact leg."""
        return initial_bearing_deg(
            self.start_lat, self.start_lon, self.impact_lat, self.impact_lon
        )

    @property
    def turn_start_h(self):
        """Implied start of the southbound track (hours UTC on Mar 8)."""
        return SEVENTH_PING_SECONDS / 3600.0 - self.cruise_hours - self.meander_h


def _arc_crossing_at_latitude(points, lat):
    """Point on the arc sample nearest a requested latitude."""
    best = None
    for p_lat, p_lon in points:
        if abs(p_lat - lat) < 0.05:
            return (p_lat, p_lon)
        if best is None or abs(p_lat - lat) < abs(best[0] - lat):
            best = (p_lat, p_lon)
    return best


def simulate_candidates(arc, n=2000, seed=7):
    arc_points = arc.surface_points(lon_min=60.0, lon_max=112.0, step_deg=0.05)
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
        crossing = _arc_crossing_at_latitude(arc_points, round(target_lat, 1))
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


def position_at(candidate, t_hours):
    """Aircraft state at time t (hours UTC on Mar 8).

    Returns (lat, lon, bearing_deg or None). During the loiter window the
    aircraft is treated as stationary at its arc crossing (bearing None) -
    a documented simplification for Doppler scoring.
    """
    turn_start_h = candidate.turn_start_h
    if t_hours < turn_start_h:
        return candidate.start_lat, candidate.start_lon, None
    cruise_end = turn_start_h + candidate.cruise_hours
    bearing = candidate.final_bearing_deg
    if t_hours >= cruise_end:
        return candidate.impact_lat, candidate.impact_lon, None
    dist_nm = haversine_nm(
        candidate.start_lat, candidate.start_lon,
        candidate.impact_lat, candidate.impact_lon,
    )
    frac = (t_hours - turn_start_h) / candidate.cruise_hours
    lat, lon = destination_point(
        candidate.start_lat, candidate.start_lon, bearing, dist_nm * frac
    )
    return lat, lon, bearing
