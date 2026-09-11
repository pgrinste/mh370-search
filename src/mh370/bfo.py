"""BFO (burst frequency offset) consistency scoring for candidate paths.

The Inmarsat ground station recorded a burst frequency offset at each
handshake: the Doppler shift of the signal relative to its expected value.
After removing a ~150 Hz bias (satellite oscillator drift plus AES
compensation; Holland 2017, arXiv:1702.02432), the residual tracks the
aircraft's line-of-sight velocity relative to the satellite:

    f_net [Hz] ~= v_los [m/s] * f_c / c      (f_c ~ 1.54 GHz L-band uplink)

The published sequence for the hourly handshakes rises monotonically,
111 -> 141 -> 168 -> 204 -> 252 Hz (arcs 2-6), consistent with an aircraft
turning south: its line-of-sight velocity rotates from approaching to
receding as it tracks toward the southern arc.

Scoring here is deliberately lightweight: for each candidate we compute the
geometric net Doppler at each handshake time and take a product of Gaussians
against the observed BFO. Residuals are re-centred per handshake (common-mode
model error removed) and sigma is kept wide (60 Hz), so this acts as a soft
consistency prior rather than a hard filter - it is not a full Holland/DSTG
five-component decomposition, which needs turn-by-turn dynamics.

Documented approximations:
- Satellite fixed at its geostationary slot; orbital drift adds ~10 Hz,
  inside the scoring sigma.
- Straight-line cruise then stationary loiter (no turn-by-turn dynamics).
- f_c = 1.54 GHz representative L-band uplink frequency.
- Handshake times for arcs 2-6 carry +/-30 min uncertainty (hourly GES
  schedule; see data/curated/mh370_satcom.csv for provenance).
"""

import math

from .flight_path import position_at
from .geodesy import geodetic_to_ecef, surface_tangent_ecef
from .ping_arc import satellite_ecef

C_MPS = 299792458.0     # speed of light, m/s
F_C_HZ = 1.54e9         # representative L-band uplink frequency
BFO_BIAS_HZ = 150.0     # Holland 2017 (Ashton et al. 2014 used 152.5)
SIGMA_HZ = 60.0         # soft filter: ~half the systematic Doppler offset of our
                        # straight-line kinematics (DSTG's own in-flight sigma was
                        # 4.3-7 Hz, but that assumes turn-by-turn dynamics)

# Observed BFO at the hourly handshakes (arcs 2-6). Times are hours UTC on
# Mar 8 (Mar 9 values expressed as >24 h), taken as published.
#
# Reading note: if the early entries were actually Malaysian time, arcs 2-5
# would fall in the pre-turn phase and contribute only a common factor to all
# candidates; under the UTC reading they fall during the southbound cruise and
# discriminate between paths. The ranking is dominated by arc-6 either way.
BFO_OBSERVED = [
    # (t_hours, bfo_hz, label)
    (19.68, 111.0, "arc2"),
    (20.68, 141.0, "arc3"),
    (21.69, 168.0, "arc4"),
    (22.69, 204.0, "arc5"),
    (24.17, 252.0, "arc6"),   # last scheduled handshake before the final log-on
]

_SAT = satellite_ecef()


def net_doppler_hz(lat, lon, bearing_deg, speed_kt):
    """Geometric net Doppler (Hz) for an aircraft at (lat, lon).

    Positive = receding from the satellite. bearing None => stationary
    (loiter), so only the ~0 Hz geometric term remains.
    """
    if bearing_deg is None or speed_kt <= 0:
        return 0.0
    p = geodetic_to_ecef(lat, lon, 0.0)
    r = math.dist(_SAT, p)
    u = tuple((p[i] - _SAT[i]) / r for i in range(3))   # satellite -> aircraft
    v_dir = surface_tangent_ecef(lat, lon, bearing_deg)
    v_los_ms = (v_dir[0] * u[0] + v_dir[1] * u[1] + v_dir[2] * u[2]) \
        * speed_kt * 0.514444
    return v_los_ms * F_C_HZ / C_MPS


def bfo_residuals(candidate):
    """(label, residual_hz) at each observed handshake for this candidate."""
    out = []
    for t_h, obs_hz, label in BFO_OBSERVED:
        lat, lon, bearing = position_at(candidate, t_h)
        pred_net = net_doppler_hz(lat, lon, bearing, candidate.speed_kt)
        out.append((label, (obs_hz - BFO_BIAS_HZ) - pred_net))
    return out


def bfo_log_weight(candidate):
    """Log-likelihood of the observed BFO sequence under this path.

    Computed in log space because the product of five Gaussians can underflow
    double precision for grossly inconsistent paths (which is the point).
    """
    lw = 0.0
    for _, resid in bfo_residuals(candidate):
        lw -= 0.5 * (resid / SIGMA_HZ) ** 2
    return lw


def bfo_weights(candidates):
    """Relative BFO likelihoods, rescaled so the best candidate has weight 1.

    Residuals are re-centred at each handshake time (common-mode model error
    removed): our straight-line kinematics carries a systematic Doppler offset
    of several hundred Hz near the arc, so only the *differential* agreement
    between paths is scored. Handshake times where all candidates predict the
    same value (e.g. the loiter phase) contribute nothing by construction.
    """
    if not candidates:
        return []
    res = [bfo_residuals(c) for c in candidates]
    n_t = len(BFO_OBSERVED)
    means = [
        sum(r[k][1] for r in res) / len(res) for k in range(n_t)
    ]
    log_ws = []
    for r in res:
        lw = 0.0
        for k in range(n_t):
            d = r[k][1] - means[k]
            lw -= 0.5 * (d / SIGMA_HZ) ** 2
        log_ws.append(lw)
    m = max(log_ws)
    return [math.exp(lw - m) for lw in log_ws]
