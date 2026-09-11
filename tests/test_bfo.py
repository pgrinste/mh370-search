from mh370.bfo import bfo_weights, net_doppler_hz
from mh370.flight_path import simulate_candidates
from mh370.ping_arc import build_seventh_arc


def test_southbound_path_recedes_from_satellite():
    # At a southern-arc position, flying south-southwest increases range to the satellite.
    f = net_doppler_hz(-35.0, 92.0, 195.0, 450.0)
    assert f > 0


def test_westbound_path_along_arc_approaches():
    # Flying due west along the southern arc closes range (satellite is north-west of the point).
    f = net_doppler_hz(-35.0, 92.0, 270.0, 450.0)
    assert f < 0


def test_stationary_loiter_has_zero_geometric_doppler():
    assert net_doppler_hz(-35.0, 92.0, None, 450.0) == 0.0


def test_bfo_weights_vary_across_candidates():
    arc = build_seventh_arc()
    cands = simulate_candidates(arc, n=200, seed=7)
    ws = bfo_weights(cands)
    assert all(0.0 <= w <= 1.0 + 1e-9 for w in ws)
    # the filter must discriminate: several distinct positive weight levels
    positives = sorted({round(w, 6) for w in ws if w > 0.0})
    assert len(positives) >= 3
    assert max(ws) == 1.0
