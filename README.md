# mh370-search

Where is Flight MH370? A probabilistic search of the southern Indian Ocean built from public data: radar fixes, satellite ATOC timestamps, the seventh-arc geometry, published BFO values, and confirmed debris locations.

[![tests](https://github.com/pgrinste/mh370-search/actions/workflows/ci.yml/badge.svg)](https://github.com/pgrinste/mh370-search/actions/workflows/ci.yml)

## What it does

1. Reconstructs the "seventh arc" — the constant-**slant-range** locus from the final ATOC handshake (2014-03-09 00:19:54 UTC = 08:19:54 MYT) as the intersection of a sphere centred on Inmarsat-3 F1 with the WGS84 ellipsoid, fitted to pass through the documented anchor point (34.13S, 93.95E).
2. Validates that geometry against JACC's published statement that the arc reaches from latitude 20S to 39S (`tests/test_ping_arc.py`), and pins down how much of a ground-circle approximation was error: with both models calibrated on the same anchor they agree to ~10 km over the southern band.
3. Samples ~2,000 candidate flight paths from the turn region (Andaman Sea / S China Sea) to points on the arc, keeping only those whose speed and implied meander time are physically plausible given fuel exhaustion at the arc (~15.6 h total flight time).
4. Scores each path against the **published BFO sequence** (arcs 2–6: 111 → 252 Hz) with a lightweight geometric Doppler model — a soft consistency prior, not a full Holland/DSTG decomposition (`src/mh370/bfo.py`).
5. Drifts floating debris from each grid cell under a simplified regional current field and scores agreement with where debris actually washed up. The field is **calibrated on the Reunion flaperon only**; the Mozambique belly panel (Mar 2017) is held out as validation — its trajectory passes within ~650 nm of the Maputo find location partway through the window, consistent with beaching before discovery.
6. Multiplies path likelihood × BFO agreement × drift agreement into a probability heatmap.

## Result (v1.5)

The top-scoring cell sits at **34S, 94E** — essentially on the ATSB unsearched-area centre (34S 93E) and within ~2° of CSIRO's final most-likely impact estimate (35.6S 92.8E), using only public data and a deliberately simple model:

![probability heatmap](output/mh370_probability_heatmap.png)

Matching those numbers is a sanity check, not an independent discovery — the drift field was calibrated to reach Reunion, so the result carries that prior. The holdout (Mozambique) and BFO checks are what give it any teeth.

## Interactive map

The same result as a browsable web page: probability heatmap, seventh arc with its documented anchor, top candidate tracks (great-circle), debris finds, and the official reference estimates — [output/mh370_interactive_map.html](output/mh370_interactive_map.html) opens in any browser (map tiles load from OSM).

Regenerate it with `pip install -e ".[viz]"` then `python -m mh370.interactive_map`.

## Layout

```
data/curated/   cleaned, documented datasets + data dictionary (committed)
src/mh370/      package: geodesy, ping arc, flight paths, BFO, drift, scoring, interactive map
tests/          unit tests incl. JACC arc-extent and holdout-drift validation
output/         generated heatmap + interactive HTML map
```

## Setup

```bash
python -m venv .venv
.venv\Scripts\pip install -e ".[dev]"        # Windows; use bin/pip on Linux/macOS
python -m mh370.run_analysis --candidates 2000
pytest tests/
```

(Or the old-fashioned way: `pip install -r requirements.txt` and set `PYTHONPATH=src`.)

## Known simplifications (v1.5)

- **Flight kinematics** are a straight track plus loiter window, not turn-by-turn dynamics with fuel burn. The BFO filter is therefore soft (σ = 60 Hz, residuals re-centred per handshake): our straight-line paths over-predict receding Doppler by several hundred Hz near the arc because real tracks had a stronger westward component.
- **Drift** uses a two-region constant current field calibrated on Reunion only. A real model would use reanalysis currents (CMEMS GLORYS12V1) with per-debris-type wind leeway — that is the v2 upgrade path, and it is also what would let us distinguish 30S from 35S source cells.
- **Satellite** treated as fixed at its geostationary slot; orbital drift adds ~±10 Hz to BFO (inside the scoring sigma) and a few km to arc geometry (absorbed by the anchor fit).
- The path likelihood is essentially uniform along the arc apart from the speed/meander gates — there is no fuel-exhaustion constraint beyond "crossed the arc before running dry".

Full source list, confidence levels, and the MYT/UTC time-base notes: [data/curated/data_dictionary.md](data/curated/data_dictionary.md).
