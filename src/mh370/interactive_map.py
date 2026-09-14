"""Interactive HTML map of the search results (folium).

Usage:  python -m mh370.interactive_map [--candidates 2000] [--top-tracks 15]
Output: output/mh370_interactive_map.html

Same pipeline as run_analysis, rendered for a browser instead of matplotlib:
probability heatmap, seventh arc + documented anchor, top candidate tracks
(great-circle interpolation), debris finds and the official reference estimates.
"""

import argparse
import os

# Pin BLAS threads BEFORE numpy import (see run_analysis docstring).
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import folium
from folium.plugins import HeatMap

from .bfo import bfo_weights
from .flight_path import simulate_candidates
from .geodesy import destination_point, haversine_nm, initial_bearing_deg
from .ping_arc import ARC_ANCHOR, arc_segment, build_seventh_arc
from .run_analysis import load_debris
from .scoring import drift_score, score_grid

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# rank -> colour for candidate tracks (top-ranked first)
TRACK_COLORS = ["#d7301f", "#e4572e", "#f38b29", "#fab833", "#ffd94a",
                "#ffea8c", "#d9ead3", "#b6d7a8", "#92c5de", "#6ca0cd"]


def track_points(cand, n=36):
    """Great-circle interpolation start -> impact (folium wants [lat, lon])."""
    d = haversine_nm(cand.start_lat, cand.start_lon, cand.impact_lat, cand.impact_lon)
    brg = initial_bearing_deg(cand.start_lat, cand.start_lon,
                              cand.impact_lat, cand.impact_lon)
    pts = [[cand.start_lat, cand.start_lon]]
    for i in range(1, n + 1):
        f = i / n
        lat, lon = destination_point(cand.start_lat, cand.start_lon, brg, d * f)
        pts.append([lat, lon])
    return pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=int, default=2000)
    ap.add_argument("--top-tracks", type=int, default=15)
    args = ap.parse_args()

    arc = build_seventh_arc()
    candidates = simulate_candidates(arc, n=args.candidates)
    weights = bfo_weights(candidates)
    debris = load_debris()
    print(f"{len(candidates)} candidate paths, {len(debris)} debris items")

    # same grid as run_analysis (1 deg resolution)
    lats = [round(-40.0 + i * 1.0, 2) for i in range(25)]   # -40 .. -16
    lons = [round(78.0 + j * 1.0, 2) for j in range(31)]    # 78 .. 108
    grid, best = score_grid(candidates, weights, debris, lats, lons)

    best_cell = max(((i, j) for i in range(len(lats)) for j in range(len(lons))),
                    key=lambda ij: grid[ij[0]][ij[1]])
    top_i, top_j = best_cell
    print(f"top cell: {lats[top_i]:.2f}, {lons[top_j]:.2f}")

    m = folium.Map(location=[-32.0, 95.0], zoom_start=4)

    # --- probability heatmap -------------------------------------------
    heat = [[lats[i], lons[j], grid[i][j] / best]
            for i in range(len(lats)) for j in range(len(lons)) if grid[i][j] > 0]
    HeatMap(heat, radius=14, blur=10, max_value=1.0, min_opacity=0.25).add_to(m)

    # --- seventh arc + documented anchor --------------------------------
    seg = arc_segment(arc)
    if seg:
        folium.PolyLine([[p[0], p[1]] for p in seg], color="#d7301f", weight=2.5,
                        dash_array="8 6", tooltip="seventh arc (constant slant-range locus)"
                        ).add_to(m)
    a_lat, a_lon = ARC_ANCHOR
    folium.CircleMarker([a_lat, a_lon], radius=6, color="#0aa", fill=True,
                        fill_color="#0aa",
                        tooltip="documented arc anchor (34.13S 93.95E)").add_to(m)

    # --- top candidate tracks -------------------------------------------
    ranked = sorted(((w * drift_score(debris, c.impact_lat, c.impact_lon), c)
                     for c, w in zip(candidates, weights) if w > 0.01),
                    reverse=True)[:args.top_tracks]
    for k, (score, cand) in enumerate(ranked):
        color = TRACK_COLORS[min(k, len(TRACK_COLORS) - 1)]
        folium.PolyLine(track_points(cand), color=color, weight=2.0 if k < 3 else 1.2,
                        opacity=0.9 if k < 3 else 0.65,
                        tooltip=(f"candidate track #{k + 1} — {cand.speed_kt:.0f} kt, "
                                 f"{cand.meander_h:.1f} h meander, impact "
                                 f"{abs(cand.impact_lat):.2f}S {cand.impact_lon:.2f}E")).add_to(m)

    # --- debris finds -----------------------------------------------------
    for item in debris:
        if not (-42.0 <= item["find_lat"] <= -12.0 and 48.0 <= item["find_lon"] <= 112.0):
            continue  # far-field finds stay off this map (as in the static plot)
        folium.Marker([item["find_lat"], item["find_lon"]], icon=folium.Icon(color="green",
                                                                              icon="triangle"),
                      tooltip=f"{item['item']} ({item['status']})").add_to(m)

    # --- official reference estimates --------------------------------------
    for lat, lon, label, color in (
        (-34.0, 93.0, "ATSB unsearched-area centre (34S 93E)", "#555"),
        (-35.6, 92.8, "CSIRO most-likely impact (35.6S 92.8E)", "#b8860b"),
    ):
        folium.CircleMarker([lat, lon], radius=7, color=color, fill=True,
                            fill_color=color, tooltip=label).add_to(m)

    # --- top cell -----------------------------------------------------------
    folium.Marker([lats[top_i], lons[top_j]], icon=folium.Icon(color="red", icon="star"),
                  tooltip=f"top-scoring cell ({lats[top_i]:.0f}S {lons[top_j]:.0f}E), "
                          f"relative probability 1.0").add_to(m)

    # --- legend ---------------------------------------------------------------
    m.get_root().html.add_child(folium.Element("""
    <div style="position:fixed;bottom:24px;left:14px;z-index:1000;background:white;
                border-radius:8px;padding:10px 12px;font:12px/1.5 sans-serif;
                box-shadow:0 1px 6px rgba(0,0,0,.3)">
      <b>MH370 search — interactive</b><br>
      heat = relative probability (path likelihood &times; BFO agreement &times; drift)<br>
      red dashed = seventh arc &nbsp;&middot;&nbsp; lines = top candidate tracks<br>
      green ^ = debris finds &nbsp;&middot;&nbsp; star = top cell (v1.5 model)
    </div>"""))

    os.makedirs(os.path.join(REPO_ROOT, "output"), exist_ok=True)
    out = os.path.join(REPO_ROOT, "output", "mh370_interactive_map.html")
    m.save(out)
    print(f"wrote {out} ({os.path.getsize(out) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
