"""Run the full MH370 analysis: candidates -> grid scoring -> heatmap.

Usage:  python -m mh370.run_analysis [--candidates 2000]
Output: output/mh370_probability_heatmap.png + top cells printed to stdout.

Note on hardware: sustained multithreaded load is a known failure mode on the
Intel i9-14900K, so we pin BLAS threads (set OMP_NUM_THREADS before numpy is
imported). Override with your own value if you know your chip is stable.
"""

import argparse
import csv
import os

# Pin BLAS threads BEFORE numpy import (see module docstring).
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .bfo import bfo_weights
from .drift import min_distance_to
from .flight_path import simulate_candidates
from .ping_arc import ARC_ANCHOR, arc_segment, build_seventh_arc
from .scoring import score_grid

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(REPO_ROOT, "data", "curated")


def load_debris():
    items = []
    with open(os.path.join(DATA, "mh370_debris.csv")) as f:
        for row in csv.DictReader(f):
            parts = row["found_date"][:10].split("-")
            # month-only dates use mid-month; days from 2014-03-09 (disappearance)
            y, m = int(parts[0]), int(parts[1])
            d = int(parts[2]) if len(parts) > 2 else 15
            drift_days = _days_since_2014_03_09(y, m, d)
            if drift_days <= 0:
                continue
            items.append(
                {
                    "item": row["item"],
                    "status": row["status"],
                    "find_lat": float(row["lat_approx"]),
                    "find_lon": float(row["lon_approx"]),
                    "drift_days": drift_days,
                }
            )
    return items


def _days_since_2014_03_09(y, m, d):
    import datetime

    delta = datetime.date(y, m, d) - datetime.date(2014, 3, 9)
    return delta.days


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=int, default=2000)
    args = ap.parse_args()

    arc = build_seventh_arc()
    candidates = simulate_candidates(arc, n=args.candidates)
    print(f"accepted {len(candidates)} candidate paths")

    # BFO consistency weights (relative likelihood vs observed handshake Doppler)
    weights = bfo_weights(candidates)
    n_surviving = sum(1 for w in weights if w > 0.01)   # >1% of best
    print(f"BFO scoring: {n_surviving} candidates retain >1% relative weight")

    debris = load_debris()
    print(f"{len(debris)} debris items with drift periods")

    # grid over the southern Indian Ocean search region (v1: 1 deg resolution)
    lats = [round(-40.0 + i * 1.0, 2) for i in range(25)]   # -40 .. -16
    lons = [round(78.0 + j * 1.0, 2) for j in range(31)]    # 78 .. 108
    grid, best = score_grid(candidates, weights, debris, lats, lons)

    top = []
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            if grid[i][j] > 0:
                top.append((grid[i][j], lat, lon))
    top.sort(reverse=True)
    print("top cells (score, lat S, lon E):")
    for s, la, lo in top[:5]:
        print(f"  {s / best:.3f} of max   {la:.2f}, {lo:.2f}")

    # ---- drift calibration vs holdout validation ----
    if top:
        _, t_lat, t_lon = top[0]
        print("\ndrift check from top cell "
              f"({t_lat:.1f}, {t_lon:.1f}) - closest approach over the window:")
        for item in debris:
            d, day = min_distance_to(
                t_lat, t_lon, item["drift_days"],
                item["find_lat"], item["find_lon"],
            )
            tag = "calibration" if "flaperon" in item["item"] else "held-out"
            print(f"  {tag:12s} {item['item'][:38]:40s} "
                  f"{d:7.0f} nm (day ~{day})")

    # ---- plot ----
    seg = arc_segment(arc)
    fig, ax = plt.subplots(figsize=(11, 8))
    LON, LAT = np.meshgrid(lons, lats)   # both shape (len(lats), len(lons))
    Z = np.array([[grid[i][j] / best for j in range(len(lons))] for i in range(len(lats))])
    im = ax.pcolormesh(LON, LAT, Z, cmap="inferno", shading="auto")
    fig.colorbar(im, ax=ax, label="relative probability (max = 1)")

    if seg:
        ax.plot([p[1] for p in seg], [p[0] for p in seg], "w--", lw=1.2, label="seventh arc")
    ax.scatter(*ARC_ANCHOR, marker="x", c="cyan", s=60, zorder=5)
    ax.annotate("documented arc anchor\n(34.13S 93.95E)", ARC_ANCHOR,
                xytext=(8, -14), textcoords="offset points", color="cyan", fontsize=8)

    for item in debris:
        if not (-42.0 <= item["find_lat"] <= -12.0 and 48.0 <= item["find_lon"] <= 112.0):
            continue  # far-field finds (e.g. Caribbean) stay in the scoring, off this map
        ax.scatter([item["find_lon"]], [item["find_lat"]], marker="^", c="lime", s=50, zorder=5)
        ax.annotate(item["item"].split(" ")[0] + " find", (item["find_lon"], item["find_lat"]),
                    xytext=(6, 4), textcoords="offset points", color="lime", fontsize=7)

    # official reference estimates for comparison
    ax.scatter([93.0], [-34.0], marker="s", c="white", s=45, zorder=5)
    ax.annotate("ATSB unsearched-area centre (34S 93E)", (93.0, -34.0),
                xytext=(-128, -6), textcoords="offset points", color="white", fontsize=7)
    ax.scatter([92.8], [-35.6], marker="s", c="yellow", s=45, zorder=5)
    ax.annotate("CSIRO most-likely impact (35.6S 92.8E)", (92.8, -35.6),
                xytext=(8, 10), textcoords="offset points", color="yellow", fontsize=7)

    ax.set_xlim(48.0, 112.0)
    ax.set_ylim(-42.0, -12.0)
    ax.set_xlabel("longitude (deg E)")
    ax.set_ylabel("latitude (deg S shown as negative)")
    ax.set_title("MH370: candidate impact probability along the seventh arc\n"
                 "(path likelihood x BFO agreement x debris-drift agreement; v1.5 model)")
    ax.legend(loc="upper right", fontsize=8)
    os.makedirs(os.path.join(REPO_ROOT, "output"), exist_ok=True)
    out = os.path.join(REPO_ROOT, "output", "mh370_probability_heatmap.png")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
