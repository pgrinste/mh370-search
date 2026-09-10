# mh370-search

Where is Flight MH370? A probabilistic search of the southern Indian Ocean built from public data: radar fixes, satellite ATOC timestamps, the 7th sonar ping arc, and confirmed debris locations.

## What it does

- Reconstructs candidate flight paths constrained by the last military radar contact and satellite handshake times
- Simulates how debris would drift on ocean currents and wind
- Scores every grid cell against where debris actually washed up
- Produces a probability heatmap of the most likely resting location

## Layout

```
data/curated/   cleaned, documented datasets (committed)
src/mh370/      package: geodesy, flight paths, drift, scoring
tests/          unit tests
notebooks/      local exploration (Kaggle version lives on Kaggle)
```

## Setup

```bash
pip install -r requirements.txt
python -m pytest
```
