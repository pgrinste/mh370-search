# Data dictionary

All files in this directory are curated from public sources. Coordinates are WGS84 decimal degrees (negative = south/west). Times are UTC. Access date for all web sources: 2026-09-10.

## mh370_timeline.csv
Key events of the flight and satellite contact, with confidence levels.
- `exact` — time/position as published by the investigating body or Inmarsat.
- `approximate` — derived from press reporting of official statements; treat as ±30 min / ±50 NM.

## mh370_debris.csv
Recovered items attributed (to varying degrees) to MH370.
- `status`: confirmed = formally identified by Boeing/Airbus or the investigating authority; probable = examined and considered consistent but not formally confirmed; unconfirmed = private claim.
- Coordinates for beach finds are approximate to ~±0.1 deg (beach location, not GPS fix of the find).

## mh370_search_areas.csv
Official search areas over time. Empty cells mean the boundary is defined along the arc rather than by a rectangular box; see `definition`.

## Primary sources
- Inmarsat: "The Search for MH370", C. Ashton, A.S. Bruce, G. Colledge, M. Dickinson — Journal of Navigation 68 (2015). https://doi.org/10.1017/S037346331400068X
- ATSB: "The Operational Search for MH370" (final report, Feb 2016). https://www.atsb.gov.au
- JACC search-area statements as reported by MarineLink, 20 Jun 2014. https://www.marinelink.com/news/delineated-seventh-search371484
- Wikipedia: "Search for Malaysia Airlines Flight 370" (aggregates CSIRO drift reports and ATSB FPR). https://en.wikipedia.org/wiki/Search_for_Malaysia_Airlines_Flight_370
- R. Godfrey-White, "The Search for MH370" (ongoing public analysis; seventh-arc geometry notes). https://www.mh370search.com

## Known simplifications
- The seventh arc is modelled as a circle on the Earth's surface centred on the Inmarsat-3 F1 sub-satellite point (0 deg, 64.5 deg E) with radius fitted to the documented anchor point (34.13 S, 93.95 E). Strictly, BTO fixes slant range to the satellite in space; the ground-circle approximation is good to a few tens of km for this use and is validated against JACC's stated arc extent (lat 20S-39S) in `tests/test_ping_arc.py`.
