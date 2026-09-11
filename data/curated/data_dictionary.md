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

## mh370_satcom.csv
Inmarsat SATCOM handshake records (BTO/BFO) curated from the Inmarsat SU log released by the Malaysian government (May 2014), cross-checked against Ashton et al. 2014 Table 6 and Holland 2017.
- `bto_us_corrected`: burst timing offset in microseconds, with the -4600 us correction of Ashton et al. 2014 note (1) applied to arcs 1 and 7. BTO values are residuals after removing the nominal satellite model; per the ATSB factsheet their tolerance is ±10 km.
- `bfo_hz`: burst frequency offset in Hz as recorded. The ~150 Hz bias (satellite oscillator drift + AES compensation, Holland 2017) must be removed before comparing with geometric Doppler.
- **Time base:** times are given in UTC. Early entries were published as Malaysian time (UTC+8) and converted; the final handshake is 00:19:29Z per Holland 2017 (= 08:19 MYT). Public reporting often quotes "08:19" without a timezone - that is the source of most confusion in secondary analyses. The `ges_handshake_6` time (published '00:10') carries ±30 min uncertainty.

## Primary sources
- Inmarsat: "The Search for MH370", C. Ashton, A.S. Bruce, G. Colledge, M. Dickinson — Journal of Navigation 68 (2015). https://doi.org/10.1017/S037346331400068X
- ATSB: "The Operational Search for MH370" (final report, Feb 2016). https://www.atsb.gov.au
- JACC search-area statements as reported by MarineLink, 20 Jun 2014. https://www.marinelink.com/news/delineated-seventh-search371484
- Wikipedia: "Search for Malaysia Airlines Flight 370" (aggregates CSIRO drift reports and ATSB FPR). https://en.wikipedia.org/wiki/Search_for_Malaysia_Airlines_Flight_370
- R. Godfrey-White, "The Search for MH370" (ongoing public analysis; seventh-arc geometry notes). https://www.mh370search.com

## Known simplifications
- The seventh arc is modelled as the exact constant-slant-range locus: intersection of a sphere centred on Inmarsat-3 F1 (geostationary slot, 42164 km from Earth centre) with the WGS84 ellipsoid, fitted to pass through the documented anchor point (34.13 S, 93.95 E). The v1 ground-circle model agreed with this to ~10 km over the southern band once both were calibrated on the same anchor (`tests/test_ping_arc.py::test_ground_circle_agreement`); the difference is purely ellipsoidal.
- BTO values are not used numerically (the arc is anchored on the documented point instead); they are kept in `mh370_satcom.csv` for provenance and future work beyond the seventh arc, which will need ATSB's BTO-R residual table.
