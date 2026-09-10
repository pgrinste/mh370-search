from mh370.geodesy import destination_point, haversine_nm, initial_bearing_deg


def test_haversine_known_distance():
    # 1 degree of latitude at the equator is ~60 nm
    d = haversine_nm(0.0, 0.0, 1.0, 0.0)
    assert abs(d - 60.0) < 0.5


def test_destination_roundtrip():
    lat, lon = destination_point(-20.0, 90.0, 245.0, 300.0)
    back = haversine_nm(lat, lon, -20.0, 90.0)
    assert abs(back - 300.0) < 1.0


def test_bearing_north():
    b = initial_bearing_deg(0.0, 50.0, 10.0, 50.0)
    assert abs(b - 0.0) < 0.1 or abs(b - 360.0) < 0.1
