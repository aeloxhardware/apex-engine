from core.surface_message import surface_error_message as message

def base(**k):
    d = dict(ambiguous_points=2, surface_holes=0, surface_flipped_edges=0, surface_nonmanifold_edges=0, surface_duplicate_faces=0)
    d.update(k); return d

def test_each_defect_agrees_in_number():
    cases = [
        (dict(surface_holes=1), "présente un trou,"), (dict(surface_holes=3), "présente des trous,"),
        (dict(surface_flipped_edges=3), "présente une face retournée,"), (dict(surface_flipped_edges=40), "présente des faces retournées,"),
        (dict(surface_nonmanifold_edges=1), "présente une arête partagée par plus de deux faces,"),
        (dict(surface_nonmanifold_edges=9), "présente des arêtes partagées par plus de deux faces,"),
        (dict(surface_duplicate_faces=1), "présente une face en double,"), (dict(surface_duplicate_faces=5), "présente des faces en double,"),
    ]
    for kw, expected in cases:
        assert expected in message(base(**kw), 66), kw

def test_self_intersections_and_list():
    assert "une auto-intersection" in message(base(surface_holes=1), 90, self_intersections=1)
    assert "des auto-intersections" in message(base(surface_holes=1), 90, self_intersections=4)
    m = message(base(surface_holes=2, surface_flipped_edges=30, surface_duplicate_faces=1), 80)
    assert "des trous, des faces retournées et une face en double" in m

def test_incoherent_points_agree():
    assert "pour un total de 1 point incohérent," in message(base(ambiguous_points=1, surface_holes=1), 99)
    assert "pour un total de 0 point incohérent," in message(base(ambiguous_points=0, surface_holes=1), 99)
    assert "2 points incohérents," in message(base(surface_holes=1), 99)
    assert "24 689 points incohérents" in message(base(ambiguous_points=24689, surface_holes=1), 66)
