from ..alanine_dipeptide import fetch_alanine_dipeptide


def test_ala():
    data = fetch_alanine_dipeptide()
    assert len(data['trajectories']) == 10
