from __future__ import absolute_import

from msmb_data import (fetch_alanine_dipeptide, fetch_fs_peptide,
                       fetch_met_enkephalin, load_doublewell, load_quadwell,
                       load_muller)


def test_alanine_dipeptide():
    data = fetch_alanine_dipeptide()
    assert len(data['trajectories']) == 10


def test_fs_peptide():
    data = fetch_fs_peptide()
    assert len(data['trajectories']) == 28


def test_met_enkephalin():
    data = fetch_met_enkephalin()
    assert len(data['trajectories']) == 10


def test_doublewell():
    data = load_doublewell()
    assert data is not None


def test_quadwell():
    data = load_quadwell()
    assert data is not None


def test_muller():
    data = load_muller()
    assert data is not None
