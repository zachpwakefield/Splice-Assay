import pandas as pd
import pytest

from splice_assay import InputError, geometry
from splice_assay.events import parse_intervals


def row(event_type, constant, variable, **kw):
    return pd.Series({**dict(event_type=event_type, constant=constant, variable=variable, psi_junctions="",
                             other_junctions=""), **kw})


def test_parse_intervals():
    assert parse_intervals("1-5; 10-20") == [(1, 5), (10, 20)]
    assert parse_intervals("") == []
    with pytest.raises(InputError):
        parse_intervals("5-1")
    with pytest.raises(InputError):
        parse_intervals("1:5")


def test_se():
    g = geometry("e", row("SE", "100-200;500-600", "300-350"))
    assert g.psi_arcs == ((200, 300), (350, 500)) and g.other_arcs == ((200, 500),)
    assert g.variable_length == 50 and g.span == (100, 600)


def test_ri_must_be_the_intron():
    g = geometry("e", row("RI", "100-200;500-600", "200-500"))
    assert g.psi_arcs == () and g.other_arcs == ((200, 500),)
    with pytest.raises(InputError, match="retained intron"):
        geometry("e", row("RI", "100-200;500-600", "250-500"))


@pytest.mark.parametrize("constant, variable, psi, other", [
    ("100-200;500-600", "470-500", ((200, 470),), ((200, 500),)),     # flank upstream, long form extends left
    ("100-200;500-600", "200-230", ((230, 500),), ((200, 500),)),     # flank downstream, long form extends right
])
def test_alternative_sites(constant, variable, psi, other):
    for t in ("A3SS", "A5SS"):
        g = geometry("e", row(t, constant, variable))
        assert g.psi_arcs == psi and g.other_arcs == other


def test_alternative_site_needs_adjacent_exon():
    with pytest.raises(InputError, match="extend exactly one"):
        geometry("e", row("A3SS", "100-200;500-600", "300-320"))


def test_mxe_psi_exon_first():
    g = geometry("e", row("MXE", "100-200;900-1000", "600-700;300-400"))
    assert g.variable == ((600, 700),) and g.other == ((300, 400),)
    assert g.psi_arcs == ((200, 600), (700, 900)) and g.other_arcs == ((200, 300), (400, 900))


def test_other_types_use_given_junctions():
    g = geometry("e", row("AFE", "500-600", "100-200", psi_junctions="200-500"))
    assert g.psi_arcs == ((200, 500),) and g.other_arcs == ()
