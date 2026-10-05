import pandas as pd
import pytest

import splice_assay as sa
from splice_assay import InputError
from splice_assay.cli import main
from splice_assay.probe import probe

FAST = sa.Settings(formats=("png",), dpi=60)
NAMES = {"tumour": "Responder", "normal": "Non-responder"}


def _renamed(tables, names=NAMES, role=True):
    s = tables["samples"].copy()
    s["group"] = s.group.map(names)
    if role:
        s["role"] = s.group.map({names["tumour"]: "case", names["normal"]: "control"})
    return dict(tables, samples=s)


def _texts(fig) -> str:
    out = [t.get_text() for t in fig.texts]
    for ax in fig.axes:
        out += [t.get_text() for t in ax.texts] + [ax.get_xlabel(), ax.get_ylabel()]
        out += [t.get_text() for t in ax.get_xticklabels()]
    return "\n".join(out)


def test_role_column_names_the_comparison(tables):
    ds = sa.Dataset.from_tables(**_renamed(tables))
    assert ds.labels == {"case": "Responder", "reference": "Non-responder"}
    a = sa.analyze(sa.Dataset.from_tables(**tables), events=["SYN1:SE:1"])
    b = sa.analyze(ds, events=["SYN1:SE:1"])
    pd.testing.assert_frame_equal(a.groups, b.groups)                 # same statistics, other names
    pd.testing.assert_frame_equal(a.survival, b.survival)
    with pytest.raises(InputError, match="role column"):                # neither flags nor a role column
        sa.Dataset.from_tables(**_renamed(tables, role=False))
    flags = sa.Dataset.from_tables(**_renamed(tables, role=False), case="Responder", reference="Non-responder")
    assert flags.labels == ds.labels
    bad = _renamed(tables)
    bad["samples"].loc[bad["samples"].index[0], "role"] = "control"   # one Responder sample marked as control
    with pytest.raises(InputError, match="has both roles"):
        sa.Dataset.from_tables(**bad)


def test_figures_and_report_use_the_names(tables, tmp_path):
    ds = sa.Dataset.from_tables(**_renamed(tables))
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH4"], "OS", settings=FAST,
                       detail=sa.CoxModel().with_clinical(("age",)))
    text = _texts(p.figure)
    for want in ("Higher in responder", "(responder − non-responder)", "Non-responder", "Responder",
                 "no responder vs non-responder"):
        assert want in text, want
    assert "umour" not in text and "ormal" not in text
    g = _texts(sa.expression_panel(ds, "SYN1", ["COH1", "COH4"], "OS", settings=FAST,
                                   model=sa.CoxModel().with_clinical(("age",))).figure)
    assert "responder vs non-responder, a KM split" in g and "umour" not in g and "ormal" not in g
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path / "pr", top=1, max_pages=1, log=lambda *_: None)
    report = res.paths["report"].read_text()
    assert "responder–non-responder" in report and "umour" not in report and "ormal" not in report


def test_names_keep_their_capitals(tables):
    ds = sa.Dataset.from_tables(**_renamed(tables, {"tumour": "IDH-mutant", "normal": "IDH-wildtype"}))
    text = _texts(sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST).figure)
    assert "Higher in IDH-mutant" in text and "(IDH-mutant − IDH-wildtype)" in text


def test_cli_names_from_flags_or_role(tables, tmp_path, capsys):
    d = tmp_path / "data"
    d.mkdir()
    t = _renamed(tables, role=False)
    for name in ("samples", "psi", "events", "survival", "expression", "clinical"):
        t[name].to_csv(d / f"{name}.csv", index=False)
    assert main(["validate", str(d)]) == 2                               # no role column, no flags
    assert "add a role column" in capsys.readouterr().err
    assert main(["validate", str(d), "--case", "Responder", "--reference", "Non-responder"]) == 0
    assert "groups: case = Responder, reference = Non-responder" in capsys.readouterr().out
    _renamed(tables)["samples"].to_csv(d / "samples.csv", index=False)  # now with the role column
    assert main(["validate", str(d)]) == 0
    assert "groups: case = Responder, reference = Non-responder" in capsys.readouterr().out
