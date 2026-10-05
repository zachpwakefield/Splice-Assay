import numpy as np
import pandas as pd
import pytest

import splice_assay as sa
from splice_assay import InputError
from splice_assay.dataset import derive_pairs, keep_ids, normalise_samples


def test_wide_and_long_psi_are_the_same(tables):
    wide = tables["psi"].pivot(index="event_id", columns="sample_id", values="psi").reset_index()
    a = sa.Dataset.from_tables(**tables)
    b = sa.Dataset.from_tables(**dict(tables, psi=wide))
    pd.testing.assert_frame_equal(a.psi, b.psi)
    assert list(a.psi.columns) == list(a.samples.sample_id)


def test_long_and_wide_expression_are_the_same(tables):
    wide = tables["expression"].pivot(index="gene", columns="sample_id", values="value").reset_index()
    a = sa.Dataset.from_tables(**tables)
    b = sa.Dataset.from_tables(**dict(tables, expression=wide))
    pd.testing.assert_frame_equal(a.expression, b.expression)


def test_from_dir_keeps_ids_as_text(tmp_path, tables):
    for name, df in tables.items():
        df.to_csv(tmp_path / f"{name}.csv", index=False)
    s = tables["samples"].copy()
    s["patient_id"] = s.patient_id.str.replace("COH1-P", "00", regex=False)   # numeric-looking IDs
    s.to_csv(tmp_path / "samples.csv", index=False)
    surv = tables["survival"].copy()
    surv["patient_id"] = surv.patient_id.str.replace("COH1-P", "00", regex=False)
    surv.to_csv(tmp_path / "survival.tsv", sep="\t", index=False)
    (tmp_path / "survival.csv").unlink()
    ds = sa.Dataset.from_dir(tmp_path)
    assert ds.samples.patient_id.str.startswith("000").any()
    assert ds.summary().set_index("cohort").loc["COH1", "OS_patients"] > 150


@pytest.mark.parametrize("change, message", [
    (lambda s: s.drop(columns="group"), "missing column"),
    (lambda s: pd.concat([s, s.iloc[:1]]), "unique"),
    (lambda s: s.assign(group="blood"), "no sample of the case group"),
])
def test_sample_errors(tables, change, message):
    with pytest.raises(InputError, match=message):
        sa.Dataset.from_tables(**dict(tables, samples=change(tables["samples"])))


def test_tumor_spelling_and_tissue_column_are_accepted(tables):
    s = tables["samples"].rename(columns={"group": "tissue"})
    s["tissue"] = s.tissue.str.replace("tumour", "Tumor")
    ds = sa.Dataset.from_tables(**dict(tables, samples=s))
    assert set(ds.samples.role) == {"case", "reference"} and ds.labels == dict(case="Tumor", reference="Normal")


def test_other_groups_are_set_aside(tables):
    s = tables["samples"].copy()
    s.loc[s.index[:3], "group"] = "metastasis"
    with pytest.warns(UserWarning, match="other groups"):
        ds = sa.Dataset.from_tables(**dict(tables, samples=s))
    assert ds.notes["other_groups"] == {"metastasis": 3} and (ds.samples.role == "other").sum() == 3


def test_psi_percentages_are_refused(tables):
    p = tables["psi"].assign(psi=tables["psi"].psi * 100)
    with pytest.raises(InputError, match="percentages"):
        sa.Dataset.from_tables(**dict(tables, psi=p))


def test_invalid_survival_rows_are_dropped_and_counted(tables):
    sv = tables["survival"].copy()
    sv.loc[sv.index[:3], "time"] = 0
    sv.loc[sv.index[3], "event"] = 2
    ds = sa.Dataset.from_tables(**dict(tables, survival=sv))
    assert sum(ds.notes["survival_rows_dropped"].values()) == 4
    same = sa.Dataset.from_tables(**dict(tables, survival=pd.concat([sv, sv.iloc[10:11]])))   # identical: merged
    assert len(same.survival) == len(sa.Dataset.from_tables(**dict(tables, survival=sv)).survival)
    clash = sv.iloc[10:11].assign(time=sv.time.iloc[10] + 5)
    with pytest.raises(InputError, match="conflicting"):
        sa.Dataset.from_tables(**dict(tables, survival=pd.concat([sv, clash])))


def test_second_tumour_needs_a_survival_choice(tables):
    s = tables["samples"]
    extra = s.iloc[[0]].assign(sample_id="COH1-P0000-T2")
    with pytest.raises(InputError, match="survival_cohort"):
        normalise_samples(pd.concat([s, extra]))
    both, _ = normalise_samples(pd.concat([s, extra]).assign(
        survival_cohort=lambda d: d.group.eq("tumour") & d.sample_id.ne("COH1-P0000-T2")))
    with pytest.raises(InputError, match="pairs table"):
        derive_pairs(both)                                  # P0000 has a normal and two tumours


def test_explicit_pairs_are_checked(tables):
    ds = sa.Dataset.from_tables(**tables)
    pairs = ds.pairs.copy()
    ok = sa.Dataset.from_tables(**dict(tables, pairs=pairs))
    assert len(ok.pairs) == len(pairs) and ok.notes["pairs_source"] == "table"
    bad = pairs.copy()
    bad.loc[0, "reference_sample"] = bad.loc[0, "case_sample"]
    with pytest.raises(InputError, match="different group"):
        sa.Dataset.from_tables(**dict(tables, pairs=bad))


def test_event_subset(tables):
    ds = sa.Dataset.from_tables(**tables, event_ids=["SYN2:MXE:1"])
    assert list(ds.psi.index) == ["SYN2:MXE:1"] and list(ds.expression.index) == ["SYN2"]
    with pytest.raises(InputError, match="missing"):
        sa.Dataset.from_tables(**tables, event_ids=["NOPE"])


def test_settings_roundtrip(tmp_path):
    s = sa.Settings(min_pairs=5, formats=("png",))
    p = tmp_path / "s.json"
    import json
    p.write_text(json.dumps(s.to_dict()))
    assert sa.Settings.from_json(p) == s
    with pytest.raises(ValueError, match="unknown"):
        sa.Settings.from_dict({"min_pair": 3})
    with pytest.raises(ValueError):
        sa.Settings(alpha=2)


def test_custom_group_names_and_columns_give_identical_results(tables):
    """Rename every column and both groups: the statistics must not change."""
    base = sa.analyze(sa.Dataset.from_tables(**tables), events=["SYN1:SE:1"])
    s = tables["samples"].assign(group=tables["samples"].group.map({"tumour": "Metastasis", "normal": "Primary"}))
    s = s.rename(columns={"sample_id": "File.ID", "patient_id": "Case.ID", "cohort": "project", "group": "site"})
    sv = tables["survival"].rename(columns={"patient_id": "Case.ID"})
    ps = tables["psi"].rename(columns={"sample_id": "File.ID"})
    ex = tables["expression"].rename(columns={"sample_id": "File.ID"})
    cols = {"sample_id": "File.ID", "patient_id": "Case.ID", "cohort": "project", "group": "site"}
    ds = sa.Dataset.from_tables(samples=s, psi=ps, events=tables["events"], survival=sv, expression=ex,
                                case="metastasis", reference="primary", columns=cols)
    assert ds.labels == dict(case="Metastasis", reference="Primary")
    alt = sa.analyze(ds, events=["SYN1:SE:1"])
    num = ["paired_p", "paired_delta_median", "unpaired_p", "unpaired_delta_median"]
    assert np.allclose(base.groups[num].to_numpy(float), alt.groups[num].to_numpy(float), equal_nan=True)
    num = ["cutoff", "km_p", "cox_beta", "cox_se", "hr_per_iqr"]
    assert np.allclose(base.survival[num].to_numpy(float), alt.survival[num].to_numpy(float), equal_nan=True)


def test_wide_survival_columns_are_detected(tables):
    sv = tables["survival"]
    wide = sv.pivot(index="patient_id", columns="endpoint", values=["time", "event"])
    wide = pd.DataFrame({"patient_id": wide.index, "OS.time": wide[("time", "OS")].to_numpy(),
                         "OS": wide[("event", "OS")].to_numpy(), "DSS_time": wide[("time", "DSS")].to_numpy(),
                         "DSS_event": wide[("event", "DSS")].to_numpy()})
    a = sa.Dataset.from_tables(**tables).survival.sort_values(["endpoint", "patient_id"]).reset_index(drop=True)
    b = sa.Dataset.from_tables(**dict(tables, survival=wide)).survival
    b = b.sort_values(["endpoint", "patient_id"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b)


def test_missing_codes(tables, tmp_path):
    """A clinical code such as 'missing' is a category unless it is declared missing."""
    cl = tables["clinical"].astype({"stage": object})
    cl.loc[cl.index[:5], "stage"] = "missing"
    kept = sa.Dataset.from_tables(**dict(tables, clinical=cl))
    assert (kept.clinical.stage == "missing").sum() == 5
    dropped = sa.Dataset.from_tables(**dict(tables, clinical=cl), na_values=["missing"])
    assert dropped.clinical.stage.isna().sum() == 5
    cl.to_csv(tmp_path / "clinical.csv", index=False)
    for k in ("samples", "psi", "events", "survival"):
        tables[k].to_csv(tmp_path / f"{k}.csv", index=False)
    assert sa.Dataset.from_dir(tmp_path, na_values=["missing"]).clinical.stage.isna().sum() == 5


def _two_tables(t):
    """The example as two tables: samples with survival and clinical columns; psi with the event columns."""
    w = t["survival"].pivot(index="patient_id", columns="endpoint", values=["time", "event"])
    w.columns = [f"{ep}.time" if what == "time" else ep for what, ep in w.columns]
    samples = t["samples"].merge(w.reset_index(), on="patient_id", how="left").merge(t["clinical"], on="patient_id",
                                                                                      how="left")
    psi = t["events"].merge(t["psi"].pivot(index="event_id", columns="sample_id", values="psi").reset_index(),
                            on="event_id")
    return samples, psi


def test_two_tables_are_enough(tables):
    samples, psi = _two_tables(tables)
    one = sa.Dataset.from_tables(**tables)
    two = sa.Dataset.from_tables(samples=samples, psi=psi, expression=tables["expression"])
    assert set(two.notes["sources"]) == {"events", "survival", "clinical"}
    pd.testing.assert_frame_equal(one.psi.sort_index(), two.psi.sort_index())      # same values; file order kept
    pd.testing.assert_frame_equal(one.events.sort_index(), two.events.sort_index())
    key = ["patient_id", "endpoint"]
    pd.testing.assert_frame_equal(one.survival.sort_values(key).reset_index(drop=True),
                                  two.survival.sort_values(key).reset_index(drop=True))
    pd.testing.assert_frame_equal(one.clinical.sort_index(), two.clinical.sort_index(), check_like=True)
    pd.testing.assert_frame_equal(one.pairs, two.pairs)
    m = sa.CoxModel().with_clinical(("age", "stage"))
    a = sa.analyze(one, events=["SYN1:SE:1"], model=m).survival
    b = sa.analyze(two, events=["SYN1:SE:1"], model=m).survival
    pd.testing.assert_frame_equal(a, b)


def test_long_psi_with_event_columns(tables):
    psi = tables["psi"].merge(tables["events"], on="event_id")
    ds = sa.Dataset.from_tables(samples=tables["samples"], psi=psi)
    pd.testing.assert_frame_equal(ds.events.sort_index(), sa.Dataset.from_tables(**tables).events.sort_index())
    bad = psi.copy()
    bad.loc[bad.index[0], "gene"] = "OTHER"
    with pytest.raises(InputError, match="event columns differ"):
        sa.Dataset.from_tables(samples=tables["samples"], psi=bad)
    with pytest.raises(InputError, match="needs the event columns"):
        sa.Dataset.from_tables(samples=tables["samples"], psi=tables["psi"])


def test_pair_id_in_the_samples_table(tables):
    s = tables["samples"].copy()
    paired = set(s.patient_id[s.group.eq("normal")]) & set(s.patient_id[s.group.eq("tumour")])
    s["pair_id"] = np.where(s.patient_id.isin(paired), s.patient_id, "")
    drop = sorted(paired)[0]
    s.loc[s.patient_id.eq(drop), "pair_id"] = ""                       # one pair left out on purpose
    ds = sa.Dataset.from_tables(**dict(tables, samples=s))
    assert ds.notes["pairs_source"] == "samples.pair_id"
    assert len(ds.pairs) == len(sa.Dataset.from_tables(**tables).pairs) - 1 and drop not in set(ds.pairs.patient_id)
    s.loc[s.patient_id.eq(sorted(paired)[1]), "pair_id"] = sorted(paired)[2]   # four samples under one pair_id
    with pytest.raises(InputError, match="pair_id must name one case and one reference"):
        sa.Dataset.from_tables(**dict(tables, samples=s))


def test_from_dir_with_two_tables_and_cli(tables, tmp_path, capsys):
    samples, psi = _two_tables(tables)
    d = tmp_path / "two"
    d.mkdir()
    samples.to_csv(d / "samples.csv", index=False)
    psi.to_csv(d / "psi.csv", index=False)
    tables["expression"].to_csv(d / "expression.csv", index=False)
    ds = sa.Dataset.from_dir(d)
    assert ds.endpoints == ["DSS", "OS"] and {"age", "sex", "stage"} <= set(ds.clinical.columns)
    from splice_assay.cli import main
    assert main(["validate", str(d)]) == 0
    out = capsys.readouterr().out
    assert "events: from the psi table" in out and "survival: from the samples table (DSS, OS)" in out
    assert main(["probe", str(d), "--gene", "SYN2", "--out", str(tmp_path / "pr"), "--max-pages", "0"]) == 0
    assert "MXE:1 (SYN2)" in capsys.readouterr().out


def test_keep_a_subset_of_patients(tables, tmp_path):
    s = tables["samples"]
    paired = sorted(set(s.patient_id[s.group.eq("normal")]) & set(s.patient_id[s.cohort.eq("COH1")]))[:12]
    lst = pd.DataFrame({"Case.ID": [f"{p}-T-01A-11R" for p in paired], "Subtype": "X"})   # barcode-style IDs
    f = tmp_path / "keep.tsv"
    lst.to_csv(f, sep="\t", index=False)
    ds = sa.Dataset.from_tables(**tables, keep=f)
    assert set(ds.samples.patient_id) == set(paired)
    assert ds.samples.role.value_counts().to_dict() == {"case": 12, "reference": 12}  # their normals are kept
    assert len(ds.pairs) == 12 and ds.psi.shape[1] == 24
    assert ds.notes["subset"]["patients"] == 12 and ds.notes["subset"]["listed"] == 12
    by_list = sa.Dataset.from_tables(**tables, keep=paired)                      # plain patient IDs work too
    pd.testing.assert_frame_equal(by_list.samples, ds.samples)
    with pytest.raises(InputError, match="none of the 2 listed IDs"):
        sa.Dataset.from_tables(**tables, keep=["nobody", "TCGA-XX-0000"])
    pairs = derive_pairs(normalise_samples(tables["samples"])[0])               # an explicit pairs table is cut too
    sub = sa.Dataset.from_tables(**dict(tables, pairs=pairs), keep=paired)
    assert len(sub.pairs) == 12


def test_keep_on_the_command_line(tables, tmp_path, capsys):
    from splice_assay.cli import main
    from splice_assay.probe import probe
    d = tmp_path / "data"
    d.mkdir()
    for name in ("samples", "psi", "events", "survival", "expression", "clinical"):
        tables[name].to_csv(d / f"{name}.csv", index=False)
    s = tables["samples"]
    keep = sorted(set(s.patient_id[s.cohort.eq("COH1")]))
    pd.DataFrame({"patient": keep}).to_csv(tmp_path / "keep.csv", index=False)
    assert main(["validate", str(d), "--keep", str(tmp_path / "keep.csv")]) == 0
    out = capsys.readouterr().out
    assert f"subset: {len(keep)} of " in out and "cohorts 1 " in out
    ds = sa.Dataset.from_dir(d, keep=tmp_path / "keep.csv")
    assert ds.cohorts == ["COH1"]
    res = probe(ds, genes=["SYN1"], settings=sa.Settings(formats=("png",), dpi=60), out_dir=tmp_path / "pr", top=1,
                max_pages=1, log=lambda *_: None)
    assert "- **Subset:**" in res.paths["report"].read_text()
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=sa.Settings(formats=("png",), dpi=60))
    assert any("subset: keep" in t.get_text() for t in p.figure.texts)


def test_where_selects_patients_by_column_values(tables):
    """where: a clinical or samples column, case-insensitive, every condition must hold, normals stay, keep combines."""
    cl, s = tables["clinical"], tables["samples"]
    fem = set(cl.patient_id[cl.sex.eq("female")]) & set(s.patient_id)
    ds = sa.Dataset.from_tables(**tables, where=["sex=FEMALE"])                    # a clinical column
    assert set(ds.samples.patient_id) == fem and ds.samples.role.eq("reference").any()
    assert ds.notes["subset"]["by"] == "sex = FEMALE" and ds.notes["subset"]["name"] == "sex=FEMALE"
    both = sa.Dataset.from_tables(**tables, where={"sex": "female", "stage": ["II", "III"]})
    want = set(cl.patient_id[cl.sex.eq("female") & cl.stage.astype(str).isin(["II", "III"])]) & set(s.patient_id)
    assert set(both.samples.patient_id) == want
    assert sa.Dataset.from_tables(**tables, where="cohort=COH1").cohorts == ["COH1"]      # a samples column
    kept = sorted(fem)[:20]
    k = sa.Dataset.from_tables(**tables, keep=kept, where=["cohort=COH1"])        # keep and where both apply
    assert set(k.samples.patient_id) == {p for p in kept if p.startswith("COH1")}
    assert k.notes["subset"]["listed"] == 20 and " and cohort = COH1" in k.notes["subset"]["by"]
    with pytest.raises(InputError, match="no column subtype"):
        sa.Dataset.from_tables(**tables, where=["subtype=Her2"])
    with pytest.raises(InputError, match="expected COLUMN=VALUE"):
        sa.Dataset.from_tables(**tables, where=["sex"])


def test_where_on_the_command_line(tables, tmp_path, capsys):
    from splice_assay.cli import main
    d = tmp_path / "data"
    d.mkdir()
    for name in ("samples", "psi", "events", "survival", "clinical"):
        tables[name].to_csv(d / f"{name}.csv", index=False)
    assert main(["validate", str(d), "--where", "sex=female", "--where", "cohort=COH1,COH2"]) == 0
    assert "selected by sex = female and cohort = COH1 or COH2" in capsys.readouterr().out


def _tiny(pids, **cols):
    """A samples table (one tumour and one normal sample per patient, in cohort C, plus `cols`) and a psi table."""
    rows = [dict(sample_id=f"{p}-{k}", patient_id=p, cohort="C", group=g, **{c: v[i] for c, v in cols.items()})
            for i, p in enumerate(pids) for k, g in (("T", "tumour"), ("N", "normal"))]
    s = pd.DataFrame(rows)
    return s, pd.DataFrame({"event_id": ["E1"], "gene": ["G"], **{x: [0.5] for x in s.sample_id}})


def test_keep_lists_are_read_as_text(tmp_path):
    """A listed ID keeps its leading zeros under any header, and a list without a header keeps its first ID."""
    s, psi = _tiny(["00123", "123"])
    (tmp_path / "keep.csv").write_text("id\n00123\n")
    assert keep_ids(tmp_path / "keep.csv") == ["00123"]
    kept = sa.Dataset.from_tables(s, psi, keep=tmp_path / "keep.csv")
    assert sorted(kept.samples.patient_id.unique()) == ["00123"]
    (tmp_path / "one.txt").write_text("00123\n")                          # no header: the first line is an ID
    with pytest.warns(UserWarning, match="no header"):
        kept = sa.Dataset.from_tables(s, psi, keep=tmp_path / "one.txt")
    assert sorted(kept.samples.patient_id.unique()) == ["00123"]
    (tmp_path / "two.txt").write_text("123\n00123\n")
    with pytest.warns(UserWarning, match="no header"):
        assert sa.Dataset.from_tables(s, psi, keep=tmp_path / "two.txt").samples.patient_id.nunique() == 2
    for name, text in (("header.csv", "patient_id\n"), ("nothing.txt", ""), ("blank.txt", "\n\n")):
        (tmp_path / name).write_text(text)
        with pytest.raises(InputError, match="the keep list is empty"):
            sa.Dataset.from_tables(s, psi, keep=tmp_path / name)
    with pytest.raises(InputError, match="keep: cannot read"):
        sa.Dataset.from_tables(s, psi, keep=tmp_path / "missing.csv")


def test_repeated_where_conditions_must_all_hold():
    s, psi = _tiny(["P1", "P2", "P3"], stage=["I", "II", "III"])
    ds = sa.Dataset.from_tables(s, psi, where=["stage=I,II", "stage=II,III"])
    assert sorted(ds.samples.patient_id.unique()) == ["P2"]
    assert ds.notes["subset"]["by"] == "stage = I or II and stage = II or III"


def test_sample_keyed_tables_with_a_subset():
    """Survival and clinical tables keyed by sample_id may list samples that a subset leaves out, and where reads a
    sample-keyed clinical table; a sample unknown to the samples table is still an error."""
    s, psi = _tiny(["P1", "P2", "P3"])
    sv = pd.DataFrame({"sample_id": ["P1-T", "P2-T", "P3-T"], "OS.time": [100, 200, 300], "OS": [1, 0, 1]})
    cl = pd.DataFrame({"sample_id": ["P1-T", "P2-T", "P3-T"], "age": [50, 60, 70]})
    ds = sa.Dataset.from_tables(s, psi, survival=sv, clinical=cl, keep=["P1", "P2"])
    assert sorted(ds.samples.patient_id.unique()) == ["P1", "P2"] and ds.clinical.loc["P2", "age"] == 60
    assert {"P1", "P2"} <= set(ds.survival.patient_id)
    by_age = sa.Dataset.from_tables(s, psi, survival=sv, clinical=cl, where=["age=50"])
    assert sorted(by_age.samples.patient_id.unique()) == ["P1"]
    with pytest.raises(InputError, match="sample_id not in samples"):
        sa.Dataset.from_tables(s, psi, survival=sv.assign(sample_id=["P1-T", "P2-T", "P9-T"]), keep=["P1"])


def test_where_matches_numbers_by_value_and_never_a_missing_value():
    s, psi = _tiny(["P1", "P2", "P3", "P4"], grade=[1.0, 2.0, np.nan, 1.0])
    assert sorted(sa.Dataset.from_tables(s, psi, where=["grade=1"]).samples.patient_id.unique()) == ["P1", "P4"]
    assert sorted(sa.Dataset.from_tables(s, psi, where=["grade=1.0,2"]).samples.patient_id.unique()) == ["P1", "P2",
                                                                                                         "P4"]
    with pytest.raises(InputError, match="no patient is selected"):
        sa.Dataset.from_tables(s, psi, where=["grade=nan"])


def test_a_listed_id_that_is_a_patient_selects_only_that_patient():
    """1-2 is a patient: it does not also select patient 1; an unknown 1-3 selects 1, its longest known part."""
    s, psi = _tiny(["1", "1-2", "TCGA-XX-0001"])
    pick = (lambda ids: sorted(sa.Dataset.from_tables(s, psi, keep=ids).samples.patient_id.unique()))
    assert pick(["1-2"]) == ["1-2"]
    assert pick(["1-3"]) == ["1"]
    assert pick(["TCGA-XX-0001-01A-11R-0000-07"]) == ["TCGA-XX-0001"]
    assert pick(["TCGA-XX-0001-T"]) == ["TCGA-XX-0001"]                   # a sample ID


def test_where_matches_long_numeric_codes_exactly():
    """Numbers past 2**53 are not exact as floats: such codes match as text only."""
    s, psi = _tiny(["P1", "P2"], code=["12345678901234567", "12345678901234568"])
    ds = sa.Dataset.from_tables(s, psi, where=["code=12345678901234567"])
    assert sorted(ds.samples.patient_id.unique()) == ["P1"]
