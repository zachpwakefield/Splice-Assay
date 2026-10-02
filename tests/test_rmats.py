import numpy as np
import pandas as pd
import pytest

import splice_assay as sa
from splice_assay import InputError
from splice_assay.rmats import import_rmats, read_events, sample_names

HEAD = "ID\tGeneID\tgeneSymbol\tchr\tstrand\t"
TAIL = "\tID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\tPValue\tFDR\t" \
       "IncLevel1\tIncLevel2\tIncLevelDifference"
COLS = {
    "SE": "exonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE",
    "RI": "riExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE",
    "A3SS": "longExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE",
    "A5SS": "longExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE",
    "MXE": "1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\t"
           "downstreamEE",
}
ROWS = {
    "SE": ['1\t"G1"\t"GENE1"\tchr1\t+\t300\t350\t100\t200\t500\t600'],
    "RI": ['1\t"G1"\t"GENE1"\tchr1\t+\t100\t600\t100\t200\t500\t600'],
    "A3SS": ['1\t"G1"\t"GENE1"\tchr1\t+\t470\t600\t500\t600\t100\t200'],        # long form extends left of short
    "A5SS": ['1\t"G2"\t"GENE2"\tchr2\t-\t470\t600\t500\t600\t100\t200'],
    "MXE": ['1\t"G2"\t"GENE2"\tchr2\t-\t300\t400\t600\t700\t100\t200\t900\t1000',
            '2\t"G1"\t"GENE1"\tchr1\t+\t300\t400\t600\t700\t100\t200\t900\t1000'],
}


@pytest.fixture()
def rmats_dir(tmp_path):
    for t, rows in ROWS.items():
        (tmp_path / f"fromGTF.{t}.txt").write_text(HEAD + COLS[t] + "\n" + "\n".join(rows) + "\n")
        mats = [r + f"\t{r.split(chr(9))[0]}\t1,2\t3,4\t5\t6\t2\t1\t0.5\t0.5\t0.25,NA\t0.75\t-0.5" for r in rows]
        (tmp_path / f"{t}.MATS.JCEC.txt").write_text(HEAD + COLS[t] + TAIL + "\n" + "\n".join(mats) + "\n")
    (tmp_path / "b1.txt").write_text("/data/s1.bam,/data/s2.bam\n")
    (tmp_path / "b2.txt").write_text("/data/s3.bam\n")
    return tmp_path


def test_events_geometry(rmats_dir):
    ev = read_events(rmats_dir).set_index("event_id")
    assert ev.loc["SE:1", ["constant", "variable", "gene"]].tolist() == ["100-200;500-600", "300-350", "GENE1"]
    assert ev.loc["RI:1", "variable"] == "200-500"
    assert ev.loc["A3SS:1", ["constant", "variable"]].tolist() == ["100-200;500-600", "470-500"]
    assert ev.loc["MXE:1", "variable"] == "600-700;300-400"      # minus strand: transcript-upstream exon first
    assert ev.loc["MXE:2", "variable"] == "300-400;600-700"      # plus strand: first-listed exon first
    for e, r in ev.iterrows():
        sa.geometry(e, r)                                        # every imported event can be drawn


def test_first_listed_mxe_option(rmats_dir):
    ev = read_events(rmats_dir, ["MXE"], mxe_psi_exon="first_listed").set_index("event_id")
    assert ev.loc["MXE:1", "variable"] == "300-400;600-700"


def test_psi_from_inclevels(rmats_dir):
    assert sample_names(rmats_dir / "b1.txt") == ["s1", "s2"]
    events, psi = import_rmats(rmats_dir, rmats_dir / "b1.txt", rmats_dir / "b2.txt")
    p = psi.set_index(["event_id", "sample_id"]).psi
    assert p[("SE:1", "s1")] == 0.25 and np.isnan(p[("SE:1", "s2")]) and p[("SE:1", "s3")] == 0.75
    assert len(events) == 6 and psi.sample_id.nunique() == 3


def test_name_count_must_match(rmats_dir):
    with pytest.raises(InputError, match="sample names"):
        import_rmats(rmats_dir, rmats_dir / "b1.txt", names1=["a", "b", "c"])


def test_import_as_one_table(rmats_dir, tmp_path):
    """By default one psi.csv holds the event columns and the PSI values; --separate writes events.csv and a long
    psi.csv."""
    from splice_assay.cli import main
    out = tmp_path / "one"
    assert main(["import-rmats", str(rmats_dir), "--b1", str(rmats_dir / "b1.txt"), "--out", str(out)]) == 0
    assert not (out / "events.csv").exists()
    t = pd.read_csv(out / "psi.csv")
    assert {"event_id", "gene", "constant", "variable", "s1", "s2"} <= set(t.columns)
    same = tmp_path / "same"                                 # --one-table, the old flag, is still accepted
    assert main(["import-rmats", str(rmats_dir), "--b1", str(rmats_dir / "b1.txt"), "--one-table", "--out",
                 str(same)]) == 0
    assert (same / "psi.csv").read_text() == (out / "psi.csv").read_text()
    two = tmp_path / "two"
    assert main(["import-rmats", str(rmats_dir), "--b1", str(rmats_dir / "b1.txt"), "--separate", "--out",
                 str(two)]) == 0
    long = pd.read_csv(two / "psi.csv").pivot(index="event_id", columns="sample_id", values="psi")
    pd.testing.assert_frame_equal(t.set_index("event_id")[["s1", "s2"]].sort_index(),
                                  long[["s1", "s2"]].sort_index(), check_names=False)
