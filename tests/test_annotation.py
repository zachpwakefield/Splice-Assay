import gzip

import pandas as pd
import pytest

from splice_assay.annotation import collapsed_model, gene_model, read_gtf, subset_gtf


def test_collapsed_model_of_the_synthetic_genes(gtf_path):
    gtf = read_gtf(gtf_path, [("chr7", 99000, 107000), ("chr12", 199000, 207000)])
    ex = gtf[gtf.feature.eq("exon") & gtf.gene_name.eq("SYN1")]
    blocks, ntx, thr = collapsed_model(ex)
    # the retained-intron transcript is excluded; the 101970 acceptor is used by one transcript only
    assert ntx == 3 and thr == 2
    assert [b[:2] for b in blocks] == [(100000, 100200), (102000, 102100), (103000, 103150), (104000, 104300),
                                       (105500, 106000)]
    m = gene_model(gtf, "SYN2", "chr12", (200000, 206000), "SYNG0002")
    assert m.strand == "-" and len(m.blocks) == 4
    assert m.blocks[0][0] == 200025                  # free 3' end: median of the transcripts' ends
    assert m.nested.gene_name.tolist() == ["SNORA900"]


def test_chromosome_aliases_and_gene_tags(gtf_path):
    gtf = read_gtf(gtf_path, [("7", 103400, 103600)], genes=["SYN1"])
    assert set(gtf.chrom) == {"7"}                   # returned as spelled in the windows
    assert (gtf.feature.eq("exon") & gtf.gene_name.eq("SYN1")).sum() == 16   # every SYN1 exon, not only the window
    assert gene_model(gtf, "SYN1", "7", (103000, 104300)).nested.gene_name.tolist() == ["SNORD900"]


def test_nested_genes_of_one_name_stay_apart(gtf_path):
    gtf = read_gtf(gtf_path, [("chr7", 99000, 107000)])
    twins = pd.concat([gtf[gtf.feature.eq("gene") & gtf.gene_name.eq("SNORD900")]] * 2, ignore_index=True)
    twins = twins.assign(gene_id=["G1", "G2"], gene_name="SNORD22", start=[100300, 105600], end=[100380, 105680])
    m = gene_model(pd.concat([gtf, twins], ignore_index=True), "SYN1", "chr7", (100000, 106000))
    assert m.nested.gene_name.tolist() == ["SNORD22", "SNORD900", "SNORD22"]   # two genes, not one bar across both


def test_missing_gene_gives_none(gtf_path):
    import warnings
    gtf = read_gtf(gtf_path, [("chr7", 99000, 107000)])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        assert gene_model(gtf, "NOPE", "chr7", (100000, 101000)) is None
    assert w


def test_subset_gtf(gtf_path, tmp_path):
    out = tmp_path / "small.gtf.gz"
    n = subset_gtf(gtf_path, out, [("chr12", 199000, 207000)])
    with gzip.open(out, "rt") as fh:
        lines = fh.read().splitlines()
    assert n == len(lines) > 0 and all(line.startswith("chr12") for line in lines)


def test_the_warning_names_why_a_gene_has_no_track():
    one = pd.DataFrame(dict(chrom="chr1", feature=["gene"] + ["exon"] * 3, start=[0, 0, 200, 400],
                            end=[500, 100, 300, 500], strand="+", gene_id="G1", gene_name="GENE1",
                            gene_type="protein_coding", transcript_id=["", "T1", "T1", "T1"],
                            transcript_type="protein_coding"))
    with pytest.warns(UserWarning, match="no exon shared by at least 2 of its 1 multi-exon transcript "):
        assert gene_model(one, "GENE1", "chr1", (150, 250)) is None
    single = one[one.feature.eq("gene") | one.start.eq(0)]                        # one single-exon transcript
    with pytest.warns(UserWarning, match="has no multi-exon transcript"):
        assert gene_model(single, "GENE1", "chr1", (50, 60)) is None
