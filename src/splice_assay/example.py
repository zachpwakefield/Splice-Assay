"""A synthetic dataset for trying the package and for its tests. Everything here is simulated.

Three genes on made-up coordinates: SYN1 (chr7, plus strand) with a cassette exon (SE), an alternative 3' splice site
(A3SS) and a retained intron hosting a snoRNA (RI); SYN2 (chr12, minus strand) with mutually exclusive exons (MXE);
SYN3 (chr3, plus strand) with two alternative first exons (AFE), two alternative last exons (ALE) and the HIT index of
two internal exons, as HITindex reports them (write() also writes those as HITindex matrices, for import-hitindex).
SYN3's values come from their own random stream, after everything else, so SYN1 and SYN2 do not change.
Four cohorts: COH1 (40 pairs, planted tumour-normal and survival effects), COH2 (6 pairs plus unpaired normals, so no
within-patient test), COH3 (24 pairs, effects on other events) and COH4 (no normals). A clinical table holds age
(which raises the hazard), sex and stage (stage II and III raise it) for the Cox covariate example. A protein
cache gives protein-coding versions of the transcripts for the protein band.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

COHORTS = {"COH1": dict(tumours=180, pairs=40, extra_normals=0),
           "COH2": dict(tumours=120, pairs=6, extra_normals=14),
           "COH3": dict(tumours=150, pairs=24, extra_normals=0),
           "COH4": dict(tumours=80, pairs=0, extra_normals=0)}

EVENTS = pd.DataFrame([
    dict(event_id="SYN1:SE:1", gene="SYN1", gene_id="SYNG0001", chrom="chr7", strand="+", event_type="SE",
         constant="102000-102100;104000-104300", variable="103000-103150", label="SE:1"),
    dict(event_id="SYN1:A3SS:1", gene="SYN1", gene_id="SYNG0001", chrom="chr7", strand="+", event_type="A3SS",
         constant="100000-100200;102000-102100", variable="101970-102000", label="A3SS:1"),
    dict(event_id="SYN1:RI:1", gene="SYN1", gene_id="SYNG0001", chrom="chr7", strand="+", event_type="RI",
         constant="103000-103150;104000-104300", variable="103150-104000", label="RI:1"),
    dict(event_id="SYN2:MXE:1", gene="SYN2", gene_id="SYNG0002", chrom="chr12", strand="-", event_type="MXE",
         constant="200000-200400;205500-206000", variable="204000-204120;203000-203090", label="MXE:1"),
])

# logit-scale baseline, tumour shift per cohort, and Cox log-HR per +0.10 PSI per cohort
EFFECTS = {
    "SYN1:SE:1": dict(base=0.3, shift=dict(COH1=1.0, COH2=0.9, COH3=0.0, COH4=0.5),
                      loghr=dict(COH1=0.30, COH2=0.15, COH3=0.0, COH4=0.0)),
    "SYN1:A3SS:1": dict(base=1.2, shift=dict(COH1=-0.7, COH2=0.0, COH3=-0.5, COH4=0.0),
                        loghr=dict(COH1=0.0, COH2=0.0, COH3=-0.2, COH4=0.0)),
    "SYN1:RI:1": dict(base=-2.8, shift=dict(COH1=1.6, COH2=0.3, COH3=0.0, COH4=0.0),
                      loghr=dict(COH1=0.2, COH2=0.0, COH3=0.0, COH4=0.0)),
    "SYN2:MXE:1": dict(base=0.8, shift=dict(COH1=0.0, COH2=0.0, COH3=0.9, COH4=0.0),
                       loghr=dict(COH1=0.0, COH2=0.0, COH3=0.25, COH4=0.0)),
}

# (gene, gene_id, chrom, strand, gene span, {transcript: (type, exons)}, nested snoRNAs); 0-based half-open
GENES = [
    ("SYN1", "SYNG0001", "chr7", "+", (100000, 106000), {
        "SYN1-201": ("lncRNA", [(100000, 100200), (102000, 102100), (103000, 103150), (104000, 104300), (105500, 106000)]),
        "SYN1-202": ("lncRNA", [(100000, 100200), (102000, 102100), (104000, 104300), (105500, 106000)]),
        "SYN1-203": ("lncRNA", [(100000, 100200), (101970, 102100), (103000, 103150), (104000, 104300), (105500, 105900)]),
        "SYN1-204": ("retained_intron", [(100000, 100200), (102000, 103150)]),
    }, [("SNORD900", (103500, 103580))]),
    ("SYN3", "SYNG0003", "chr3", "+", (300000, 310000), {
        "SYN3-201": ("lncRNA", [(300000, 300200), (303000, 303100), (305000, 305150), (308000, 308300)]),
        "SYN3-202": ("lncRNA", [(301000, 301150), (303000, 303100), (305000, 305150), (309000, 309400)]),
        "SYN3-203": ("lncRNA", [(300000, 300200), (303000, 303100), (309000, 309400)]),
        "SYN3-204": ("lncRNA", [(301000, 301150), (305000, 305150), (308000, 308300)]),
    }, []),
    ("SYN2", "SYNG0002", "chr12", "-", (200000, 206000), {
        "SYN2-201": ("lncRNA", [(200000, 200400), (204000, 204120), (205500, 206000)]),
        "SYN2-202": ("lncRNA", [(200000, 200400), (203000, 203090), (205500, 206000)]),
        "SYN2-203": ("lncRNA", [(200100, 200400), (204000, 204120), (205500, 205950)]),
        "SYN2-204": ("lncRNA", [(200050, 200400), (203000, 203090), (205500, 206000)]),
    }, [("SNORA900", (202000, 202130))]),
]


def gtf_text() -> str:
    """A small GENCODE-style GTF of the two synthetic genes (1-based, inclusive coordinates)."""
    out = ["##description: splice-assay synthetic annotation"]

    def line(chrom, feat, a, b, strand, attrs):
        return "\t".join([chrom, "SYNTH", feat, str(a + 1), str(b), ".", strand, ".",
                          " ".join(f'{k} "{v}";' for k, v in attrs.items())])

    for gene, gid, chrom, strand, (g0, g1), txs, nested in GENES:
        ga = dict(gene_id=f"{gid}.1", gene_type="lncRNA", gene_name=gene)
        out.append(line(chrom, "gene", g0, g1, strand, ga))
        for tid, (ttype, exons) in txs.items():
            ta = dict(ga, transcript_id=f"{tid}.1", transcript_type=ttype, transcript_name=tid)
            out.append(line(chrom, "transcript", exons[0][0], exons[-1][1], strand, ta))
            for i, (a, b) in enumerate(exons if strand == "+" else exons[::-1], start=1):
                out.append(line(chrom, "exon", a, b, strand, dict(ta, exon_number=i)))
        for k, (name, (a, b)) in enumerate(nested, start=1):
            na = dict(gene_id=f"{gid}S{k}.1", gene_type="snoRNA", gene_name=name)
            out.append(line(chrom, "gene", a, b, strand, na))
            ta = dict(na, transcript_id=f"{gid}S{k}T.1", transcript_type="snoRNA", transcript_name=name)
            out.append(line(chrom, "transcript", a, b, strand, ta))
            out.append(line(chrom, "exon", a, b, strand, dict(ta, exon_number=1)))
    return "\n".join(out) + "\n"


# A synthetic protein cache for the same genes (see protein.py): protein-coding versions of the transcripts, with
# proteins translated from a made-up sequence that has no stop codon in any frame (no T on the plus strand, no A on
# the minus strand), and a few made-up features. 0-based half-open exons; CDS from CDS_START (1-based) to near the end.
PROTEIN_TX = {
    "SYN1": ("SYNG0001", "chr7", "+", 100101, {
        "SYN1-201": ("protein_coding", "1", [(100000, 100200), (102000, 102100), (103000, 103150), (104000, 104300),
                                             (105500, 106000)]),
        "SYN1-202": ("protein_coding", "1", [(100000, 100200), (102000, 102100), (104000, 104300), (105500, 106000)]),
        "SYN1-203": ("protein_coding", "2", [(100000, 100200), (101970, 102100), (103000, 103150), (104000, 104300),
                                             (105500, 106000)]),
        "SYN1-205": ("retained_intron", "2", [(100000, 100200), (102000, 102100), (103000, 104300), (105500, 106000)]),
    }),
    "SYN2": ("SYNG0002", "chr12", "-", 205900, {
        "SYN2-201": ("protein_coding", "1", [(200000, 200400), (204000, 204120), (205500, 206000)]),
        "SYN2-202": ("protein_coding", "2", [(200000, 200400), (203000, 203090), (205500, 206000)]),
    }),
}
# (gene, database, name, genomic span 1-based) -> drawn on every coding transcript that encodes part of the span
PROTEIN_FEATURES = [
    ("SYN1", "interpro", "Synthetic kinase domain", (104031, 104290)),
    ("SYN1", "mobidblite", "mobidb-lite", (103041, 103120)),
    ("SYN1", "elm", "LIG_SYN_1", (103061, 103075)),
    ("SYN1", "interpro", "Synthetic repeat", (100121, 100200)),
    ("SYN2", "interpro", "Synthetic Ig-like domain", (203001, 205900)),
    ("SYN2", "tmhmm", "TMhelix", (200301, 200360)),
]
_CODE = {a + b + c: aa for (a, b, c), aa in zip(((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
                                                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def protein_cache(seed: int = 11) -> dict[str, pd.DataFrame]:
    """exons, proteins and features tables of a protein cache for the synthetic genes."""
    rng = np.random.default_rng(seed)
    ex_rows, pr_rows, ft_rows = [], [], []
    for gene, (gid, chrom, strand, cds_start, txs) in PROTEIN_TX.items():
        lo = min(a for _, _, ex in txs.values() for a, _ in ex)
        hi = max(b for _, _, ex in txs.values() for _, b in ex)
        sense = dict(zip(range(lo + 1, hi + 1), rng.choice(list("ACG" if strand == "+" else "TGC"), hi - lo)))
        for tid, (ttype, tsl, exons) in txs.items():
            order = exons if strand == "+" else exons[::-1]
            pos = [x for a, b in order for x in (range(a + 1, b + 1) if strand == "+" else range(b, a, -1))]
            coding = ttype == "protein_coding"
            cds = []
            if coding:
                cds = pos[pos.index(cds_start):len(pos) - 120]
                cds = cds[:len(cds) - len(cds) % 3]
            rel = {g: i + 1 for i, g in enumerate(cds)}
            for k, (a, b) in enumerate(order, start=1):
                inside = [g for g in range(a + 1, b + 1) if g in rel]
                ex_rows.append(dict(transcript_id=tid, transcript_name=tid, transcript_type=ttype,
                                    transcript_support_level=tsl, gene_id=gid, gene_name=gene, chr=chrom,
                                    strand=strand, start=a + 1, end=b, exon_id=f"{tid}-E{k}", exon_number=k,
                                    protein_id=f"{tid}-P" if coding else "",
                                    cds_gen_start=min(inside) if inside else np.nan,
                                    cds_gen_stop=max(inside) if inside else np.nan,
                                    cds_rel_start=min(rel[g] for g in inside) if inside else np.nan,
                                    cds_rel_stop=max(rel[g] for g in inside) if inside else np.nan))
            if not coding:
                continue
            nt = "".join(sense[g] for g in cds)
            pr_rows.append(dict(transcript_id=tid, protein_id=f"{tid}-P",
                                protein_seq="".join(_CODE[nt[i:i + 3]] for i in range(0, len(nt), 3))))
            for g_, db, name, (g0, g1) in PROTEIN_FEATURES:
                hit = sorted(rel[g] for g in range(g0, g1 + 1) if g in rel) if g_ == gene else []
                if hit:
                    ft_rows.append(dict(transcript_id=tid, database=db, feature_id=name, name=name,
                                        start=(hit[0] + 2) // 3, stop=(hit[-1] + 2) // 3, g_start=g0, g_end=g1))
    return dict(exons=pd.DataFrame(ex_rows), proteins=pd.DataFrame(pr_rows), features=pd.DataFrame(ft_rows))


# SYN3's HITindex events, as import-hitindex numbers them: (type, exon 0-based half-open, other first/last exon)
HIT_EVENTS = [("AFE", (300000, 300200), (301000, 301150)), ("AFE", (301000, 301150), (300000, 300200)),
              ("ALE", (308000, 308300), (309000, 309400)), ("ALE", (309000, 309400), (308000, 308300)),
              ("HIT", (303000, 303100), None), ("HIT", (305000, 305150), None)]


def _hitindex_events() -> pd.DataFrame:
    rows, n = [], {}
    for kind, (a, b), other in HIT_EVENTS:
        n[kind] = n.get(kind, 0) + 1
        rows.append(dict(event_id=f"SYN3:{kind}:{n[kind]:04d}", gene="SYN3", gene_id="SYNG0003", chrom="chr3",
                         strand="+", event_type=kind, constant="" if other is None else f"{other[0]}-{other[1]}",
                         variable=f"{a}-{b}", label=f"{kind}:{n[kind]:04d}",
                         source_id=f"SYNG0003.1;chr3:{a + 1}-{b};{kind.lower()}"))
    return pd.DataFrame(rows)


def _hitindex_values(samples, survival, seed) -> tuple[pd.DataFrame, pd.DataFrame]:
    """SYN3 values (long) and expression, from their own random stream: AFE:0001 higher in COH1 tumours; ALE:0001
    lower in COH3 tumours; the HIT index of HIT:0001 rises by 0.35 in COH1 tumours; HIT:0002 goes with shorter OS.
    The second first and last exons take the rest of the usage (1 - the first, plus noise)."""
    rng = np.random.default_rng(seed + 303)
    tum = samples.group.eq("tumour").to_numpy()
    coh = samples.cohort.to_numpy()
    n = len(samples)
    u = rng.normal(0, 0.6, n)
    logistic = (lambda z: 1 / (1 + np.exp(-z)))
    afe = logistic(0.4 + u + rng.normal(0, 0.4, n) + 1.0 * (tum & (coh == "COH1")))
    ale = logistic(0.2 + rng.normal(0, 0.5, n) - 0.8 * (tum & (coh == "COH3")))
    hit1 = -0.05 + 0.08 * rng.normal(size=n) + 0.35 * (tum & (coh == "COH1"))
    os_ = survival[survival.endpoint.eq("OS")].set_index("patient_id").time
    z = np.zeros(n)
    for c in np.unique(coh):                         # shorter OS -> higher HIT:0002 (standardised log time per cohort)
        m = (coh == c) & tum
        t = np.log(samples.patient_id[m].map(os_).to_numpy(float))
        ok = np.isfinite(t)
        if ok.sum() > 2:
            zz = np.zeros(m.sum())
            zz[ok] = (t[ok] - t[ok].mean()) / t[ok].std()
            z[np.flatnonzero(m)] = zz
    hit2 = 0.10 - 0.12 * z + 0.10 * rng.normal(size=n)
    vals = {"SYN3:AFE:0001": afe, "SYN3:AFE:0002": 1 - afe + rng.normal(0, 0.02, n),
            "SYN3:ALE:0001": ale, "SYN3:ALE:0002": 1 - ale + rng.normal(0, 0.02, n),
            "SYN3:HIT:0001": hit1, "SYN3:HIT:0002": hit2}
    rows = []
    for e, v in vals.items():
        lo = -1.0 if ":HIT:" in e else 0.0
        v = np.round(np.clip(v, lo, 1.0), 3)
        v[rng.uniform(size=n) < 0.03] = np.nan
        rows += [dict(event_id=e, sample_id=sid, psi=float(x)) for sid, x in zip(samples.sample_id, v)]
    expr = pd.DataFrame(dict(gene="SYN3", sample_id=samples.sample_id,
                             value=np.round(5 + 0.8 * rng.normal(size=n), 4)))
    return pd.DataFrame(rows), expr


def hitindex_matrices(tables: dict) -> dict[str, pd.DataFrame]:
    """SYN3's values as HITindex matrices (afe, ale, hit): rows '<gene_id>;<chrom>:<start>-<end>;<kind>' (1-based),
    one column per sample."""
    ev = tables["events"].dropna(subset=["source_id"]) if "source_id" in tables["events"] else tables["events"].iloc[:0]
    ev = ev[ev.source_id.astype(str).ne("")]
    wide = tables["psi"][tables["psi"].event_id.isin(ev.event_id)].pivot(index="event_id", columns="sample_id",
                                                                           values="psi")
    wide = wide.reindex(index=ev.event_id, columns=tables["samples"].sample_id)
    wide.index = ev.source_id.to_numpy()
    kind = ev.event_type.str.lower().to_numpy()
    return {k: wide[kind == k].rename_axis("id") for k in ("afe", "ale", "hit")}


def make(seed: int = 7) -> dict[str, pd.DataFrame]:
    """samples, psi (long), events, survival (OS and DSS, days), expression (log2(TPM + 1), long) and clinical."""
    rng = np.random.default_rng(seed)
    srows = []
    for c, k in COHORTS.items():
        for i in range(k["tumours"]):
            pid = f"{c}-P{i:04d}"
            srows.append(dict(sample_id=f"{pid}-T", patient_id=pid, cohort=c, group="tumour"))
            if i < k["pairs"]:
                srows.append(dict(sample_id=f"{pid}-N", patient_id=pid, cohort=c, group="normal"))
        for i in range(k["extra_normals"]):
            pid = f"{c}-Q{i:04d}"
            srows.append(dict(sample_id=f"{pid}-N", patient_id=pid, cohort=c, group="normal"))
    samples = pd.DataFrame(srows)
    patients = samples.patient_id.unique()
    u = {e: dict(zip(patients, rng.normal(0, 0.7, len(patients)))) for e in EFFECTS}   # patient effect per event
    psi_rows = []
    for e, eff in EFFECTS.items():
        for r in samples.itertuples():
            z = eff["base"] + u[e][r.patient_id] + rng.normal(0, 0.4)
            if r.group == "tumour":
                z += eff["shift"][r.cohort]
            val = round(float(1 / (1 + np.exp(-z))), 3)
            if rng.uniform() < 0.03:
                val = np.nan
            psi_rows.append(dict(event_id=e, sample_id=r.sample_id, psi=val))
    psi = pd.DataFrame(psi_rows)
    wide = psi.pivot(index="sample_id", columns="event_id", values="psi")
    expr_rows = []
    for gene in ("SYN1", "SYN2"):
        for r in samples.itertuples():
            v = 5 + 0.8 * rng.normal() + 0.3 * u["SYN1:SE:1"][r.patient_id]
            expr_rows.append(dict(gene=gene, sample_id=r.sample_id, value=round(float(v), 4)))
    expression = pd.DataFrame(expr_rows)
    ex = expression[expression.gene.eq("SYN1")].set_index("sample_id").value
    surv_rows, clin_rows = [], []
    tum = samples[samples.group.eq("tumour")]
    means = {e: float(wide[e].mean()) for e in EFFECTS}
    for r in tum.itertuples():
        age = round(float(rng.normal(62, 11)), 1)
        sex = "female" if rng.uniform() < 0.5 else "male"
        stage = str(rng.choice(["I", "II", "III"], p=[0.45, 0.35, 0.2]))
        clin_rows.append(dict(patient_id=r.patient_id, age=age if rng.uniform() > 0.02 else np.nan, sex=sex,
                              stage=stage))
        if rng.uniform() < 0.05:                      # no follow-up
            continue
        lp = 0.2 * (ex[r.sample_id] - 5) / 0.8 + 0.03 * (age - 62) + {"I": 0.0, "II": 0.4, "III": 0.9}[stage]
        for e, eff in EFFECTS.items():
            x = wide.at[r.sample_id, e]
            if np.isfinite(x):
                lp += eff["loghr"][r.cohort] * (x - means[e]) / 0.10
        rate = np.exp(lp) / (7 * 365.25)
        t_event = rng.exponential(1 / rate)
        t_cens = rng.uniform(0.3, 12) * 365.25
        t = max(1.0, round(min(t_event, t_cens)))
        died = t_event <= t_cens
        surv_rows.append(dict(patient_id=r.patient_id, endpoint="OS", time=t, event=int(died)))
        disease = rng.uniform() < 0.75                # drawn for everyone, so the random stream never depends on outcomes
        surv_rows.append(dict(patient_id=r.patient_id, endpoint="DSS", time=t, event=int(died and disease)))
    survival = pd.DataFrame(surv_rows)
    hit_psi, hit_expr = _hitindex_values(samples, survival, seed)
    return dict(samples=samples, psi=pd.concat([psi, hit_psi], ignore_index=True),
                events=pd.concat([EVENTS, _hitindex_events()], ignore_index=True), survival=survival,
                expression=pd.concat([expression, hit_expr], ignore_index=True), clinical=pd.DataFrame(clin_rows))


def write(out_dir, seed: int = 7) -> dict[str, Path]:
    """Write the synthetic tables (CSV) and the GTF into out_dir/data/, a protein cache into out_dir/proteins/, and
    SYN3's AFE/ALE/HIT values as HITindex matrices into out_dir/hitindex/."""
    d = Path(out_dir) / "data"
    d.mkdir(parents=True, exist_ok=True)
    paths = {}
    tables = make(seed)
    for name, df in tables.items():
        paths[name] = d / f"{name}.csv"
        df.to_csv(paths[name], index=False, lineterminator="\n")
    hd = Path(out_dir) / "hitindex"                        # SYN3 again, as HITindex matrices (for import-hitindex)
    hd.mkdir(parents=True, exist_ok=True)
    for kind, m in hitindex_matrices(tables).items():
        m.to_csv(hd / f"{kind}.csv", lineterminator="\n")
    paths["hitindex"] = hd
    paths["gtf"] = d / "annotation.gtf"
    paths["gtf"].write_text(gtf_text())
    pc = Path(out_dir) / "proteins"
    pc.mkdir(parents=True, exist_ok=True)
    for name, df in protein_cache().items():
        df.to_csv(pc / f"{name}.csv", index=False, lineterminator="\n")
    paths["proteins"] = pc
    return paths
