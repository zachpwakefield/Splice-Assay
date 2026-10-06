import json
import subprocess

import pandas as pd
import pytest

import splice_assay as sa
from splice_assay import agent
from splice_assay.agent import PROMPT, check, numbers, summarize
from splice_assay.cli import main
from splice_assay.dataset import InputError
from splice_assay.probe import probe

FAST = sa.Settings(formats=("png",), dpi=60)


@pytest.fixture(scope="module")
def probe_dir(ds, tmp_path_factory):
    out = tmp_path_factory.mktemp("probe")
    probe(ds, genes=["SYN1"], settings=FAST, out_dir=out, max_pages=0, correlation=True, log=lambda *_: None)
    return out


def _fake(result="", code=0, is_error=False, seen=None, cost=0.01):
    """A stand-in for subprocess.run that answers like `claude -p --output-format json`."""
    def run(cmd, **kw):
        if "--version" in cmd:
            return subprocess.CompletedProcess(cmd, 0, stdout="2.1.0 (Claude Code)\n", stderr="")
        if seen is not None:
            seen.update(cmd=cmd, **kw)
        out = dict(type="result", is_error=is_error, result=result, total_cost_usd=cost, duration_ms=1200,
                   modelUsage={"claude-opus-5-5": {}})
        return subprocess.CompletedProcess(cmd, code, stdout=json.dumps(out), stderr="")
    return run


def test_the_probe_writes_what_an_agent_needs(ds, probe_dir):
    text = (probe_dir / PROMPT).read_text()
    for part in ("## Rules for reading a probe", "## What to write", "### Bottom line", "### Next checks",
                 "### splice-assay probe · SYN1 · OS", "#### How much is chance", "## The best-ranked events in detail",
                 "### 1. SE:1 (SYN1) · SE, PSI = inclusion", "| Cohort | Adjusted Cox | Base Cox |", "ρ with expression"):
        assert part in text, part
    assert not any(sid in text for sid in ds.samples.sample_id)                   # nothing patient-level
    assert "claude -p <" not in text and "pages/*.csv hold per-patient values" in text   # no run with file access
    assert "| `agent_prompt.md` |" in (probe_dir / "report.md").read_text()


def test_the_rules_follow_the_run(probe_dir):
    r = agent.rules(sa.Settings())
    assert "p < 0.05 chance alone gives" in r and "fewer than 20 events" in r and "The Cox models adjust for it" in r
    r = agent.rules(sa.Settings(alpha=0.01, cox_low_power_events=30, ph_note_below=0.01), expression=False,
                    adjusted=False)
    assert "p < 0.01 chance alone gives" in r and "fewer than 30 events" in r and "test p < 0.01" in r
    assert "do not include it" in r and "adjust" not in r
    cells, ranked = pd.read_csv(probe_dir / "cells.csv"), pd.read_csv(probe_dir / "events.csv")
    text = agent.prompt((probe_dir / "report.md").read_text(), ranked, cells, FAST, expression=False)
    assert "The Cox models here do not include it" in text and "in a survival test (adjusted Cox, base Cox or KM)" \
        in text


def test_event_table_marks_each_design_and_counts_the_other_cohorts(probe_dir):
    cells, ranked = pd.read_csv(probe_dir / "cells.csv"), pd.read_csv(probe_dir / "events.csv")
    r, u = next(ranked.assign(host="SYN1").itertuples()), FAST.psi_hr_unit
    d = cells[cells.event_id.eq(r.event_id)].copy()
    a, b, k = [c for c in d.cohort if c != r.best_cohort]
    best = d.cohort.eq(r.best_cohort)
    d.loc[best, ["paired_status", "unpaired_status", "adj_cox_status", "cox_status"]] = "tested"
    d.loc[best, ["paired_delta_median", "paired_p", "unpaired_delta_median", "unpaired_p"]] = [0.01, 0.5, 0.03, 0.01]
    d.loc[best, ["paired_hit", "unpaired_hit", "group_hit"]] = [False, True, True]       # an unpaired hit only
    d.loc[best, ["adj_cox_n", "cox_n"]] = [180, 240]
    d.loc[~best, ["adj_cox_p", "cox_p", "km_p", "group_hit"]] = [0.5, 0.5, 0.5, False]   # nothing notable
    d.loc[d.cohort.eq(a), ["adj_cox_status", "cox_status", "km_status", f"adj_hr_per_{u}"]] = \
        ["tested", "tested", "tested", 1.2]
    d.loc[d.cohort.eq(b), ["adj_cox_status", "cox_status", "km_status", f"hr_per_{u}"]] = \
        ["too_few_events", "tested", "tested", 0.8]                                  # a base fit only
    d.loc[d.cohort.eq(k), ["adj_cox_status", "cox_status", "km_status"]] = ["too_few_events", "too_few_events",
                                                                            "tested"]   # KM only
    rows = agent._event_table(r, d, None, FAST)
    row = next(x for x in rows if x.startswith(f"| {r.best_cohort} |")).split(" | ")
    assert row[4] == "paired Δ +0.01, p 0.50; unpaired Δ +0.03, p 0.010 (hit)"
    assert row[1].endswith(", n 180 of 240") or ", n 180 of 240 " in row[1]          # the adjusted fit's column
    assert rows[-1] == "Other cohorts tested: 3 (HR above 1 in 1 of the 2 with a Cox fit)."


RESULTS = """## The probe's report

### splice-assay probe · SYN1 · OS

Tests: 24; p < 0.05: 3.

---

## The best-ranked events in detail

### 1. SE:1 (SYN1) · SE, PSI = inclusion

| Cohort | Adjusted Cox | Base Cox | KM | Case vs reference | SYN1 expression Cox |
|---|---|---|---|---|---|
| COH1 | HR 1.36 (1.08–1.72), p 0.009, q 0.04 | HR 1.30 (1.01–1.66), p 0.041, q 0.12 | p 0.02 | paired Δ −0.236, p 0.01 (hit) | HR 0.90 (0.70–1.15), p 0.40 |

### 2. ALE:2 (SYN1) · ALE, PSI = use of the distal exon

| Cohort | Adjusted Cox | Base Cox | KM | Case vs reference | SYN1 expression Cox |
|---|---|---|---|---|---|
| COH2 | HR 0.70 (0.55–0.89), p 6.1e-5, q 0.03 | HR 0.75 (0.60–0.94), p 0.012, q 0.05 | p 0.03 | – | HR 1.12 (0.86–1.45), p 0.47 |

### 3. SE:1 (SYN2) · SE, PSI = inclusion

| Cohort | Adjusted Cox | Base Cox | KM | Case vs reference | SYN2 expression Cox |
|---|---|---|---|---|---|
| COH1_B | HR 2.46 (1.42–4.26), p 0.001, q 0.26 | HR 2.10 (1.30–3.39), p 0.002, q 0.20 | p 0.01 | – | – |
"""


def test_numbers_and_the_check():
    text = "HR 1.42 (1.10–1.83), p 0.003, q 0.04. ρ −0.82 and .42; range 1.10-1.83; v0.2.0; 1.2e-5"
    assert [p for p, _ in numbers(text)] == ["1.42", "1.10", "1.83", "0.003", "0.04", "−0.82", ".42", "1.10",
                                              "1.83", "1.2e-5"]
    assert check("HR 1.42, p 0.0031, 0.1 and −0.82", "HR 1.42 p 0.003 0.10 ρ −0.82") == (4, ["0.0031"], 0, [])
    assert check("ρ −0.82", "ρ 0.82").missing == ["−0.82"]                         # a minus must be printed
    assert check("a decrease of 0.82, p 6.1×10⁻⁵", "Δ −0.82, p 6.1e-5").missing == []


def test_the_number_check_matches_each_hr_to_its_event_and_cohort():
    ok = check("1. **SE:1 (SYN1)**, adjusted, OS\n   - COH1: HR 1.36 (1.08–1.72), p 0.009, q 0.04\n"
               "- ALE:2 (SYN1) in COH2: HR 0.70 (0.55–0.89), p 6.1×10⁻⁵; SYN1 expression in COH2: "
               "HR 1.12 (0.86–1.45).\n\nPSI fell (a decrease of 0.236).\n\n#### SE:1 (SYN2)\n"
               "- COH1\\_B (base model): HR 2.10 (1.30–3.39), P = 0.002\n"     # a markdown escape, a capital P
               "| SE:1 (SYN1) | COH1 | HR 1.30 (1.01–1.66), p 0.041 |", RESULTS)  # a table; no model named
    assert ok.missing == [] and ok.hrs == 5 and ok.unmatched == []
    bad = check("- SE:1 (SYN1) in COH2: HR 1.36 (1.08–1.72), p 0.009.\n"          # another cohort's
                "- ALE:2 (SYN1) in COH1: HR 1.36 (1.08–1.72).\n"                   # another event's
                "- SE:1 (SYN1) in COH1: HR 1.36 (1.08–1.72), p 0.041.\n"           # the base model's p
                "- SE:1 (SYN1) in COH1: HR 1.42 (1.10–1.83), p 0.003, q 0.04.\n"   # the instructions' example
                "- SE:1 (SYN1) in COH1: HR 2.46 (1.42–4.26), p 0.001.\n"           # SE:1 of SYN2's
                "- SE:1 (SYN2) in COH1: HR 2.46 (1.42–4.26), p 0.001.\n"           # COH1_B's, not COH1's
                "- SE:1 (SYN1) in COH1: HR 1.36 (1.08–1.72), p 0.04, q 0.009.\n"   # p and q swapped
                "- SE:1 (SYN1) in COH1, adjusted: HR 0.90 (0.70–1.15), p 0.40.\n"  # the expression HR
                "- SE:1 (SYN1) in COH1, adjusted: HR 1.30 (1.01–1.66), p 0.041.", RESULTS)   # the base HR
    assert bad.hrs == 9 and len(bad.unmatched) == 8 and bad.missing == ["1.10", "1.83", "0.003"]   # each text once


def test_summarize_writes_the_narrative_and_checks_its_numbers(probe_dir):
    text = (probe_dir / PROMPT).read_text()
    real = numbers(text.split("## The best-ranked events in detail")[1])[0][0]
    seen = {}
    paths = summarize(probe_dir, executable="/x/claude", runner=_fake(
        f"### Bottom line\nSE:1 has {real}, and HR 9.99 nowhere.", seen=seen), log=lambda *_: None)
    assert seen["input"] == text and seen["cmd"][:3] == ["/x/claude", "-p", agent.INSTRUCTION]
    assert {"--tools", "", "--strict-mcp-config", "--no-session-persistence"} <= set(seen["cmd"])
    assert seen["cwd"] != str(probe_dir)                                           # an empty folder
    md = paths["summary"].read_text()
    assert "Machine-written and unreviewed" in md and "# Agent summary · SYN1 · OS" in md
    assert "Number check: 1 decimal number quoted (of 2) is not in the probe's results: 9.99. Treat the sentences " \
        "that quote it" in md
    rec = json.loads(paths["record"].read_text())
    assert rec["models"] == ["claude-opus-5-5"] and rec["numbers_not_found"] == ["9.99"] and rec["cost_usd"] == 0.01
    assert rec["agent_version"] == "2.1.0 (Claude Code)" and len(rec["prompt_sha256"]) == 64
    assert rec["hrs_quoted"] == 0 and rec["hrs_not_matched"] == []
    ok = summarize(probe_dir, executable="/x/claude", runner=_fake(f"Only {real}."), log=lambda *_: None)
    assert "Number check: every decimal number quoted (1) is in the probe's results. Whole numbers" in \
        ok["summary"].read_text()
    none = summarize(probe_dir, executable="/x/claude", runner=_fake("No numbers."), log=lambda *_: None)
    assert "Number check: the summary quotes no decimal numbers." in none["summary"].read_text()


def test_summarize_failures(probe_dir, tmp_path, monkeypatch):
    with pytest.raises(InputError, match="sign in first"):
        summarize(probe_dir, executable="/x/claude", runner=_fake("Not logged in · Please run /login", code=1),
                  log=lambda *_: None)
    with pytest.raises(InputError, match="failed"):
        summarize(probe_dir, executable="/x/claude", runner=_fake("boom", is_error=True), log=lambda *_: None)
    monkeypatch.setattr(agent.shutil, "which", lambda _: None)
    with pytest.raises(InputError, match="not installed"):
        summarize(probe_dir, log=lambda *_: None)
    with pytest.raises(InputError, match="no agent_prompt.md"):
        summarize(tmp_path, log=lambda *_: None)

    def slow(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, kw.get("timeout"))
    with pytest.raises(InputError, match="no answer within 5 s"):
        summarize(probe_dir, executable="/x/claude", runner=slow, timeout=5, log=lambda *_: None)

    def broken(cmd, **kw):
        raise PermissionError(13, "Permission denied")
    with pytest.raises(InputError, match="could not be run"):
        summarize(probe_dir, executable="/x/claude", runner=broken, log=lambda *_: None)
    with pytest.raises(InputError, match="sign in first"):
        summarize(probe_dir, executable="/x/claude", runner=_fake("OAuth token has expired", is_error=True),
                  log=lambda *_: None)

    def raw(stdout):                                     # output that is not one result object
        return lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, stdout="2.1.0" if "--version" in cmd else stdout,
                                                             stderr="")
    for out in ("[]", "null", '"text"', json.dumps({"result": {"text": "x"}}), "not json"):
        with pytest.raises(InputError, match="failed"):
            summarize(probe_dir, executable="/x/claude", runner=raw(out), log=lambda *_: None)
    verbose = json.dumps([{"type": "system"}, {"type": "result", "is_error": False, "result": "Fine."}])
    assert "Fine." in summarize(probe_dir, executable="/x/claude", runner=raw(verbose),
                                log=lambda *_: None)["summary"].read_text()


def test_cli_summarize_and_probe_agent_summary(probe_dir, tmp_path, capsys, monkeypatch):
    assert main(["summarize", str(probe_dir), "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "would send agent_prompt.md" in out and "no per-patient values" in out and '--tools ""' in out
    assert main(["summarize", str(tmp_path)]) == 2
    assert "no agent_prompt.md" in capsys.readouterr().err
    from splice_assay import example
    data = example.write(tmp_path / "ex")["samples"].parent
    monkeypatch.setattr(agent.shutil, "which", lambda _: "/x/claude")
    monkeypatch.setattr(agent.subprocess, "run", _fake("### Bottom line\nNothing stands out."))
    (tmp_path / "fast.json").write_text('{"formats": ["png"], "dpi": 60}')
    assert main(["probe", str(data), "--gene", "SYN2", "--no-adjust", "--max-pages", "0", "--agent-summary",
                 "--out", str(tmp_path / "pr"), "--settings", str(tmp_path / "fast.json")]) == 0
    assert "Nothing stands out." in (tmp_path / "pr" / "agent_summary.md").read_text()
    monkeypatch.setattr(agent.subprocess, "run", _fake("Not logged in · Please run /login", code=1))
    capsys.readouterr()
    assert main(["probe", str(data), "--gene", "SYN2", "--no-adjust", "--max-pages", "0", "--agent-summary",
                 "--endpoint", "all", "--out", str(tmp_path / "all"), "--settings", str(tmp_path / "fast.json")]) == 1
    out, err = capsys.readouterr()
    assert err.count("agent summary not written") == 1 and "sign in first" in err and "splice-assay summarize" in err
    assert sorted(p.name for p in (tmp_path / "all").iterdir()) == ["DSS", "OS"]  # every endpoint is still probed
    assert all((tmp_path / "all" / e / "report.md").exists() for e in ("DSS", "OS")) and "  1. " in out
