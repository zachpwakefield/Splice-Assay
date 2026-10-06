"""An agent's narrative of a probe (opt-in: `probe --agent-summary`, `splice-assay summarize PROBE_DIR`).

Every probe writes agent_prompt.md: what to write, the package's rules for reading a probe and what not to claim
(AGENT_GUIDE.md), and the probe's aggregate results (its report and, for the best-ranked events, the statistics of
their notable cohorts). Nothing patient-level goes in: no page's plotted values, no input table.

`summarize` hands that file to Claude Code (`claude -p` with no tools, run from an empty folder), checks the numbers
the narrative quotes against the results, and writes agent_summary.md (the narrative, marked as machine-written, with
the check's result) and agent_summary.json (the agent, model, date, prompt hash and cost). Any other agent can be
given agent_prompt.md by hand.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd

from . import __version__
from .config import Settings
from .dataset import InputError
from .events import PSI_MEANING

PROMPT, SUMMARY, RECORD = "agent_prompt.md", "agent_summary.md", "agent_summary.json"
TOP = 15                    # best-ranked events given in detail
ROWS = 12                   # notable cohorts shown per event
INSTRUCTION = "Write the summary that the document on standard input asks for. Reply with the summary only."
RESULTS = "## The probe's report"        # where the results start: the number check reads from here on
DETAIL = "## The best-ranked events in detail"


def rules(s: Settings, expression: bool = True, adjusted: bool = True) -> str:
    """The package's rules for reading a probe (AGENT_GUIDE.md), with the run's thresholds. `expression`: the Cox
    models include host expression; `adjusted`: there is an adjusted model."""
    a = f"{s.alpha:g}"
    host = ("The Cox models adjust for it (the report's notes name the fits where it was left out); the KM split "
            "does not." if expression else "The Cox models here do not include it, so any of the tests may echo it.")
    return "\n".join([
        "## Rules for reading a probe", "",
        "- A probe ranks candidates; it does not test a hypothesis, and every p is nominal. Say how many tests were "
        f"run and how many p < {a} chance alone gives (the report's \"How much is chance\"), and weigh the hits "
        "against that.",
        "- An event is worth a closer look when several of these hold: its hits clearly exceed chance; the "
        "significant cohorts agree in direction; " + ("the association holds in the adjusted model (compare the "
                                                       "adjusted and the base fit); " if adjusted else "")
        + "the same cohort shows a case-vs-reference hit, ideally within patients (paired).",
        f"- Treat an event with caution when: a single cohort has p just under {a} among many tested; directions "
        "flip between cohorts; " + ("the association disappears after adjustment; " if adjusted else "")
        + f"its fit is marked ‡ (fewer than {s.cox_low_power_events} events: low power) or † (proportional-hazards "
        f"test p < {s.ph_note_below:g}: the HR averages an effect that changes over follow-up); the report's notes "
        "flag it (a narrow PSI range, " + ("host expression left out, " if expression else "")
        + "few events per model term).",
        *(["- The adjusted model keeps only patients with every clinical variable recorded. Where it kept fewer than "
           "the base model (\"n … of …\" in its column), a p that rises may reflect the patients lost rather than "
           "the adjustment."] if adjusted else []),
        "- Host-gene expression: where the gene's own expression is prognostic in the same cohort, a splicing "
        f"association may echo it. {host}",
        "- Events that are strongly correlated (the report's correlation section) or that share exons may be one "
        "piece of evidence, not several. Two alternative first (or last) exons of one gene sum to 1, so their ρ near "
        "−1 is built in, not a finding; with three or more, a pair's ρ is not fixed by construction.",
        "- Protein changes are suggestions read from annotation, never measurements.",
        "- These are observational cohorts. Say \"is associated with\" or \"a candidate\"; never call an event a "
        "biomarker, and never say it predicts, drives or causes anything.", ""])


def task(adjusted: bool = True) -> str:
    """What to write, in a form the number check can follow. `adjusted`: there is an adjusted model."""
    return "\n".join([
        "## What to write", "", "Markdown, at most about 350 words, with exactly these four sections:", "",
        "### Bottom line",
        "Two or three sentences: is anything here worth a closer look, and how strong is it against chance?", "",
        "### Candidates",
        "Up to five events, best first, one list item each: the event by its label and its gene, the cohort or "
        "cohorts by name, " + ("the model (adjusted or base), " if adjusted else "") + "the endpoint, the HR with "
        "its 95% CI, p and q, and in one clause why it is or is not convincing.", "",
        "### Caveats", "The rules above that apply to these candidates.", "",
        "### Next checks",
        "What to look at next: which pages, the expression page, the gene map, or which analysis to run.", "",
        "Copy every number exactly as it is printed below, for example \"HR 1.42 (1.10–1.83), p 0.003, q 0.04\", and "
        "name its event (label and gene), " + ("model, " if adjusted else "") + "and cohort in the same paragraph or "
        "list item (for the host gene's own HR: the gene, \"expression\" and the cohort). Do not compute new "
        "numbers, round them differently or average them across cohorts: a script checks the numbers you quote "
        "against the results.", ""])


def _fmt_hr(hr, lo, hi) -> str:
    if hr is None or not np.isfinite(hr):
        return "–"
    ci = f" ({lo:.2f}–{hi:.2f})" if lo is not None and hi is not None and np.isfinite(lo) and np.isfinite(hi) else ""
    return f"HR {hr:.2f}{ci}"


def _flag(v) -> bool:
    """A True/False cell that may be missing (NaN, None) or a numpy bool."""
    return False if v is None or (isinstance(v, float) and np.isnan(v)) else bool(v)


def _num(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def _event_table(r, cells: pd.DataFrame, gex: pd.DataFrame | None, s: Settings) -> list[str]:
    """One event's notable cohorts: p < alpha in a survival test, a case-vs-reference hit, or its best cohort."""
    from .plot.style import fd, fp
    u = s.psi_hr_unit
    d = cells[cells.event_id.eq(r.event_id)]
    has_adj = "adj_cox_status" in d
    d = d.assign(_p=d.get("adj_cox_p", d.cox_p).fillna(d.cox_p)).sort_values("_p", na_position="last", kind="stable")
    with np.errstate(invalid="ignore"):
        hit = (d.get("adj_cox_p", pd.Series(np.nan, index=d.index)) < s.alpha) | (d.cox_p < s.alpha) | \
            (d.km_p < s.alpha) | d.group_hit.map(_flag).astype(bool) | d.cohort.eq(r.best_cohort)
    shown, rest = d[hit], d[~hit]
    head = ["| Cohort | " + ("Adjusted Cox | " if has_adj else "") + "Base Cox | KM | Case vs reference | "
            + ("ρ with expression | " if "expr_rho" in d else "") + (f"{r.host} expression Cox |" if gex is not None
                                                                         else ""),
            "|---|" + ("---|" if has_adj else "") + "---|---|---|" + ("---|" if "expr_rho" in d else "")
            + ("---|" if gex is not None else "")]
    rows = []
    for x in shown.head(ROWS).itertuples():
        def fit(prefix):
            st = getattr(x, f"{prefix}cox_status", None)
            if st != "tested":
                return f"not fitted ({st})" if isinstance(st, str) and st else "–"
            marks = (" ‡" if _flag(getattr(x, f"{prefix}cox_low_power", False)) else "") + (
                " †" if _num(getattr(x, f"{prefix}ph_p", np.nan)) < s.ph_note_below else "")
            q = _num(getattr(x, f"{prefix}cox_q", np.nan))
            n_a, n_b = _num(getattr(x, "adj_cox_n", np.nan)), _num(getattr(x, "cox_n", np.nan))
            fewer = prefix and x.cox_status == "tested" and np.isfinite(n_a) and np.isfinite(n_b) and n_a < n_b
            return (_fmt_hr(_num(getattr(x, f"{prefix}hr_per_{u}")), _num(getattr(x, f"{prefix}ci_low_{u}")),
                            _num(getattr(x, f"{prefix}ci_high_{u}")))
                    + f", p {fp(_num(getattr(x, f'{prefix}cox_p')))}" + (f", q {fp(q)}" if np.isfinite(q) else "")
                    + (f", n {n_a:.0f} of {n_b:.0f}" if fewer else "") + marks)
        km = f"p {fp(x.km_p)}" if x.km_status == "tested" else "–"
        grp = "; ".join(f"{design} Δ {fd(_num(getattr(x, f'{design}_delta_median')))}, "
                        f"p {fp(_num(getattr(x, f'{design}_p')))}"
                        + (" (hit)" if _flag(getattr(x, f"{design}_hit", False)) else "")
                        for design in ("paired", "unpaired") if getattr(x, f"{design}_status", None) == "tested")
        cols = [x.cohort] + ([fit("adj_")] if has_adj else []) + [fit(""), km, grp or "–"]
        if "expr_rho" in d:
            rho = _num(x.expr_rho)
            cols.append(f"ρ {rho:.2f}".replace("-", "−") if np.isfinite(rho) else "–")
        if gex is not None:
            g = gex[gex.cohort.eq(x.cohort)]
            g = g.iloc[0] if len(g) else None
            cols.append(f"{_fmt_hr(_num(g.hr_per_sd), _num(g.ci_low_sd), _num(g.ci_high_sd))}, p {fp(_num(g.cox_p))}"
                        if g is not None and g.cox_status == "tested" else "–")
        rows.append("| " + " | ".join(cols) + " |")
    more = len(shown) - min(len(shown), ROWS)
    t_a = rest.adj_cox_status.eq("tested") if has_adj else pd.Series(False, index=rest.index)
    t_b = rest.cox_status.eq("tested")
    n_tested = int((t_a | t_b | rest.km_status.eq("tested")).sum())
    hr = (rest[f"adj_hr_per_{u}"].where(t_a, rest[f"hr_per_{u}"]) if has_adj else rest[f"hr_per_{u}"])[t_a | t_b]
    tail = [f"Other cohorts tested: {n_tested}" + (f" (HR above 1 in {int((hr > 1).sum())}" + (
        f" of the {len(hr)} with a Cox fit" if len(hr) < n_tested else "") + ")" if len(hr) else "")
        + (f"; {more} more notable cohorts are in cells.csv" if more else "") + "."]
    return head + rows + [""] + tail


def prompt(report: str, ranked: pd.DataFrame, cells: pd.DataFrame, s: Settings, gex_cells: pd.DataFrame | None = None,
           top: int = TOP, expression: bool = True) -> str:
    """The text of agent_prompt.md: the instructions, the rules, the probe's report and the notable cohorts of its
    `top` best-ranked measurable events. `expression`: the Cox models include host expression."""
    U = s.psi_hr_unit.upper()
    adjusted = "adj_cox_status" in cells
    best = ranked[ranked.measurable].head(top) if len(ranked) else ranked
    L = ["# Summarise this splice-assay probe", "",
         f"<!-- Written by splice-assay {__version__}. `splice-assay summarize PROBE_DIR` gives this file to Claude "
         "Code with no tools. It holds the probe's report and aggregate statistics only. Give it to any other agent "
         "without access to the probe's folder, whose pages/*.csv hold per-patient values. -->",
         "",
         "You are writing a short narrative summary of a splice-assay probe for the researcher who ran it. Use only "
         "this document; there is nothing else to read.", "", rules(s, expression, adjusted), task(adjusted), "---", "",
         RESULTS, "",
         "\n".join("##" + x if x.startswith("#") else x for x in report.strip().splitlines()), "", "---", "",
         DETAIL, "",
         f"For each of the {len(best)} best-ranked measurable events: the cohorts with p < {s.alpha:g} in a survival "
         f"test ({'adjusted Cox, base Cox or KM' if adjusted else 'Cox or KM'}), a case-vs-reference hit, or its "
         f"best cohort. HRs are per {U} of the event's value (expression: per SD); q is Benjamini–Hochberg within "
         f"the gene, as in the report. ‡ fewer than {s.cox_low_power_events} events (low power); † "
         f"proportional-hazards test p < {s.ph_note_below:g}"
         + ("; \"n … of …\": the adjusted fit kept fewer patients than the base fit" if adjusted else "")
         + ". \"Other cohorts tested\" counts the cohorts not listed, and how many of their HRs ("
         + ("adjusted where fitted, else base" if adjusted else "Cox") + ") are above 1.", ""]
    for r in best.itertuples():
        gx = None
        if gex_cells is not None and len(gex_cells) and "gene" in gex_cells:
            g = gex_cells[gex_cells.gene.eq(r.host)]
            gx = g if len(g) else None
        meaning = PSI_MEANING.get(str(r.event_type).upper(), "PSI")
        L += [f"### {r.rank}. {r.label} ({r.gene}) · {r.event_type}, {meaning}", ""]
        L += _event_table(r, cells, gx, s) + [""]
    return "\n".join(L).rstrip() + "\n"


def write_prompt(out_dir, report: str, ranked: pd.DataFrame, cells: pd.DataFrame, s: Settings, hosts: dict,
                 gex_cells: pd.DataFrame | None = None, expression: bool = True) -> Path:
    """Write agent_prompt.md into a probe folder. `hosts`: event ID -> host gene (Dataset.events.expression_gene)."""
    ranked = ranked.assign(host=ranked.event_id.map(hosts).fillna(ranked.gene))
    path = Path(out_dir) / PROMPT
    path.write_text(prompt(report, ranked, cells, s, gex_cells, expression=expression), encoding="utf-8")
    return path


# ------------------------------------------------------------------------------------------------ the number check
_NUM = re.compile(r"\d*\.\d+(?:[eE][−+-]?\d+)?")
_P = r"(\d*\.\d+(?:[eE][−+-]?\d+)?)"
_HR = re.compile(r"(\d*\.\d+)\s*\(\s*(?:95\s*%\s*CI[:,]?\s*)?(\d*\.\d+)\s*(?:[–—-]|to)\s*(\d*\.\d+)\s*\)"
                 rf"(?:\s*[,;]?\s*[pP]\s*[=<]?\s*{_P})?(?:\s*[,;]?\s*[qQ]\s*[=<]?\s*{_P})?")  # HR (lo–hi), p …, q …
_POW = re.compile(r"(?<![\d.])(\d+(?:\.\d+)?|\.\d+)\s*[×x·]\s*10(?:\^\s*([−-]?\d+)|([⁻]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+))")
_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻−", "0123456789--")
_ESC = re.compile(r"\\([\\`*_{}\[\]()#+\-.!|<>])")


def _plain(text: str) -> str:
    """Markdown escapes dropped (COHORT\\_1 is COHORT_1), and powers of ten as e-notation (6.1×10⁻⁵, 6.1 x 10^-5)."""
    return _POW.sub(lambda m: f"{m.group(1)}e{(m.group(2) or m.group(3)).translate(_SUP)}", _ESC.sub(r"\1", text))


def numbers(text: str) -> list[tuple[str, float]]:
    """The decimal numbers in a text, as (printed, value): '1.42', '−0.82', '.42', '1.2e-5'; a minus counts only when
    it does not join two numbers (1.10-1.83 is a range)."""
    out = []
    digit = (lambda i: 0 <= i < len(text) and text[i].isdigit())                    # noqa: E731
    for m in _NUM.finditer(text):
        a, b = m.span()
        if (a and text[a - 1].isalnum()) or (a and text[a - 1] == "." and digit(a - 2)) or \
                (b < len(text) and text[b].isalnum()) or (b < len(text) and text[b] == "." and digit(b + 1)):
            continue                                     # part of a word, an ID or a version (0.2.0)
        neg = a and text[a - 1] in "−-" and not (a > 1 and text[a - 2].isalnum())
        printed = ("−" if neg else "") + m.group(0)
        out.append((printed, -float(m.group(0).replace("−", "-")) if neg else float(m.group(0).replace("−", "-"))))
    return out


def _found(v: float, have: set) -> bool:
    """A quoted value among the printed ones; a positive one also matches its negative ("a decrease of 0.24")."""
    return v in have or (v > 0 and -v in have)


def _float(x: str | None) -> float | None:
    return None if x is None else float(x.replace("−", "-"))


class _Fit(NamedTuple):
    values: tuple               # HR, CI low, CI high, p, q (None where not printed)
    kind: str                   # "adjusted", "base", "expression" or "splicing"; "" outside the detail tables
    label: str                  # the event's label ("" for expression)
    gene: str                   # the event's gene, or the gene of an expression fit
    cohort: str
    both: bool                  # its table has an adjusted and a base column


def _fits(results: str) -> list[_Fit]:
    """Every HR printed with a CI in the results; in the detail tables with its column, event, gene and cohort."""
    out, label, gene, head, detail = [], "", "", [], False
    for line in results.splitlines():
        detail = detail or line.startswith(DETAIL)
        m = re.match(r"### \d+\. (.+) \(([^()]*)\) · ", line) if detail else None
        if m:
            label, gene, head = m.group(1), m.group(2), []
        if detail and line.startswith("| Cohort |"):
            head = [c.strip() for c in line.strip().strip("|").split("|")]
        elif detail and head and line.startswith("| "):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            both = "Adjusted Cox" in head and "Base Cox" in head
            for h, c in zip(head[1:], cells[1:]):
                x = h.endswith(" expression Cox")
                kind = "expression" if x else {"Adjusted Cox": "adjusted", "Base Cox": "base"}.get(h, "splicing")
                out += [_Fit(tuple(map(_float, g.groups())), kind, "" if x else label,
                             h[:-len(" expression Cox")] if x else gene, cells[0], both) for g in _HR.finditer(c)]
        else:
            out += [_Fit(tuple(map(_float, g.groups())), "", "", "", "", False) for g in _HR.finditer(line)]
    return out


def _blocks(text: str) -> list[str]:
    """A narrative's paragraphs and top-level list items (with their nested and indented lines), each after the
    heading above it."""
    out, head, state = [], "", "new"                 # new: a block starts next; gap: after a blank line
    for line in text.splitlines():
        if not line.strip():
            state = "gap" if state == "in" else state
        elif re.match(r"#{1,6}\s", line):
            head, state = line, "new"
        elif state == "in" and not re.match(r"(?:[-*+]|\d+[.)])\s", line) or state == "gap" and line[:1].isspace():
            out[-1] += "\n" + line
            state = "in"
        else:
            out.append(head + "\n" + line)
            state = "in"
    return out


def _mentions(text: str, names) -> set:
    """The names a text mentions; one seen only inside a longer name (BRCA in BRCA_LumA) does not count."""
    spans = {n: [m.span() for m in re.finditer(rf"(?<![A-Za-z0-9]){re.escape(n)}(?![A-Za-z0-9])", text)]
             for n in names}

    def inside(a, b, n):
        return any(a2 <= a and b <= b2 and b2 - a2 > b - a for o, ss in spans.items() if o != n for a2, b2 in ss)
    return {n for n, ss in spans.items() if any(not inside(a, b, n) for a, b in ss)}


def _same(quoted: tuple, printed: tuple) -> bool:
    """An HR and CI, and the p and q when quoted, as one fit printed them, in that order."""
    return quoted[:3] == printed[:3] and all(a is None or a == b for a, b in zip(quoted[3:], printed[3:]))


def _attributed(f: _Fit, named: set, text: str) -> bool:
    """Whether a paragraph names a fit's cohort and event (label and gene; for expression the gene and the word),
    and does not name only the other model when its table has both."""
    if not f.kind:
        return True                                  # outside the detail tables: nothing to attribute to
    if f.cohort not in named or f.gene not in named:
        return False
    if f.kind == "expression":
        return "expression" in text.lower()
    if f.label not in named:
        return False
    adj = re.search(r"(?<!un)adjust", text, re.I) is not None
    base = re.search(r"\bbase\b|unadjusted", text, re.I) is not None
    return not (f.both and (f.kind == "adjusted" and base and not adj or f.kind == "base" and adj and not base))


class Checked(NamedTuple):
    quoted: int                 # decimal numbers the narrative quotes
    missing: list               # ... that the results do not print
    hrs: int                    # HRs quoted with a CI
    unmatched: list             # ... not printed so for the event and cohort (and model) named beside them


def check(narrative: str, results: str) -> Checked:
    """The narrative's numbers against the probe's results (agent_prompt.md from RESULTS on): the decimal numbers the
    results do not print, and each HR quoted with a CI that no fit prints with that CI (and the p and q quoted after
    it, in that order) for the event and cohort named in the same paragraph or list item, and for its model where the
    paragraph names only one. Whole numbers are not checked."""
    narrative = _plain(narrative)
    have = {v for _, v in numbers(results)}
    quoted = numbers(narrative)
    missing = [p for p, v in quoted if not _found(v, have)]
    fits = _fits(results)
    names = {x for f in fits for x in (f.label, f.gene, f.cohort)} - {""}
    hrs, unmatched = 0, []
    for block in _blocks(narrative):
        named = None
        for m in _HR.finditer(block):
            hrs += 1
            named = _mentions(block, names) if named is None else named
            q = tuple(map(_float, m.groups()))
            if not any(_same(q, f.values) and _attributed(f, named, block) for f in fits):
                unmatched.append(m.group(0).strip())
    return Checked(len(quoted), list(dict.fromkeys(missing)), hrs, list(dict.fromkeys(unmatched)))


def _note(c: Checked) -> str:
    """The check's result, for agent_summary.md."""
    def n(k, word):
        return f"{k} {word}{'' if k == 1 else 's'}"
    if not c.quoted:
        return "*Number check: the summary quotes no decimal numbers.*"
    if not (c.missing or c.unmatched):
        return (f"*Number check: every decimal number quoted ({c.quoted}) is in the probe's results" + (
            f", and every HR quoted with a CI ({c.hrs}) is printed with that CI, p and q for the event and cohort "
            "(and model, where named) beside it" if c.hrs else "") + ". Whole numbers are not checked, nor which event or cohort "
            "the other numbers belong to.*")
    out = []
    if c.missing:
        one = len(c.missing) == 1
        out.append(f"{n(len(c.missing), 'decimal number')} quoted (of {c.quoted}) {'is' if one else 'are'} not in "
                   f"the probe's results: {', '.join(c.missing)}. Treat the sentences that quote "
                   f"{'it' if one else 'them'} as wrong until you have checked {'it' if one else 'them'}.")
    if c.unmatched:
        one = len(c.unmatched) == 1
        out.append(f"{n(len(c.unmatched), 'HR')} quoted with a CI (of {c.hrs}) {'is' if one else 'are'} not printed "
                   f"with those numbers, in that order, for the event and cohort (and model, where named) beside "
                   f"{'it' if one else 'them'}: {'; '.join(c.unmatched)}. Check {'it' if one else 'them'} against "
                   "the tables.")
    return "*Number check: " + " ".join(out) + "*"


# ------------------------------------------------------------------------------------------------ running the agent
def command(model: str | None = None, executable: str = "claude") -> list[str]:
    """Claude Code in print mode: no tools, no MCP servers, nothing asked, nothing saved, one JSON result."""
    return [executable, "-p", INSTRUCTION, "--tools", "", "--strict-mcp-config", "--permission-mode", "dontAsk",
            "--no-session-persistence", "--output-format", "json"] + (["--model", model] if model else [])


def _result(stdout) -> dict:
    """The result object of `claude -p --output-format json`: the object itself or, in verbose mode, the last of a
    list of messages (on the last line when other lines precede it); {} when there is none."""
    text = str(stdout or "").strip()
    for chunk in dict.fromkeys((text, text.splitlines()[-1] if text else "")):
        try:
            out = json.loads(chunk or "{}")
            break
        except json.JSONDecodeError:
            continue
    else:
        return {}
    if isinstance(out, list):
        out = next((m for m in reversed(out) if isinstance(m, dict) and m.get("type") == "result"), {})
    return out if isinstance(out, dict) else {}


def summarize(probe_dir, model: str | None = None, timeout: float = 600, dry_run: bool = False, log=print,
              executable: str | None = None, runner=None) -> dict:
    """Have Claude Code write a narrative of a probe from its agent_prompt.md. Writes agent_summary.md and
    agent_summary.json beside the report and returns their paths (with dry_run: the command, and nothing is sent).
    `model`: passed to `claude --model` (default: the CLI's own). `runner`: subprocess.run, or a stand-in (tests)."""
    runner = runner or subprocess.run
    d = Path(probe_dir)
    src = d / PROMPT
    if not src.exists():
        raise InputError(f"{src}: no {PROMPT} (run `splice-assay probe` with this version to write it)")
    text = src.read_text(encoding="utf-8")
    exe = executable or shutil.which("claude")
    cmd = command(model, exe or "claude")
    what = (f"{PROMPT} ({len(text.encode()) / 1024:.0f} KB: the probe's report and aggregate statistics, no "
            "per-patient values)")
    if dry_run:
        log(f"would send {what} to Claude Code: " + " ".join(f'"{c}"' if not c or " " in c else c for c in cmd))
        return dict(command=cmd)
    if exe is None:
        raise InputError("Claude Code is not installed (no `claude` on PATH): see https://code.claude.com, or give "
                         f"{src} to another agent by hand")
    log(f"sending {what} to Claude Code ({Path(exe).name} -p, no tools); waiting for the whole answer, which can "
        f"take a few minutes (giving up after {timeout:g} s)")
    with tempfile.TemporaryDirectory() as empty:          # no project CLAUDE.md, settings or files to find
        try:
            proc = runner(cmd, input=text, capture_output=True, text=True, encoding="utf-8", errors="replace",
                          timeout=timeout, cwd=empty)
        except subprocess.TimeoutExpired:
            raise InputError(f"Claude Code gave no answer within {timeout:g} s") from None
        except OSError as e:
            raise InputError(f"Claude Code could not be run ({exe}): {e}") from None
    out = _result(proc.stdout)
    result = out.get("result")
    if proc.returncode != 0 or out.get("is_error") or not isinstance(result, str) or not result.strip():
        why = ((result.strip() if isinstance(result, str) else "")
               or str(proc.stderr or proc.stdout or "").strip()[-400:] or "no output")
        login = any(w in why.lower() for w in ("log in", "logged in", "/login", "api key", "authenticat", "oauth"))
        raise InputError(f"Claude Code failed (exit {proc.returncode}): {why}" + (
            " (sign in first: run `claude` in a terminal and use /login, or set ANTHROPIC_API_KEY)" if login else ""))
    checked = check(result, text.split(RESULTS, 1)[-1])
    models = sorted(out["modelUsage"]) if isinstance(out.get("modelUsage"), dict) else []
    version = ""
    try:
        version = str(runner([exe, "--version"], capture_output=True, text=True, encoding="utf-8",
                             errors="replace", timeout=30).stdout or "").strip()
    except (OSError, subprocess.SubprocessError):
        pass
    when = datetime.now(timezone.utc).replace(microsecond=0)
    first = text.split(RESULTS, 1)[-1].strip().splitlines()
    title = next((x[4:] for x in first if x.startswith("### splice-assay probe · ")), "")
    by = f"Claude Code ({', '.join(models)})" if models else "Claude Code"
    body = [f"<!-- Written by {by} from {PROMPT} on {when:%Y-%m-%d}. Machine-written and unreviewed: check it "
            "against report.md and the pages before you rely on it. -->", "",
            f"# Agent summary{' · ' + title.replace('splice-assay probe · ', '') if title else ''}", "",
            result.strip(), "", "---", "", _note(checked), ""]
    paths = dict(summary=d / SUMMARY, record=d / RECORD)
    paths["summary"].write_text("\n".join(body), encoding="utf-8")
    record = dict(splice_assay=__version__, agent="claude", agent_version=version, models=models,
                  command=cmd, created=when.isoformat(), prompt=PROMPT,
                  prompt_sha256=hashlib.sha256(text.encode()).hexdigest(), prompt_bytes=len(text.encode()),
                  numbers_quoted=checked.quoted, numbers_not_found=checked.missing, hrs_quoted=checked.hrs,
                  hrs_not_matched=checked.unmatched, cost_usd=out.get("total_cost_usd"),
                  duration_ms=out.get("duration_ms"))
    paths["record"].write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    flagged = len(checked.missing) + len(checked.unmatched)
    log(f"agent summary: {paths['summary']}" + (f" (the number check flags {flagged}: see the end of it)"
                                                  if flagged else ""))
    return paths
