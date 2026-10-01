"""Event geometry: the constant exons, the region PSI measures, and the junctions of each splice form.

Coordinates are 0-based and half-open (the BED and rMATS convention), written "start-end" and joined with ";".

    constant          the flanking constant exons (two for SE, RI, A3SS, A5SS and MXE)
    variable          the region PSI measures:
                        SE    the cassette exon
                        RI    the retained intron (it must span exactly the gap between the constant exons)
                        A3SS  the extension of the long form (between the two alternative sites)
                        A5SS  as A3SS
                        MXE   two exons: the one PSI measures first, then the other one
    psi_junctions     optional: introns "donor-acceptor" of the form PSI counts
    other_junctions   optional: introns of the other form

For SE, RI, A3SS, A5SS and MXE the junctions follow from the exons and are derived when the columns are empty. Any
other event_type is drawn from its constant and variable regions, with arcs only where junctions are given.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from .dataset import InputError

EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")
PSI_MEANING = {"SE": "PSI = inclusion", "RI": "PSI = retention", "A3SS": "PSI = long form",
               "A5SS": "PSI = long form", "MXE": "PSI = filled exon"}
_IV = re.compile(r"\s*(\d+)\s*-\s*(\d+)\s*")

Interval = tuple[int, int]


def parse_intervals(text, what: str = "interval") -> list[Interval]:
    """'100-200;300-400' -> [(100, 200), (300, 400)]; empty -> []."""
    if text is None or str(text).strip() in ("", "nan", "None"):
        return []
    out = []
    for part in str(text).split(";"):
        m = _IV.fullmatch(part)
        if not m:
            raise InputError(f"cannot read {what} {part.strip()!r}: expected start-end (0-based, half-open)")
        a, b = int(m.group(1)), int(m.group(2))
        if b <= a:
            raise InputError(f"{what} {part.strip()!r}: end must be greater than start")
        out.append((a, b))
    return out


def format_intervals(ivs) -> str:
    return ";".join(f"{a}-{b}" for a, b in ivs)


@dataclass(frozen=True)
class Geometry:
    event_type: str
    constant: tuple[Interval, ...]
    variable: tuple[Interval, ...]          # the region PSI measures
    other: tuple[Interval, ...]             # MXE: the exon of the other form (drawn outlined)
    psi_arcs: tuple[Interval, ...]          # junctions of the PSI form
    other_arcs: tuple[Interval, ...]        # junctions of the other form

    @property
    def span(self) -> Interval:
        pieces = self.constant + self.variable + self.other
        return min(p[0] for p in pieces), max(p[1] for p in pieces)

    @property
    def variable_length(self) -> int:
        return sum(b - a for a, b in self.variable)


def geometry(event_id: str, row) -> Geometry:
    """The geometry of one event from its events-table row (event_type, constant, variable, optional junctions)."""
    t = str(row["event_type"]).upper()
    if not t:
        raise InputError(f"event {event_id}: event_type is empty (needed to draw it)")
    const = sorted(parse_intervals(row["constant"], f"event {event_id} constant"))
    var = parse_intervals(row["variable"], f"event {event_id} variable")
    given_psi = parse_intervals(row.get("psi_junctions", ""), f"event {event_id} psi_junctions")
    given_oth = parse_intervals(row.get("other_junctions", ""), f"event {event_id} other_junctions")
    if not var:
        raise InputError(f"event {event_id}: the variable region is empty (needed to draw it)")

    def need(ok: bool, msg: str):
        if not ok:
            raise InputError(f"event {event_id} ({t}): {msg}")

    other: list[Interval] = []
    if t in EVENT_TYPES:
        n_var = 2 if t == "MXE" else 1
        need(len(const) == 2, f"expected 2 constant exons, got {len(const)}")
        need(len(var) == n_var, f"expected {n_var} variable region(s), got {len(var)}")
        (u, d) = const
        need(u[1] <= d[0], "the constant exons overlap")
    if t == "SE":
        c = var[0]
        need(u[1] <= c[0] and c[1] <= d[0], "the cassette exon must lie between the constant exons")
        psi, oth = [(u[1], c[0]), (c[1], d[0])], [(u[1], d[0])]
    elif t == "RI":
        need(var[0] == (u[1], d[0]), f"the retained intron must be {u[1]}-{d[0]}, the gap between the constant exons")
        psi, oth = [], [(u[1], d[0])]
    elif t in ("A3SS", "A5SS"):
        v = var[0]
        adjacent = [c for c in const if c[0] == v[1] or c[1] == v[0]]
        need(len(adjacent) == 1, "the variable region must extend exactly one constant exon (the short form)")
        short = adjacent[0]
        flank = const[1] if short == const[0] else const[0]
        long_ = (min(short[0], v[0]), max(short[1], v[1]))
        if flank[1] <= long_[0]:
            psi, oth = [(flank[1], long_[0])], [(flank[1], short[0])]
        elif long_[1] <= flank[0]:
            psi, oth = [(long_[1], flank[0])], [(short[1], flank[0])]
        else:
            raise InputError(f"event {event_id} ({t}): the flanking exon overlaps the alternative exon")
    elif t == "MXE":
        a, b = var
        need(max(a[0], b[0]) >= min(a[1], b[1]), "the two exons overlap")
        need(u[1] <= min(a[0], b[0]) and max(a[1], b[1]) <= d[0], "both exons must lie between the constant exons")
        psi, oth = [(u[1], a[0]), (a[1], d[0])], [(u[1], b[0]), (b[1], d[0])]
        var, other = [a], [b]
    elif t not in EVENT_TYPES:
        psi, oth = [], []
    psi = given_psi or psi
    oth = given_oth or oth
    return Geometry(t, tuple(const), tuple(var), tuple(other), tuple(psi), tuple(oth))
