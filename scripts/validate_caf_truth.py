#!/usr/bin/env python3
"""
Validate CAF truth-particle parent/ancestor bookkeeping in one file.

This is intentionally a PyROOT-side diagnostic script: it opens a CAF file,
walks rec.mc.nu[*].sec[*], and checks that every secondary can be resolved
through the stored G4 parent chain to the primary referenced by ancestor_id.
It can also warn when secondary start positions are not near the immediate
parent's endpoint, which is a useful sanity check for inelastic/cascade-like
secondary production.

Examples
--------
# Validate every event in a CAF file
./validate_caf_truth.py file.root

# Look at one entry with verbose per-particle diagnostics
./validate_caf_truth.py file.root --event 12 --verbose

# Treat endpoint mismatches as failures instead of warnings
./validate_caf_truth.py file.root --strict-endpoint --endpoint-tol 0.05
"""

import argparse
import math
import os
import sys
from collections import defaultdict


TYPE_UNKNOWN = 0
TYPE_PRIMARY = 1
TYPE_PREFSI = 2
TYPE_SECONDARY = 3

TYPE_NAMES = {
    TYPE_UNKNOWN: "unknown",
    TYPE_PRIMARY: "primary",
    TYPE_PREFSI: "prefsi",
    TYPE_SECONDARY: "secondary",
}


parser = argparse.ArgumentParser(
    description="Validate CAF truth secondary parent chains and ancestor_id references"
)
parser.add_argument("caf_file", help="Input CAF ROOT file")
parser.add_argument("-t", "--tree-name", default="cafTree", help="Tree name (default: cafTree)")
parser.add_argument("-r", "--rec-tree", default="rec", help="Record branch name; use empty string for flat CAFs")
parser.add_argument("-l", "--sr-lib-name", default="libduneanaobj_StandardRecord_dict.so", help="StandardRecord dictionary library name")
parser.add_argument("-L", "--sr-lib-dir", help="Colon-separated directories to prepend to LD_LIBRARY_PATH")
parser.add_argument("-e", "--event", type=int, action="append", help="Specific tree entry to validate (repeatable)")
parser.add_argument("-n", "--nevents", type=int, default=-1, help="Maximum number of entries to scan")
parser.add_argument("--endpoint-tol", type=float, default=0.1, help="Parent endpoint vs child start tolerance in cm (default: 0.1)")
parser.add_argument("--strict-endpoint", action="store_true", help="Count endpoint mismatches as errors instead of warnings")
parser.add_argument("--allow-cross-ixn", action="store_true", help="Allow ancestor_id.ixn to point to a different interaction")
parser.add_argument("--max-report", type=int, default=50, help="Maximum issues to print before suppressing details")
parser.add_argument("--verbose", action="store_true", help="Print per-secondary trace details")
args = parser.parse_args()


try:
    import ROOT
    from ROOT import TFile, gSystem
except ImportError:
    print("ERROR: PyROOT is not available in this Python environment.", file=sys.stderr)
    print("Set up your ROOT/duneanaobj environment first, then rerun.", file=sys.stderr)
    sys.exit(2)


if args.sr_lib_dir:
    old = os.environ.get("LD_LIBRARY_PATH", "")
    os.environ["LD_LIBRARY_PATH"] = f"{args.sr_lib_dir}:{old}" if old else args.sr_lib_dir

if gSystem.Load(args.sr_lib_name) < 0:
    print("ERROR: Could not load StandardRecord library", file=sys.stderr)
    print(f"  library: {args.sr_lib_name}", file=sys.stderr)
    print(f"  LD_LIBRARY_PATH: {os.environ.get('LD_LIBRARY_PATH', '')}", file=sys.stderr)
    sys.exit(2)


def open_tree(path):
    tf = TFile.Open(path, "READ")
    if not tf or tf.IsZombie():
        raise RuntimeError(f"Could not open file: {path}")
    tree = tf.Get(args.tree_name)
    if not tree:
        raise RuntimeError(f"Could not load tree '{args.tree_name}' from file: {path}")
    return tf, tree


def get_record_root(tree):
    return getattr(tree, args.rec_tree) if args.rec_tree else tree


def safe_size(container):
    try:
        return int(container.size())
    except Exception:
        return len(container)


def safe_at(container, idx):
    try:
        if idx < 0:
            return None
        if hasattr(container, "size"):
            if idx >= container.size():
                return None
            return container.at(idx)
        if idx >= len(container):
            return None
        return container[idx]
    except Exception:
        return None


def vec3_tuple(v):
    return (float(v.x), float(v.y), float(v.z))


def dist(a, b):
    if a is None or b is None:
        return float("nan")
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def is_finite_vec(v):
    return v is not None and all(math.isfinite(x) for x in v)


def type_name(v):
    return TYPE_NAMES.get(int(v), str(int(v)))


def collection_for_type(ixn, ptype):
    ptype = int(ptype)
    if ptype == TYPE_PRIMARY:
        return ixn.prim, "prim"
    if ptype == TYPE_PREFSI:
        return ixn.prefsi, "prefsi"
    if ptype == TYPE_SECONDARY:
        return ixn.sec, "sec"
    return None, "unknown"


def iter_collection(coll):
    for i in range(safe_size(coll)):
        yield i, coll.at(i) if hasattr(coll, "at") else coll[i]


def build_g4_index(ixn):
    """Return G4ID -> (collection name, index, particle) for one interaction."""
    out = {}
    duplicates = defaultdict(list)
    for coll_name in ("prim", "prefsi", "sec"):
        coll = getattr(ixn, coll_name)
        for i, part in iter_collection(coll):
            g4id = int(part.G4ID)
            if g4id in out:
                duplicates[g4id].append((coll_name, i))
            else:
                duplicates[g4id].append((coll_name, i))
                out[g4id] = (coll_name, i, part)
    return out, duplicates


def issue(reports, level, entry, ixn_idx, coll, part_idx, message):
    reports.append((level, entry, ixn_idx, coll, part_idx, message))
    if len(reports) <= args.max_report:
        loc = f"entry {entry} ixn {ixn_idx} {coll}[{part_idx}]"
        print(f"{level}: {loc}: {message}")
    elif len(reports) == args.max_report + 1:
        print(f"... suppressing further issue details after --max-report={args.max_report}")


def resolve_ancestor_ref(root, this_ixn_idx, sec, reports, entry, sec_idx):
    anc = sec.ancestor_id
    anc_ixn = int(anc.ixn)
    anc_type = int(anc.type)
    anc_part = int(anc.part)

    if anc_ixn < 0 or anc_ixn >= safe_size(root.mc.nu):
        issue(reports, "ERROR", entry, this_ixn_idx, "sec", sec_idx,
              f"ancestor_id.ixn={anc_ixn} is outside rec.mc.nu range")
        return None

    if anc_ixn != this_ixn_idx and not args.allow_cross_ixn:
        issue(reports, "ERROR", entry, this_ixn_idx, "sec", sec_idx,
              f"ancestor_id.ixn={anc_ixn} does not match containing interaction")

    anc_ixn_obj = root.mc.nu.at(anc_ixn)
    anc_coll, anc_coll_name = collection_for_type(anc_ixn_obj, anc_type)
    if anc_coll is None:
        issue(reports, "ERROR", entry, this_ixn_idx, "sec", sec_idx,
              f"ancestor_id.type={anc_type} ({type_name(anc_type)}) is not a particle collection")
        return None

    anc_part_obj = safe_at(anc_coll, anc_part)
    if anc_part_obj is None:
        issue(reports, "ERROR", entry, this_ixn_idx, "sec", sec_idx,
              f"ancestor_id points to {anc_coll_name}[{anc_part}], but collection size is {safe_size(anc_coll)}")
        return None

    return anc_ixn, anc_coll_name, anc_part, anc_part_obj


def trace_to_primary(sec, g4_index):
    """Trace a secondary's CAF parent chain until a primary is reached.

    Returns (status, trace, primary_tuple). trace contains tuples
    (G4ID, collection, index, parent_G4ID). primary_tuple is
    (collection, index, particle), normally for collection 'prim'.
    """
    current_parent = int(sec.parent)
    trace = [(int(sec.G4ID), "sec", None, current_parent)]
    visited = {int(sec.G4ID)}

    while current_parent >= 0:
        if current_parent in visited:
            return "cycle", trace, None
        visited.add(current_parent)

        found = g4_index.get(current_parent)
        if found is None:
            return "missing-parent", trace, None

        coll_name, idx, part = found
        parent = int(part.parent)
        trace.append((int(part.G4ID), coll_name, idx, parent))
        if coll_name == "prim":
            return "ok", trace, found
        current_parent = parent

    return "no-primary", trace, None


def validate_entry(root, entry, reports, counters):
    n_ixn = safe_size(root.mc.nu)
    for ixn_idx in range(n_ixn):
        ixn = root.mc.nu.at(ixn_idx)
        g4_index, duplicates = build_g4_index(ixn)

        for g4id, places in duplicates.items():
            if len(places) > 1:
                counters["duplicate_g4id"] += 1
                issue(reports, "ERROR", entry, ixn_idx, "ixn", -1,
                      f"G4ID {g4id} appears multiple times: {places}")

        children_by_parent = defaultdict(list)
        for sec_idx, sec in iter_collection(ixn.sec):
            children_by_parent[int(sec.parent)].append((sec_idx, sec))

        for sec_idx, sec in iter_collection(ixn.sec):
            counters["secondaries"] += 1
            anc_ref = resolve_ancestor_ref(root, ixn_idx, sec, reports, entry, sec_idx)
            status, trace, primary = trace_to_primary(sec, g4_index)

            if args.verbose:
                trace_txt = " -> ".join(
                    f"{coll}[{idx}] G4ID={g4id} parent={parent}"
                    for g4id, coll, idx, parent in trace
                )
                if anc_ref is None:
                    anc_txt = "unresolved"
                else:
                    _, anc_coll, anc_idx, anc_part = anc_ref
                    anc_txt = f"{anc_coll}[{anc_idx}] G4ID={int(anc_part.G4ID)} pdg={int(anc_part.pdg)}"
                print(f"TRACE: entry {entry} ixn {ixn_idx} sec[{sec_idx}] G4ID={int(sec.G4ID)}: {trace_txt}; ancestor={anc_txt}; status={status}")

            if status != "ok":
                counters[f"trace_{status}"] += 1
                issue(reports, "ERROR", entry, ixn_idx, "sec", sec_idx,
                      f"could not trace parent chain to a primary: status={status}, trace={trace}")
                continue

            primary_coll, primary_idx, primary_part = primary
            if anc_ref is None:
                continue

            _, anc_coll, anc_idx, anc_part = anc_ref
            if anc_coll != "prim":
                issue(reports, "ERROR", entry, ixn_idx, "sec", sec_idx,
                      f"ancestor_id resolves to {anc_coll}[{anc_idx}], expected prim[{primary_idx}] from parent chain")
            elif anc_idx != primary_idx or int(anc_part.G4ID) != int(primary_part.G4ID):
                issue(reports, "ERROR", entry, ixn_idx, "sec", sec_idx,
                      f"ancestor_id resolves to prim[{anc_idx}] G4ID={int(anc_part.G4ID)}, "
                      f"but parent chain reaches prim[{primary_idx}] G4ID={int(primary_part.G4ID)}")
            else:
                counters["ancestor_ok"] += 1

            parent_tuple = g4_index.get(int(sec.parent))
            if parent_tuple is not None:
                parent_coll, parent_idx, parent_part = parent_tuple
                parent_end = vec3_tuple(parent_part.end_pos)
                child_start = vec3_tuple(sec.start_pos)
                d = dist(parent_end, child_start)
                if is_finite_vec(parent_end) and is_finite_vec(child_start) and d > args.endpoint_tol:
                    counters["endpoint_mismatch"] += 1
                    level = "ERROR" if args.strict_endpoint else "WARN"
                    issue(reports, level, entry, ixn_idx, "sec", sec_idx,
                          f"child start is {d:.6g} cm from immediate parent endpoint "
                          f"({parent_coll}[{parent_idx}] G4ID={int(parent_part.G4ID)}); "
                          f"parent_end={parent_end}, child_start={child_start}")
                else:
                    counters["endpoint_ok"] += 1

        # A compact parent-level endpoint check for inelastic/cascade-like vertices:
        # if a particle has multiple secondary children, their starts should be
        # mutually clustered near the parent's endpoint. This is warning-only
        # unless --strict-endpoint is set, same as the per-child check above.
        for parent_g4id, children in children_by_parent.items():
            if parent_g4id < 0 or len(children) < 2 or parent_g4id not in g4_index:
                continue
            parent_coll, parent_idx, parent_part = g4_index[parent_g4id]
            parent_end = vec3_tuple(parent_part.end_pos)
            far = []
            for child_idx, child in children:
                d = dist(parent_end, vec3_tuple(child.start_pos))
                if math.isfinite(d) and d > args.endpoint_tol:
                    far.append((child_idx, int(child.G4ID), d))
            if far:
                counters["multi_child_endpoint_mismatch"] += 1
                level = "ERROR" if args.strict_endpoint else "WARN"
                issue(reports, level, entry, ixn_idx, parent_coll, parent_idx,
                      f"parent G4ID={parent_g4id} has {len(children)} secondary children; "
                      f"{len(far)} are farther than {args.endpoint_tol} cm from endpoint: {far[:8]}")


def main():
    try:
        tf, tree = open_tree(args.caf_file)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    nentries = int(tree.GetEntries())
    if args.nevents >= 0:
        nentries = min(nentries, args.nevents)

    if args.event:
        entries = args.event
    else:
        entries = list(range(nentries))

    reports = []
    counters = defaultdict(int)

    for entry in entries:
        if entry < 0 or entry >= int(tree.GetEntries()):
            print(f"WARN: skipping out-of-range entry {entry}", file=sys.stderr)
            continue
        tree.GetEntry(entry)
        validate_entry(get_record_root(tree), entry, reports, counters)

    n_error = sum(1 for r in reports if r[0] == "ERROR")
    n_warn = sum(1 for r in reports if r[0] == "WARN")

    print("\n=== validation summary ===")
    print(f"file: {args.caf_file}")
    print(f"entries checked: {len(entries)}")
    print(f"secondaries checked: {counters['secondaries']}")
    print(f"ancestor chain matches: {counters['ancestor_ok']}")
    print(f"endpoint checks passed: {counters['endpoint_ok']}")
    if counters.get("endpoint_mismatch"):
        print(f"endpoint mismatches: {counters['endpoint_mismatch']} ({'errors' if args.strict_endpoint else 'warnings'})")
    if counters.get("multi_child_endpoint_mismatch"):
        print(f"multi-child endpoint mismatches: {counters['multi_child_endpoint_mismatch']} ({'errors' if args.strict_endpoint else 'warnings'})")
    for key in sorted(k for k in counters if k.startswith("trace_") or k == "duplicate_g4id"):
        print(f"{key}: {counters[key]}")
    print(f"errors: {n_error}")
    print(f"warnings: {n_warn}")

    if n_error:
        print("FAIL")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
