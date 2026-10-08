#!/usr/bin/env python3
"""
tile_spine_events.py — make a longer SPINE HDF5 input by tiling its trigger
                       and event lists N times, to stress-test per-event memory.

    makeCAF's record loop is driven by the *trigger* dataset (one record per
    trigger), and each trigger is paired to an event by position — so to get
    N x the records you must tile BOTH the trigger and events datasets, in
    lockstep. Two wrinkles this handles:

      * The events dataset holds region references into the per-product
        datasets, not the data itself. This makes a byte copy of the file,
        which keeps every referenced dataset's identity, so the duplicated
        refs stay valid and point at the original product regions.
      * Records are matched/deduped by trigger timestamp, so identical copies
        would collapse back together. Each trigger copy therefore gets its
        time field bumped by k*offset (k = copy index) to keep the copies
        distinct. The read path keys off position and refs, not time, so this
        doesn't change which event a trigger reads.

    Usage:
      tile_spine_events.py -f 100 in.h5 out.h5
      tile_spine_events.py -f 100 --trigger-dataset trigger --time-field time_s \\
                           --time-offset 86400 in.h5 out.h5
      tile_spine_events.py -f 100 --skip-trigger in.h5 out.h5   # events only

    The events dataset is auto-detected by its region-reference fields; the
    trigger dataset defaults to 'trigger'. Needs h5py + numpy. Start with -f 2
    and confirm makeCAF makes 2*N records without erroring before scaling up.

    Written by Claude Code (Anthropic, Claude Opus 4.8) at J. Wolcott's request.
"""

import argparse
import shutil
import sys

import h5py
import numpy as np


def has_region_refs(dt):
    if dt.names:
        return any(h5py.check_dtype(ref=dt[n]) is h5py.RegionReference for n in dt.names)
    return h5py.check_dtype(ref=dt) is h5py.RegionReference


def find_ref_datasets(f):
    hits = []

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and has_region_refs(obj.dtype):
            hits.append(name)

    f.visititems(visit)
    return hits


def tile_dataset(f, name, factor, forge_field=None, forge_offset=0):
    """Rewrite dataset `name` as `factor` circular copies of its rows. If
    forge_field is given, add k*forge_offset to that field in copy k so the
    copies stay distinguishable. Returns (old_len, new_len)."""
    ds = f[name]
    old = ds[...]
    n0 = old.shape[0]
    dt = ds.dtype
    attrs = dict(ds.attrs)

    tiled = np.concatenate([old] * factor)
    if forge_field is not None:
        if dt.names is None or forge_field not in dt.names:
            sys.exit("ERROR: no field '%s' in dataset '%s' (fields: %s)"
                     % (forge_field, name, dt.names))
        col = tiled[forge_field].copy()
        for k in range(factor):
            col[k * n0:(k + 1) * n0] = old[forge_field] + k * forge_offset
        tiled[forge_field] = col

    del f[name]
    nd = f.create_dataset(name, shape=(tiled.shape[0],), dtype=dt)
    nd[...] = tiled
    for k, v in attrs.items():
        nd.attrs[k] = v
    return n0, tiled.shape[0]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("-f", "--factor", type=int, required=True,
                    help="number of circular copies to make")
    ap.add_argument("--dataset", default=None,
                    help="events dataset path (default: auto-detect by region-ref fields)")
    ap.add_argument("--trigger-dataset", default="trigger",
                    help="trigger dataset path (default: 'trigger')")
    ap.add_argument("--time-field", default="time_s",
                    help="trigger time field to forge (default: 'time_s')")
    ap.add_argument("--time-offset", type=int, default=86400,
                    help="value added per copy to keep triggers distinct (default: 86400)")
    ap.add_argument("--skip-trigger", action="store_true",
                    help="tile only the events dataset (no trigger tiling/forging)")
    args = ap.parse_args()

    if args.factor < 1:
        ap.error("--factor must be >= 1")

    print("copying %s -> %s ..." % (args.input, args.output))
    shutil.copyfile(args.input, args.output)

    with h5py.File(args.output, "r+") as f:
        # events: the region-reference dataset
        ev_name = args.dataset
        if ev_name is None:
            hits = find_ref_datasets(f)
            if not hits:
                sys.exit("ERROR: no region-ref dataset found; pass --dataset.")
            if len(hits) > 1:
                sys.exit("ERROR: multiple ref datasets %s; pass --dataset." % hits)
            ev_name = hits[0]
        if ev_name not in f:
            sys.exit("ERROR: events dataset '%s' not in file." % ev_name)
        ev0, ev1 = tile_dataset(f, ev_name, args.factor)
        print("events  '%s' : %d -> %d rows" % (ev_name, ev0, ev1))

        # trigger: drives the record count; forge its time so copies stay distinct
        if not args.skip_trigger:
            tn = args.trigger_dataset
            if tn not in f:
                sys.exit("ERROR: trigger dataset '%s' not in file; pass "
                         "--trigger-dataset or --skip-trigger." % tn)
            tg0, tg1 = tile_dataset(f, tn, args.factor,
                                    forge_field=args.time_field,
                                    forge_offset=args.time_offset)
            print("trigger '%s' : %d -> %d rows  (%s += k*%d)"
                  % (tn, tg0, tg1, args.time_field, args.time_offset))
            if ev0 != tg0:
                print("  WARNING: events and trigger had different original lengths "
                      "(%d vs %d); the position pairing may be off." % (ev0, tg0))

    print("done: %s" % args.output)


if __name__ == "__main__":
    main()
