#!/usr/bin/env python3
"""
mem_sample.py — launch a command and sample its resident memory on an interval,
                so you see the memory *trajectory*, not just a single peak.

    Built to watch makeCAF and distinguish a leak (RSS climbing across events)
    from a bounded footprint (RSS ~flat). Writes a TSV of (elapsed_s, rss_kb)
    you can plot, and prints a summary — including the true peak (VmHWM) — at
    the end.

    Usage:
      mem_sample.py [-i SEC] [-o FILE] [--tree] -- <command> [args...]

      -i / --interval   sampling interval in seconds (default 0.5)
      -o / --output     trajectory TSV (default: mem_sample.<pid>.tsv)
      --tree            sum RSS of the whole process tree, not just the launched
                        pid (use if makeCAF spawns children)

    Example — one run per build, then compare:
      mem_sample.py -o base.tsv -- makeCAF -c job.fcl ...
      mem_sample.py -o snap.tsv -- makeCAF -c job.fcl ...

    Standard library only (reads /proc), so Linux-only — which is where makeCAF
    runs anyway. No numpy/uproot needed.

    Written by Claude Code (Anthropic, Claude Opus 4.8) at J. Wolcott's request.
"""

import argparse
import os
import signal
import subprocess
import sys
import time


def _status_field_kb(pid, key):
    """Value (kB) of a /proc/<pid>/status line like 'VmRSS:', or None if gone."""
    try:
        with open("/proc/%d/status" % pid) as f:
            for line in f:
                if line.startswith(key):
                    return int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError, ValueError):
        return None
    return None


def descendants(root):
    """root and every pid beneath it, via /proc PPID links."""
    children = {}
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open("/proc/%s/stat" % entry) as f:
                after_comm = f.read().rsplit(")", 1)[1].split()
            ppid = int(after_comm[1])          # fields after "(comm)": state, ppid, ...
        except (FileNotFoundError, ProcessLookupError, IndexError, ValueError):
            continue
        children.setdefault(ppid, []).append(int(entry))
    out, stack = [], [root]
    while stack:
        p = stack.pop()
        out.append(p)
        stack.extend(children.get(p, []))
    return out


def total_rss_kb(root, tree):
    if not tree:
        return _status_field_kb(root, "VmRSS:")
    total, seen = 0, False
    for p in descendants(root):
        r = _status_field_kb(p, "VmRSS:")
        if r is not None:
            total += r
            seen = True
    return total if seen else None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--interval", type=float, default=0.5)
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("--tree", action="store_true")
    ap.add_argument("command", nargs=argparse.REMAINDER)
    args = ap.parse_args()

    cmd = args.command
    if cmd and cmd[0] == "--":
        cmd = cmd[1:]
    if not cmd:
        ap.error("no command given; put it after `--`, e.g. `-- makeCAF -c job.fcl`")

    proc = subprocess.Popen(cmd)
    pid = proc.pid
    out_path = args.output or ("mem_sample.%d.tsv" % pid)

    # make Ctrl-C / TERM take makeCAF down with us
    def _forward(signum, _frame):
        try:
            proc.send_signal(signum)
        except ProcessLookupError:
            pass
    signal.signal(signal.SIGINT, _forward)
    signal.signal(signal.SIGTERM, _forward)

    samples = []
    peak_hwm = 0                 # true peak for a single process (VmHWM is monotonic)
    peak_sampled = 0             # fallback / tree-mode peak
    t0 = time.monotonic()
    with open(out_path, "w") as out:
        out.write("elapsed_s\trss_kb\n")
        while proc.poll() is None:
            rss = total_rss_kb(pid, args.tree)
            if rss is not None:
                t = time.monotonic() - t0
                out.write("%.3f\t%d\n" % (t, rss))
                out.flush()
                samples.append((t, rss))
                peak_sampled = max(peak_sampled, rss)
                if not args.tree:
                    hwm = _status_field_kb(pid, "VmHWM:")
                    if hwm is not None:
                        peak_hwm = max(peak_hwm, hwm)
            time.sleep(args.interval)

    rc = proc.wait()
    dur = time.monotonic() - t0

    def mb(kb):
        return kb / 1024.0

    w = sys.stderr.write
    w("\n=== mem_sample summary ===\n")
    w("  command exit code : %d\n" % rc)
    w("  wall time         : %.1f s\n" % dur)
    w("  samples           : %d @ %ss\n" % (len(samples), args.interval))
    w("  trajectory TSV    : %s\n" % out_path)
    if samples:
        peak = peak_hwm if (peak_hwm and not args.tree) else peak_sampled
        label = "VmHWM" if (peak_hwm and not args.tree) else "sampled"
        first_rss, last_rss = samples[0][1], samples[-1][1]
        w("  peak RSS (%-7s): %d kB (%.1f MB)\n" % (label, peak, mb(peak)))
        w("  first / last RSS  : %d / %d kB  (Δ %+d kB)\n"
          % (first_rss, last_rss, last_rss - first_rss))
        if dur > 0:
            slope = (last_rss - first_rss) / dur
            w("  mean slope        : %+.1f kB/s  "
              "(steadily rising ⇒ likely leak; ~flat ⇒ bounded)\n" % slope)
    sys.exit(rc)


if __name__ == "__main__":
    main()