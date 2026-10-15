import argparse
import hashlib
import json
import math
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def json_safe(value, vector_entries=None):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        if vector_entries is not None:
            value = value[:vector_entries]
        return json_safe(value.tolist(), vector_entries)
    if isinstance(value, np.generic):
        return json_safe(value.item(), vector_entries)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, (list, tuple)):
        if vector_entries is not None:
            value = value[:vector_entries]
        return [json_safe(item, vector_entries) for item in value]
    if isinstance(value, dict):
        return {
            str(key): json_safe(item, vector_entries)
            for key, item in value.items()
        }
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Dump readable ROOT tree branches to JSON using uproot and Awkward."
    )
    parser.add_argument("input", help="Input ROOT file")
    parser.add_argument(
        "-o",
        "--output",
        help="Output JSON file (default: INPUT.json)",
    )
    parser.add_argument(
        "-n",
        "--entries",
        type=int,
        default=1,
        help="Entries to dump per tree; use 0 to dump all entries (default: 1)",
    )
    parser.add_argument(
        "--vector-entries",
        type=int,
        help="Maximum elements to keep from each vector (default: all)",
    )
    parser.add_argument(
        "--fingerprint-unreadable",
        action="store_true",
        help="Add per-entry SHA-1 fingerprints for branches that cannot be decoded",
    )
    return parser.parse_args()


def record_fingerprints(branch, max_entries):
    target_entries = branch.num_entries
    if max_entries > 0:
        target_entries = min(max_entries, target_entries)

    fingerprints = []
    try:
        for basket_index in range(branch.num_baskets):
            if len(fingerprints) >= target_entries:
                break

            basket = branch.basket(basket_index)
            data = np.asarray(basket.data, dtype=np.uint8)
            offsets = basket.byte_offsets
            if offsets is None:
                basket_entries = basket.num_entries
                if basket_entries <= 0 or len(data) % basket_entries:
                    return None
                entry_size = len(data) // basket_entries
                if entry_size == 0:
                    return None
                offsets = range(0, len(data) + 1, entry_size)
            else:
                offsets = np.asarray(offsets).tolist()

            for start, stop in zip(offsets, offsets[1:]):
                if len(fingerprints) >= target_entries:
                    break
                payload = data[start:stop].tobytes()
                fingerprints.append(hashlib.sha1(payload).hexdigest())
    except Exception:
        return None

    if len(fingerprints) != target_entries:
        return None
    return fingerprints


def dump_root_file(
    input_path,
    output_path,
    max_entries,
    max_vector_entries,
    fingerprint_unreadable,
):
    results = {"trees": {}, "unreadable_branches": {}}
    if fingerprint_unreadable:
        results["raw_fingerprints"] = {
            "algorithm": "sha1",
            "trees": {},
            "unavailable": {},
        }
    readable_count = 0
    unreadable_count = 0

    with uproot.open(input_path) as root_file:
        for path, tree in root_file.items(recursive=True):
            if not hasattr(tree, "num_entries") or not hasattr(tree, "iteritems"):
                continue

            tree_name = str(path).split(";")[0]
            entry_stop = tree.num_entries
            if max_entries > 0:
                entry_stop = min(max_entries, entry_stop)

            branches = {}
            unreadable = {}
            for branch_name, branch in tree.iteritems(
                recursive=True, full_paths=True
            ):
                if len(branch.branches):
                    continue

                try:
                    values = branch.array(
                        entry_start=0,
                        entry_stop=entry_stop,
                        library="ak",
                    )
                    values = ak.to_list(values)
                    if not values:
                        continue

                    safe_values = [
                        json_safe(value, max_vector_entries) for value in values
                    ]
                    branches[branch_name] = (
                        safe_values[0] if len(safe_values) == 1 else safe_values
                    )
                    readable_count += 1
                except Exception as error:
                    unreadable[branch_name] = (
                        f"{type(error).__name__}: {str(error).splitlines()[0]}"
                    )
                    unreadable_count += 1
                    if fingerprint_unreadable:
                        hashes = record_fingerprints(branch, entry_stop)
                        if hashes is None:
                            results["raw_fingerprints"]["unavailable"].setdefault(
                                tree_name, {}
                            )[branch_name] = "Could not split the basket data by entry"
                        else:
                            results["raw_fingerprints"]["trees"].setdefault(
                                tree_name, {}
                            )[branch_name] = hashes

            results["trees"][tree_name] = {
                "entries": entry_stop,
                "branches": branches,
            }
            if unreadable:
                results["unreadable_branches"][tree_name] = unreadable

    with open(output_path, "w", encoding="utf-8") as output_file:
        json.dump(results, output_file, indent=2, allow_nan=False)

    return readable_count, unreadable_count


def main():
    args = parse_args()
    if args.vector_entries is not None and args.vector_entries < 0:
        raise SystemExit("--vector-entries must be zero or greater")

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else Path(f"{input_path}.json")
    readable_count, unreadable_count = dump_root_file(
        input_path,
        output_path,
        args.entries,
        args.vector_entries,
        args.fingerprint_unreadable,
    )
    print(
        f"Wrote {output_path}: {readable_count} readable branches, "
        f"{unreadable_count} unreadable branches skipped"
    )


if __name__ == "__main__":
    main()


