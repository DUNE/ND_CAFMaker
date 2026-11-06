import argparse
import html
import json
import os
from pathlib import Path
from urllib.parse import urlparse


MISSING = object()
STAGES = {
    "resolve-inputs": "Find prepared CAFMaker inputs",
    "prepare-pr-inputs": "Prepare inputs for PR branch",
    "build": "Build (prof)",
    "run-cafmaker": "Run CAFMaker",
    "verify": "Verify CAF output",
}


def changes(reference, candidate, path=""):
    if isinstance(reference, dict) and isinstance(candidate, dict):
        for key in sorted(reference.keys() | candidate.keys()):
            pointer = key.replace("~", "~0").replace("/", "~1")
            yield from changes(
                reference.get(key, MISSING), candidate.get(key, MISSING),
                f"{path}/{pointer}",
            )
    elif isinstance(reference, list) and isinstance(candidate, list):
        for position in range(max(len(reference), len(candidate))):
            yield from changes(
                reference[position] if position < len(reference) else MISSING,
                candidate[position] if position < len(candidate) else MISSING,
                f"{path}/{position}",
            )
    elif type(reference) is not type(candidate) or reference != candidate:
        kind = "Added" if reference is MISSING else "Removed" if candidate is MISSING else "Changed"
        yield path or "/", kind, reference, candidate


def cell(value, limit=160):
    text = str(value)
    if len(text) > limit:
        text = text[:limit] + "..."
    return "<code>" + html.escape(text).replace("|", "&#124;").replace("\n", "&#10;") + "</code>"


def value_cell(value):
    if value is MISSING:
        return cell("(missing)")
    return cell(json.dumps(value, ensure_ascii=True, sort_keys=True))


def raw_diff(diff):
    lines = diff.splitlines()
    excerpt = "\n".join(lines[:200])[:40000]
    truncated = len(lines) > 200 or len("\n".join(lines[:200])) > 40000
    note = "\n\nPreview truncated; the diagnostic artifact contains the full diff." if truncated else ""
    return (
        "<details><summary>Raw unified diff</summary>\n\n<pre>"
        + html.escape(excerpt) + "</pre>" + note + "\n\n</details>\n"
    )


def render_regression(reference, candidate, diff):
    counts = {"Added": 0, "Removed": 0, "Changed": 0}
    rows = []
    for path, kind, before, after in changes(reference, candidate):
        counts[kind] += 1
        if len(rows) < 50:
            rows.append(f"| {cell(path)} | {kind} | {value_cell(before)} | {value_cell(after)} |")
    total = sum(counts.values())
    result = [
        "### Regression Differences", "",
        f"**{total} changed fields**: {counts['Changed']} changed, "
        f"{counts['Added']} added, {counts['Removed']} removed.", "",
    ]
    if rows:
        result += [
            "Fields use JSON Pointer paths; array indices are zero-based.", "",
            "| Field | Change | Reference | Generated |",
            "| :--- | :--- | :--- | :--- |", *rows, "",
        ]
        if total > len(rows):
            result += ["Showing the first 50 fields. Full values and differences are in the diagnostic artifact.", ""]
    else:
        result += ["JSON values are identical; the exact comparison detected formatting differences.", ""]
    result.append(raw_diff(diff))
    return "\n".join(result)


def link(url, label):
    if urlparse(url).scheme != "https":
        return label
    return f'<a href="{html.escape(url, quote=True)}">{html.escape(label)}</a>'


def render_summary(needs, jobs, regression="", preparation_only=False, run_url="", artifact_url="", reporting_errors=()):
    active = ["prepare-pr-inputs"] if preparation_only else ["resolve-inputs", "build", "run-cafmaker", "verify"]
    comparison = needs.get("verify", {}).get("outputs", {}).get("comparison", "")
    mismatch = not preparation_only and comparison == "mismatch"
    by_name = {job["name"]: job for job in jobs}
    unexpected = []
    failed_steps = []
    for stage in active:
        result = needs.get(stage, {}).get("result", "unknown")
        if result != "success" and not (stage == "verify" and mismatch and result == "failure"):
            unexpected.append(stage)
        job = by_name.get(STAGES[stage], {})
        for step in job.get("steps", []):
            if step.get("conclusion") not in ("failure", "cancelled"):
                continue
            if stage == "verify" and mismatch and step["name"] == "Dump and compare CAF output":
                continue
            failed_steps.append((STAGES[stage], step["name"], job.get("html_url", run_url)))
    if unexpected or failed_steps or reporting_errors:
        title = "CI Failure"
        description = "One or more required stages failed, were cancelled, or did not run."
    elif mismatch:
        title = "Regression Mismatch"
        description = "CAFMaker completed, but the generated output differs from the committed reference."
    elif preparation_only:
        title = "Input Preparation Passed"
        description = "The input artifact is ready. This branch-only dispatch did not build or run the regression test."
    elif comparison == "match":
        title = "Passed"
        description = "All required stages passed. The generated CAF matches the committed reference exactly."
    else:
        title = "CI Failure"
        description = "No confirmed regression comparison result was reported."

    report = ["# ND_CAFMaker CI", "", f"## {title}", "", description, "",
              "| Stage | Result | Logs |", "| :--- | :--- | :--- |"]
    for stage in active:
        result = needs.get(stage, {}).get("result", "unknown")
        label = "Mismatch" if stage == "verify" and mismatch and result == "failure" else result.capitalize()
        job = by_name.get(STAGES[stage], {})
        report.append(f"| {STAGES[stage]} | **{label}** | {link(job.get('html_url', run_url), 'View logs')} |")
    report += [""]
    if not jobs:
        report += ["Job-step details are unavailable; use the run logs to diagnose failures.", ""]
    if failed_steps:
        report += ["### Failed Steps", ""]
        for stage_name, step_name, url in failed_steps:
            report.append(f"- {cell(stage_name)}: {cell(step_name)} ({link(url, 'logs')})")
        report += [""]
    if unexpected:
        report += ["### Other Failures", ""]
        for stage in unexpected:
            result = needs.get(stage, {}).get("result", "unknown")
            report.append(f"- **{STAGES[stage]}**: {result}. See the linked job logs.")
        report += [""]
    if reporting_errors:
        report += ["### Reporting Failures", ""]
        report += [f"- {cell(error)}" for error in reporting_errors]
        report += [""]
    if not preparation_only and not mismatch and comparison != "match":
        explanation = "Verification failed before a usable comparison was completed." if comparison == "error" else "Regression was not completed; check the upstream stages and verification logs."
        report += [explanation, ""]
    if mismatch:
        report += [regression or "Detailed regression preview is unavailable. See the comparison diagnostics and verification logs.", "",
                   "### Intentional Reference Changes", "",
                   "Only if these changes are intentional and justified, run from the repository root:", "",
                   "```bash", "./scripts/caf-regression --update-reference", "```", "",
                   "Review and commit the resulting reference diff. Do not update the reference to hide an unrelated failure.", ""]
    if artifact_url:
        report += [link(artifact_url, "Download comparison diagnostics (generated JSON and full diff)"), ""]
    if run_url:
        report += [link(run_url, "Open workflow run"), ""]
    return "\n".join(report)


def main():
    parser = argparse.ArgumentParser(description="Render ND_CAFMaker regression and workflow summaries.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    regression_parser = subparsers.add_parser("regression")
    for argument in ("reference", "candidate", "diff", "output"):
        regression_parser.add_argument(f"--{argument}", required=True, type=Path)
    summary_parser = subparsers.add_parser("summary")
    summary_parser.add_argument("--jobs", required=True, type=Path)
    summary_parser.add_argument("--regression", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "regression":
        report = render_regression(
            json.loads(args.reference.read_text()), json.loads(args.candidate.read_text()),
            args.diff.read_text(),
        )
        args.output.write_text(report, encoding="utf-8")
    else:
        needs = json.loads(os.environ["CI_NEEDS"])
        regression = args.regression.read_text() if args.regression.is_file() else ""
        report = render_summary(
            needs, json.loads(args.jobs.read_text()), regression,
            os.environ.get("CI_PREPARATION_ONLY") == "true",
            os.environ["CI_RUN_URL"], os.environ.get("CI_ARTIFACT_URL", ""),
            ["Comparison diagnostic artifact could not be downloaded."]
            if os.environ.get("CI_DOWNLOAD_OUTCOME") == "failure" else [],
        )
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as summary:
            summary.write(report)


if __name__ == "__main__":
    main()