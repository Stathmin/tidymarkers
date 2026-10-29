#!/usr/bin/env python3
"""Submit the R dependency tree to GitHub's dependency graph.

GitHub cannot read `DESCRIPTION` or `renv.lock`, so the graph stays empty
("No dependencies found") until a snapshot is pushed through the
dependency-submission API. Versions and the dependency edges come from
`renv.lock` (the pinned R 4.6.1 / Bioconductor 3.23 snapshot); which packages
are *direct* comes from `DESCRIPTION`.

Usage:
    python3 .github/scripts/dependency-snapshot.py --dry-run   # print the JSON
    python3 .github/scripts/dependency-snapshot.py             # POST it

Exit status is non-zero when the API rejects the snapshot, so a broken payload
fails the workflow instead of silently leaving the graph empty.
"""

import argparse
import datetime
import json
import os
import subprocess
import sys
import urllib.error
import urllib.request

#: Packages shipped with R. They have no version of their own and no purl.
BASE_PACKAGES = {
    "base", "compiler", "datasets", "grDevices", "graphics", "grid", "methods",
    "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils",
}

#: Fields of a lockfile entry that carry hard dependencies.
DEP_FIELDS = ("Depends", "Imports", "LinkingTo")


def read_dcf(path):
    """Parse a Debian-control-style file (`DESCRIPTION`)."""
    fields, key = {}, None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                key = None
            elif line[0] in " \t":
                if key:
                    fields[key] += " " + line.strip()
            elif ":" in line:
                key, _, value = line.partition(":")
                key = key.strip()
                fields[key] = value.strip()
    return fields


def dep_names(value):
    """`"dplyr, R (>= 4.1.0), Biostrings"` -> `["dplyr", "Biostrings"]`."""
    items = value if isinstance(value, list) else str(value or "").split(",")
    names = []
    for item in items:
        name = str(item).split("(")[0].strip()
        if name and name != "R":
            names.append(name)
    return names


def package_url(name, version, repository):
    """A purl GitHub accepts.

    GitHub's dependency-submission API rejects `pkg:bioconductor/...`
    ("invalid package url type"), so every package is reported as `cran`.
    The name and version are exact; only the ecosystem label is coarser.
    """
    return f"pkg:cran/{name}@{version}"


def build_snapshot(sha, ref, run_id):
    lock = json.load(open("renv.lock", encoding="utf-8"))
    packages = lock["Packages"]
    description = read_dcf("DESCRIPTION")

    runtime = set()
    for field in ("Depends", "Imports", "LinkingTo"):
        runtime.update(dep_names(description.get(field)))
    development = set(dep_names(description.get("Suggests"))) - runtime

    url_of = {
        name: package_url(name, meta["Version"], meta.get("Repository"))
        for name, meta in packages.items()
        if name not in BASE_PACKAGES
    }

    resolved = {}
    for name, meta in packages.items():
        if name in BASE_PACKAGES:
            continue
        dependencies = []
        for field in DEP_FIELDS:
            for dependency in dep_names(meta.get(field)):
                if dependency in url_of:
                    dependencies.append(url_of[dependency])
        resolved[url_of[name]] = {
            "package_url": url_of[name],
            "relationship": "direct" if name in runtime or name in development else "indirect",
            "scope": "development" if name in development else "runtime",
            "dependencies": sorted(set(dependencies)),
        }

    scanned = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    return {
        "version": 0,
        "sha": sha,
        "ref": ref,
        "scanned": scanned,
        "job": {"correlator": "tidymarkers-dependency-submission", "id": run_id},
        "detector": {
            "name": "tidymarkers-dependency-snapshot",
            "version": "1.0.0",
            "url": "https://github.com/Stathmin/tidymarkers",
        },
        "manifests": {
            "renv.lock": {
                "name": "tidymarkers",
                "file": {"source_location": "renv.lock"},
                "resolved": resolved,
            }
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="print, do not POST")
    parser.add_argument("--sha", default=os.environ.get("GITHUB_SHA"))
    parser.add_argument("--ref", default=os.environ.get("GITHUB_REF"))
    parser.add_argument("--run-id", default=os.environ.get("GITHUB_RUN_ID", "local"))
    args = parser.parse_args()

    sha = args.sha or subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True
    ).strip()
    ref = args.ref or "refs/heads/main"

    snapshot = build_snapshot(sha, ref, args.run_id)
    body = json.dumps(snapshot)
    if args.dry_run:
        print(body)
        return 0

    repository = os.environ["GITHUB_REPOSITORY"]
    token = os.environ["GITHUB_TOKEN"]
    request = urllib.request.Request(
        f"https://api.github.com/repos/{repository}/dependency-graph/snapshots",
        data=body.encode("utf-8"),
        method="POST",
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2022-11-28",
            "Content-Type": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(request) as response:
            print(f"submitted {len(snapshot['manifests']['renv.lock']['resolved'])} packages "
                  f"({response.status})")
    except urllib.error.HTTPError as error:
        print(error.read().decode("utf-8", "replace"), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
