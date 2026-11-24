#!/usr/bin/env python3
"""Write a tidymarkers trace from a `marker_trace_payload()` file.

This is the second implementation of the record `marker_artefact()` writes: it
reads the payload R emits, checks it against `inst/extdata/trace-schema.json`,
and re-renders `artefact.json`, `decisions.csv`, `checks.csv` and `why.csv`
byte for byte.  It deliberately does not re-derive anything: the manifest comes
from the payload (where R's atomic-vs-list distinction survives) and the table
rows come from the payload already read with the schema's types, so a schema
drift fails here instead of passing silently.

Usage:
    python3 marker_trace.py PAYLOAD.json OUTDIR [--solved solved.rds] [--schema PATH]

Standard library only.  `solved.rds` is an R serialisation that Python cannot
author, so `--solved` copies the file verbatim; the conformance test then
replays both directories through `marker_replay()`, which re-solves from that
file, so the copy is not the part that is being trusted.
"""

import argparse
import csv
import json
import math
import shutil
import sys
from pathlib import Path

PAYLOAD_SCHEMA = "tidymarkers/trace-payload@1"
COLUMN_TYPES = ("string", "integer", "number", "logical")


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def default_schema_path():
    return Path(__file__).resolve().parent.parent / "extdata" / "trace-schema.json"


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------
def payload_problems(payload, schema):
    """Every divergence between the payload and the schema, as one line each."""
    problems = []
    if payload.get("schema") != PAYLOAD_SCHEMA:
        problems.append(
            "payload schema is %r, expected %r" % (payload.get("schema"), PAYLOAD_SCHEMA)
        )
    if payload.get("trace_schema") != schema.get("trace_schema"):
        problems.append(
            "payload trace_schema is %r, schema says %r"
            % (payload.get("trace_schema"), schema.get("trace_schema"))
        )

    known = {table["name"]: table for table in schema["tables"]}
    for name, table in payload.get("tables", {}).items():
        if name not in known:
            problems.append("table %r is not in the schema" % name)
            continue
        if table.get("file") != known[name]["file"]:
            problems.append(
                "%s: file is %r, schema says %r"
                % (name, table.get("file"), known[name]["file"])
            )
        want = [(col["name"], col["type"]) for col in known[name]["columns"]]
        got = [(col.get("name"), col.get("type")) for col in table.get("columns", [])]
        if got != want:
            problems.append(
                "%s: columns diverge; schema %s, payload %s" % (name, want, got)
            )
        unknown = sorted({kind for _, kind in got if kind not in COLUMN_TYPES})
        if unknown:
            problems.append("%s: unknown column type(s) %s" % (name, unknown))
        for index, row in enumerate(table.get("rows", [])):
            if len(row) != len(want):
                problems.append(
                    "%s: row %d has %d cells, schema has %d"
                    % (name, index + 1, len(row), len(want))
                )
                break

    missing = [
        key for key in schema["manifest"]["required"] if key not in payload.get("manifest", {})
    ]
    if missing:
        problems.append("manifest is missing required key(s) %s" % missing)
    inline = payload.get("inline_arrays", [])
    if not isinstance(inline, list):
        problems.append("inline_arrays must be a list of dot paths")
    return problems


# --------------------------------------------------------------------------
# CSV cell formatting
# --------------------------------------------------------------------------
def _digits_and_power(value):
    """Significant digits and the power of ten of the first one.

    Python's `repr(float)` is the shortest string that round-trips, which is
    the same set of digits readr uses; only the fixed-versus-exponent choice
    differs, so only that choice is re-implemented below.
    """
    text = repr(abs(value))
    if "e" in text:
        mantissa, _, exponent = text.partition("e")
        exponent = int(exponent)
    else:
        mantissa, exponent = text, 0
    whole, _, fraction = mantissa.partition(".")
    combined = whole + fraction
    first = len(combined) - len(combined.lstrip("0"))
    digits = combined[first:].rstrip("0")
    power = exponent + len(whole) - first - 1
    return digits, power


def format_number(value):
    """A double exactly as readr::write_csv() spells it."""
    if isinstance(value, bool):
        raise TypeError("logical in a number column")
    if isinstance(value, int):
        return str(value)
    if math.isnan(value):
        return "NA"
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    if value == 0.0:
        return "-0" if math.copysign(1.0, value) < 0 else "0"
    digits, power = _digits_and_power(value)
    sign = "-" if value < 0 else ""
    count = len(digits)

    def fixed():
        if power >= 0:
            if count <= power + 1:
                return digits + "0" * (power + 1 - count)
            return digits[: power + 1] + "." + digits[power + 1 :]
        return "0." + "0" * (-power - 1) + digits

    def scientific():
        if power < 0:
            mantissa = digits[0] + ("." + digits[1:] if count > 1 else "")
            return mantissa + "e" + str(power)
        return digits + "e" + str(power - (count - 1))

    if power < -3:
        return sign + scientific()
    if power >= 15 and len(scientific()) < len(fixed()):
        return sign + scientific()
    return sign + fixed()


def cell_text(value, kind):
    if value is None:
        return "NA"
    if kind == "string":
        return str(value)
    if kind == "logical":
        return "TRUE" if value else "FALSE"
    if kind == "integer":
        return str(int(value))
    if kind == "number":
        return format_number(value)
    raise ValueError("unknown column type %r" % kind)


# --------------------------------------------------------------------------
# Manifest rendering
# --------------------------------------------------------------------------
def _scalar_text(value):
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        # jsonlite digits = NA is C's %.15g, including its zero-padded exponent.
        return "%.15g" % value
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    raise TypeError("cannot render %r as JSON" % (value,))


def render_manifest(value, inline, path="", indent=0):
    """The layout jsonlite::toJSON(pretty = TRUE) writes.

    A dict and a list are multiline unless the list's dot path is named in
    `inline`; an atomic vector is inline because R hands it to jsonlite as one
    value, and JSON cannot carry that distinction, which is why R names it.
    """
    pad = "  " * indent
    if isinstance(value, dict):
        if not value:
            return "{}"
        lines = []
        for key, item in value.items():
            child = path + "." + key if path else key
            lines.append(
                '%s  %s: %s'
                % (pad, json.dumps(key, ensure_ascii=False), render_manifest(item, inline, child, indent + 1))
            )
        return "{\n" + ",\n".join(lines) + "\n" + pad + "}"
    if isinstance(value, list):
        if not value:
            return "[]"
        if path in inline:
            return "[" + ", ".join(_scalar_text(item) for item in value) + "]"
        lines = [pad + "  " + render_manifest(item, inline, path, indent + 1) for item in value]
        return "[\n" + ",\n".join(lines) + "\n" + pad + "]"
    return _scalar_text(value)


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------
def write_manifest(payload, schema, outdir):
    text = render_manifest(payload["manifest"], set(payload.get("inline_arrays", [])))
    path = outdir / schema["manifest"]["file"]
    path.write_text(text + "\n", encoding="utf-8")


def write_table(table, outdir):
    path = outdir / table["file"]
    columns = table["columns"]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        writer.writerow([col["name"] for col in columns])
        for row in table["rows"]:
            writer.writerow(
                [cell_text(cell[0] if isinstance(cell, list) else cell, col["type"])
                 for cell, col in zip(row, columns)]
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("payload", type=Path, help="JSON from marker_trace_payload()")
    parser.add_argument("outdir", type=Path, help="directory to write the four files into")
    parser.add_argument("--solved", type=Path, help="solved.rds to copy verbatim")
    parser.add_argument("--schema", type=Path, default=None, help="conformance schema")
    args = parser.parse_args(argv)

    schema_path = args.schema or default_schema_path()
    if not schema_path.is_file():
        print("no conformance schema at %s" % schema_path, file=sys.stderr)
        return 2
    schema = read_json(schema_path)
    payload = read_json(args.payload)

    problems = payload_problems(payload, schema)
    if problems:
        for problem in problems:
            print(problem, file=sys.stderr)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_manifest(payload, schema, args.outdir)
    for table in payload["tables"].values():
        write_table(table, args.outdir)
    if args.solved:
        shutil.copyfile(args.solved, args.outdir / "solved.rds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
