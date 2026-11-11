#!/usr/bin/env python3
"""Restore archived logs and child telemetry for independent re-summarization.

Usage: unpack.py RAW_JSON_GZ EMPTY_OUTPUT_DIRECTORY
Only basename entries from our archive format are accepted. Existing files are
not overwritten; choose an empty directory so stale logs cannot enter summaries.
"""
import gzip
import json
import pathlib
import sys

archive = json.loads(gzip.decompress(pathlib.Path(sys.argv[1]).read_bytes()))
out = pathlib.Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
if any(out.iterdir()):
    raise SystemExit('output directory must be empty')
for name, contents in archive['files'].items():
    if pathlib.Path(name).name != name or name in ('.', '..'):
        raise SystemExit('archive entry must be a basename')
    (out/name).write_text(contents)
for log, telemetry in archive['telemetry_for_log'].items():
    if pathlib.Path(log).name != log or pathlib.Path(telemetry).name != telemetry:
        raise SystemExit('archive mapping must use basenames')
    target = out/(pathlib.Path(log).stem+'-cpu.jsonl')
    target.write_text((out/telemetry).read_text())
