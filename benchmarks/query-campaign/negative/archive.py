#!/usr/bin/env python3
"""Export stage-scoped summaries and a compact raw log/telemetry archive.

Run from the workspace root after every controller has exited. The archive maps
child log names to a single common telemetry stream per stage; unpacking must
restore those per-child *-cpu.jsonl files before using summarize_filters.py.
Large reference/key/query files and binaries remain in the local campaign area.
"""
import gzip
import hashlib
import importlib.util
import json
import pathlib
import shutil
import sys

sys.dont_write_bytecode = True
base = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('summaries', base / 'summarize_filters.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
logs = pathlib.Path('campaign/logs')
results = module.summarize(logs)
out = base / 'results/negative'
out.mkdir(parents=True, exist_ok=True)
raw = {'files': {}, 'telemetry_for_log': {}}
for stage, prefixes in [(1, ('filter-negative-pipeline-',)),
                        (2, ('filter-qmer-oracle-',)),
                        (3, ('filter-real-gencode-q16-', 'filter-real-gencode-q18-')),
                        (4, ('filter-real-gencode-q0-', 'filter-real-human-', 'filter-confirm-')),
                        (5, ('filter-net-',)),
                        (6, ('filter-direct-',))]:
    def belongs(row):
        return row['run'].startswith(prefixes)
    subset = {key: [row for row in values if belongs(row)] for key, values in results.items()}
    groups = {(r['run'], r['layout'], r['workload'], r['method']) for r in subset['trials']}
    accepted = {(r['run'], r['layout'], r['workload'], r['method']) for r in subset['summaries'] if r['sufficient']}
    subset['coverage'] = dict(trials=len(subset['trials']), accepted=sum(r['quiet'] for r in subset['trials']),
                              cells=len(groups), insufficient_cells=[list(x) for x in sorted(groups-accepted)])
    (out / f'stage{stage}-results.json').write_text(json.dumps(subset, indent=2)+'\n')
    telemetry = f'filter-negative-stage{stage}-cpu.jsonl'
    for suffix in ['cpu.jsonl', 'gate.jsonl', '.log']:
        name = f'filter-negative-stage{stage}' + (suffix if suffix.startswith('.') else '-'+suffix)
        path = logs/name
        if path.exists(): raw['files'][name] = path.read_text()
    for path in logs.glob('*.log'):
        if path.stem.startswith(prefixes):
            raw['files'][path.name] = path.read_text()
            raw['telemetry_for_log'][path.name] = telemetry
    provenance = pathlib.Path(f'campaign/negative-query/stage{stage}-provenance.json')
    if provenance.exists(): shutil.copy2(provenance, out/provenance.name)
    print(stage, subset['coverage'])
for path in logs.glob('extract-*.log'):
    if 'gencode-q' in path.name or 'human-q' in path.name: raw['files'][path.name] = path.read_text()
(out/'raw.json.gz').write_bytes(gzip.compress(json.dumps(raw).encode(), mtime=0))
shutil.copy2('campaign/negative-query/trace-validation.json', out/'trace-validation.json')
shutil.copy2('campaign/negative-query/validation.json', out/'validation.json')
inputs = pathlib.Path('campaign/negative-query/inputs.json')
if inputs.exists(): shutil.copy2(inputs, out/'inputs.json')
shutil.copy2('campaign/negative-query/plan.json', out/'campaign.json')
manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.name != 'SHA256.json'}
(out/'SHA256.json').write_text(json.dumps(manifest, indent=2)+'\n')
