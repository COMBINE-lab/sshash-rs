#!/usr/bin/env python3
"""Fingerprint local campaign inputs after timing, without copying large data.

Run from the original workspace root. This records file contents, sizes and
resolved paths, including full reference key sets and untimed mapper traces.
Do not run disk/hash scans concurrently with performance measurements.
"""
import hashlib
import json
import pathlib

paths = [pathlib.Path('campaign/index/gencode'+suffix) for suffix in ('.ssi','.ssi.mphf')]
paths += [pathlib.Path('/scratch1/alevin-fry-ecosystem/tmp/sux-phastr/idx/human_r0152_t32'+suffix) for suffix in ('.ssi','.ssi.mphf')]
paths += [pathlib.Path('campaign/reads/rna100k_1.fq')]
for pattern in ('*.keys','*.trace','*.fq'):
    paths += sorted(pathlib.Path('campaign/negative-query').glob(pattern))
records=[]
for path in paths:
    digest=hashlib.sha256()
    with path.open('rb') as source:
        while block:=source.read(8*1024*1024):digest.update(block)
    records.append(dict(path=str(path),resolved_path=str(path.resolve()),bytes=path.stat().st_size,sha256=digest.hexdigest()))
    print(path,flush=True)
pathlib.Path('campaign/negative-query/inputs.json').write_text(json.dumps(records,indent=2)+'\n')
