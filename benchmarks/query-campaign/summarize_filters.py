#!/usr/bin/env python3
"""Summarize filter/SBWT ablations and retain quiet-sample acceptance decisions.

Usage: summarize_filters.py LOG_DIRECTORY OUTPUT_JSON
Requires the per-second *-cpu.jsonl records produced by the campaign runner.
No telemetry means no accepted timing, rather than an assumed quiet machine.
Diagnostic rows are untimed and always retained. Minimum accepted cell size: 5.
"""
import collections
import json
import pathlib
import statistics
import sys


def accepted(telemetry, start, end):
    """Require complete coverage and use actual consecutive sample intervals.

    The first record's start is unknown, so it cannot establish coverage. Missing
    leading/trailing telemetry rejects a trial, including a still-running log's
    unmonitored tail. Monitor delays widen intervals rather than leave gaps.
    """
    if (len(telemetry) < 2 or start < telemetry[0]['time']
            or end > telemetry[-1]['time'] or end < start):
        return False
    overlap = [current for previous, current in zip(telemetry, telemetry[1:])
               if current['time'] >= start and previous['time'] <= end]
    return bool(overlap) and all(
        x['outside_busy_cpus'] < 2 and x['sibling']['busy'] <= 0.01 for x in overlap
    )


def summarize(log_directory):
    trials, diagnostics, builds, work, exact = [], [], [], [], []
    records = []
    for log in sorted(log_directory.glob('filter-*.log')):
        telemetry_path = log.with_name(log.stem + '-cpu.jsonl')
        telemetry = [json.loads(s) for s in telemetry_path.read_text().splitlines()] if telemetry_path.exists() else []
        current = None
        for line in log.read_text().splitlines():
            p = line.split(',')
            if p[0] == 'LAYOUT':
                _, round_, layout, workload, method, rep, n, ns, start, end = p
                trials.append(dict(run=log.stem, round=int(round_), layout=layout,
                    workload=workload, method=method, rep=int(rep), queries=int(n),
                    ns=float(ns), start=float(start), end=float(end),
                    quiet=accepted(telemetry, float(start), float(end))))
            elif p[0] in ('OPPORTUNITY', 'ORACLE_SET', 'BASELINE_SEEDS', 'ORACLE_BUILD',
                           'FILTER_BUILD', 'FILTER_PROBES', 'ZOR_REMAINDER', 'EXTRACT', 'HITS'):
                # Preserve named raw diagnostic records; schema is documented
                # in the campaign README. These never count as timed samples.
                records.append(dict(run=log.stem, record=p[0], fields=p[1:]))
            elif p[0] == 'DIAG':
                _, workload, head, found, kind, n, p4, b7, b255, w16 = p
                diagnostics.append(dict(run=log.stem, workload=workload,
                    run_head=head == 'true', found=found == 'true', kind=kind,
                    windows=int(n), packed4_pass=int(p4), bits7_pass=int(b7),
                    byte255_pass=int(b255), wide16_pass=int(w16)))
            elif p[0] == 'BUILD':
                builds.append(dict(run=log.stem, round=int(p[1]), layout=p[2],
                    bytes=int(p[3]), seconds=float(p[4])))
            elif p[0] == 'EXACT' or line.startswith('CHECK ablation'):
                exact.append(dict(run=log.stem, check=line))
            elif p[0] == 'WORK':
                work.append(dict(run=log.stem, workload=p[1], skip=p[2] == 'true',
                    **dict(zip(['chars', 'continuations', 'continuations_found',
                        'restarts', 'table_hits', 'restart_extensions'], map(int, p[3:])))))
            elif p[0] == 'TRIAL':
                # TRIAL precedes the untimed warmup. Accepting that entire interval
                # is conservative; it does not alter the reported timed duration.
                if current is not None and 'ns' in current:
                    current['end'] = float(p[4])
                    current['quiet'] = accepted(telemetry, current['start'], current['end'])
                    trials.append(current)
                current = dict(run=log.stem, rep=int(p[1]), workload=p[2],
                    start=float(p[4]), layout='sbwt-ablation', round=0)
            elif p[0] == 'TIME' and current is not None:
                current.update(method=p[2], queries=int(p[3]), found=int(p[4]), ns=float(p[6]))
        if current is not None and 'ns' in current:
            current['end'] = log.stat().st_mtime
            current['quiet'] = accepted(telemetry, current['start'], current['end'])
            trials.append(current)
    groups = collections.defaultdict(list)
    for row in trials:
        if row['quiet']:
            groups[(row['run'], row['layout'], row['workload'], row['method'])].append(row['ns'])
    summaries = []
    for key, values in sorted(groups.items()):
        summaries.append(dict(zip(['run', 'layout', 'workload', 'method'], key),
            n=len(values), sufficient=len(values) >= 5, median=statistics.median(values),
            minimum=min(values), maximum=max(values)))
    return dict(trials=trials, summaries=summaries, diagnostics=diagnostics,
                builds=builds, work=work, exact=exact, records=records)


if __name__ == '__main__':
    results = summarize(pathlib.Path(sys.argv[1]))
    pathlib.Path(sys.argv[2]).write_text(json.dumps(results, indent=2) + '\n')
    print(f"{len(results['trials'])} trials; "
          f"{sum(x['quiet'] for x in results['trials'])} accepted; "
          f"{len(results['summaries'])} summary cells")
