#!/usr/bin/env python3
"""Gate a pinned filter benchmark and pause between trials during host contention.

Usage: run_filter_quiet.py --logs DIR [--cpu 24] LABEL COMMAND [ARGS...]
Never signals or changes other processes. The Rust harness cooperatively waits
on QUERY_QUIET_FILE outside each warmup/timing interval. Trials overlapping a new
busy sample are still excluded by summarize_filters.py, rather than repaired.
"""
import argparse
import json
import os
import pathlib
import subprocess
import time


def snapshot():
    return {p[0]: list(map(int, p[1:9]))
            for line in pathlib.Path('/proc/stat').read_text().splitlines()
            if (p := line.split()) and p[0].startswith('cpu')}


def usage(before, after):
    result = {}
    for key, old in before.items():
        delta = [b-a for a, b in zip(old, after[key])]
        total = max(1, sum(delta))
        result[key] = dict(busy=(total-delta[3]-delta[4])/total,
                           iowait=delta[4]/total, steal=delta[7]/total)
    return result


def write_status(path, ready):
    temp = path.with_suffix('.tmp')
    temp.write_text('ready' if ready else 'wait')
    temp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--logs', required=True, type=pathlib.Path)
    parser.add_argument('--cpu', type=int, default=24)
    parser.add_argument('label')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.command:
        parser.error('a benchmark command is required')
    siblings = []
    topology = pathlib.Path(f'/sys/devices/system/cpu/cpu{args.cpu}/topology/thread_siblings_list')
    for item in topology.read_text().strip().split(','):
        lo, _, hi = item.partition('-')
        siblings.extend(range(int(lo), int(hi or lo)+1))
    other = [x for x in siblings if x != args.cpu]
    args.logs.mkdir(parents=True, exist_ok=True)
    prefix = args.logs / args.label
    before = snapshot()
    streak = 0
    with prefix.with_name(prefix.name+'-gate.jsonl').open('w', buffering=1) as gate:
        while streak < 3:
            time.sleep(10)
            now = snapshot()
            current = usage(before, now)
            before = now
            quiet = (current['cpu']['busy'] < .04 and current['cpu']['iowait'] < .01
                     and all(current[f'cpu{x}']['busy'] < .05
                             and current[f'cpu{x}']['steal'] < .01 for x in siblings))
            streak = streak+1 if quiet else 0
            gate.write(json.dumps(dict(time=time.time(), host=current['cpu'],
                siblings={x: current[f'cpu{x}'] for x in siblings},
                quiet=quiet, streak=streak))+'\n')
    status = prefix.with_name(prefix.name+'-status').resolve()
    write_status(status, True)
    environment = dict(os.environ, QUERY_QUIET_FILE=str(status))
    with prefix.with_name(prefix.name+'.log').open('w') as log, prefix.with_name(prefix.name+'-cpu.jsonl').open('w', buffering=1) as telemetry:
        before = snapshot()
        child = subprocess.Popen(['taskset', '-c', str(args.cpu), *args.command],
            stdout=log, stderr=subprocess.STDOUT, env=environment)
        streak = 3
        while child.poll() is None:
            time.sleep(1)
            now = snapshot()
            current = usage(before, now)
            before = now
            outside = sum(x['busy'] for k, x in current.items()
                          if k not in ('cpu', f'cpu{args.cpu}'))
            sibling = max((current[f'cpu{x}'] for x in other), key=lambda x: x['busy'],
                          default=dict(busy=0, iowait=0, steal=0))
            quiet = outside < 2 and sibling['busy'] <= .01
            streak = streak+1 if quiet else 0
            write_status(status, streak >= 3)
            telemetry.write(json.dumps(dict(time=time.time(), outside_busy_cpus=outside,
                sibling=sibling, host=current['cpu'], ready=streak >= 3))+'\n')
        if child.returncode:
            raise SystemExit(child.returncode)


if __name__ == '__main__':
    main()
