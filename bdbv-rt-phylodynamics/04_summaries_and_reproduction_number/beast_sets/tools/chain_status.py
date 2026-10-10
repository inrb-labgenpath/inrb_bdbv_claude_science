#!/usr/bin/env python3
"""chain_status.py - reads the state of every chain of the queue from its directory (reads only; changes nothing).
usage: python tools/chain_status.py <root> [--append-progress]
prints one line per chain that was started: label, state (running/finished/interrupted), exit code, last state logged, newest saved state.
--append-progress adds one row per running chain to sched/progress.csv (utc, chain, last_state_logged, mtime of the log)."""
import csv, fcntl, glob, os, re, sys, time, datetime as dt

STEPS = int(os.environ.get('STEPS', 40_000_000))


def lock_held(path):
    if not os.path.exists(path):
        return False
    with open(path, 'a') as fh:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return True
        fcntl.flock(fh, fcntl.LOCK_UN)
        return False

def last_state(log):
    if not os.path.exists(log):
        return None
    with open(log, 'rb') as fh:
        fh.seek(0, 2); n = fh.tell(); fh.seek(max(0, n - 4000))
        lines = fh.read().decode(errors='replace').splitlines()
    for l in reversed(lines):
        m = re.match(r'^(\d+)\t', l)
        if m and l.count('\t') > 5:
            return int(m.group(1))
    return None

def saved_state(statefile):
    if not os.path.exists(statefile):
        return None
    with open(statefile, errors='replace') as fh:
        for i, l in enumerate(fh):
            if l.startswith('state\t'):
                return int(l.split('\t')[1])
            if i > 5:
                break
    return None

def read_queue(root):
    with open(os.path.join(root, 'sched', 'queue.tsv')) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

def status(root):
    out = []
    for q in read_queue(root):
        d = os.path.join(root, 'runs', q['analysis'], q['label'])
        if not os.path.isdir(d):
            out.append(dict(q, state='waiting')); continue
        log = os.path.join(d, q['xml_stem'] + '.log')
        ec = None
        p = os.path.join(d, 'exit_code.txt')
        if os.path.exists(p) and os.path.getsize(p):
            ec = open(p).read().strip()
        # the lock is probed only for a chain that was started and has no exit code (a wrapper takes its lock before it
        # writes started.txt, so that a probe cannot meet the first start of a chain)
        if ec is not None and ec not in ('129', '130', '137', '143'):
            st = 'finished'
        elif ec is None and last_state(log) == STEPS and os.path.exists(os.path.join(d, 'started.txt')) and not lock_held(os.path.join(d, 'chain.lock')):
            st = 'finished'; ec = 'not recorded'          # judged by its files
        elif not os.path.exists(os.path.join(d, 'started.txt')):
            st = 'waiting'
        else:
            st = 'running' if lock_held(os.path.join(d, 'chain.lock')) else 'interrupted'
        out.append(dict(q, state=st, exit_code=ec, last_state_logged=last_state(log), newest_saved_state=saved_state(os.path.join(d, q['label'] + '.state')),
                        log_mtime=(os.path.getmtime(log) if os.path.exists(log) else None), dir=d,
                        started=(open(os.path.join(d, 'started.txt')).read().strip() if os.path.exists(os.path.join(d, 'started.txt')) else '')))
    return out

if __name__ == '__main__':
    root = os.path.abspath(sys.argv[1])
    S = status(root)
    now = dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    n = {k: sum(1 for s in S if s['state'] == k) for k in ('running', 'finished', 'interrupted', 'waiting')}
    print(now, n)
    for s in S:
        if s['state'] != 'waiting':
            print(f"{s['label']:22s} {s['state']:11s} exit {str(s.get('exit_code')):5s} last state {str(s.get('last_state_logged')):>9s} saved state {str(s.get('newest_saved_state')):>9s} started {s.get('started')}")
    if '--append-progress' in sys.argv:
        p = os.path.join(root, 'sched', 'progress.csv')
        new = not os.path.exists(p)
        with open(p, 'a') as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(['utc', 'chain', 'analysis', 'state_of_chain', 'last_state_logged', 'log_mtime_utc', 'newest_saved_state'])
            for s in S:
                if s['state'] in ('running', 'finished', 'interrupted'):
                    mt = dt.datetime.fromtimestamp(s['log_mtime'], dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ') if s.get('log_mtime') else ''
                    w.writerow([now, s['label'], s['analysis'], s['state'], s.get('last_state_logged'), mt, s.get('newest_saved_state')])
