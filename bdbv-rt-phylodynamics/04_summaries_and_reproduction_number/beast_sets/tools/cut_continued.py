#!/usr/bin/env python3
"""cut_continued.py - input of the summary CUT AT THE SAVED STATES.
Reads the chain directories; changes nothing in them.
usage: python tools/cut_continued.py <root> <analysis> <outdir> [--steps 40000000 --log-every 20000] [--runs runs]
For every FINISHED chain of the analysis it writes <outdir>/<analysis>/<label>/<stem>.log:
  chain without a continuation: a copy of the log of the chain directory, unchanged (SHA-256 compared);
  continued chain: the log of the FIRST part (the chain as it ran before it was ended from outside), cut at the state N of the saved state
      from which its first continuation was started: the head and every row with a state of at most N, that state included. Nothing that
      a continuation logged is used. No line is written by this tool itself.
  No tree file is written: the estimates are derived from the logs.
The script STOPS (exit code 1, no file for the analysis is left) unless the cut log of a continued chain holds every state once, in steps
of --log-every from 0 to N, and N is the state of the file loaded_state_<N>.state of the first continuation (read from the file itself).
<outdir>/<analysis>/cut_report.csv   one row per chain
prints one line per chain."""
import argparse, csv, hashlib, os, re, shutil, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chain_status as CS


def read_log(path):
    """head (comment lines and the line of the column names), rows as (state, line): complete lines only"""
    head, cols, rows, cut = [], None, [], 0
    with open(path, newline='') as fh:
        for line in fh:
            if line.startswith('#'):
                head.append(line); continue
            if cols is None:
                cols = line.rstrip('\n').split('\t'); head.append(line); continue
            if not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != len(cols) or not f[0].isdigit() or not line.endswith('\n'):
                cut += 1; continue                         # a last line that was cut when the chain was ended
            rows.append((int(f[0]), line))
    return head, cols, rows, cut


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root'); ap.add_argument('analysis'); ap.add_argument('outdir')
    ap.add_argument('--steps', type=int, default=40_000_000); ap.add_argument('--log-every', type=int, default=20_000); ap.add_argument('--runs', default='runs')
    ap.add_argument('--burnin', type=float, default=0.30, help='share of the logged states that the stored function discards (used for the columns of the report only; the stored function itself is applied by the derivation)')
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    out = os.path.join(os.path.abspath(a.outdir), a.analysis); tmp = out + '.tmp'
    shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
    q = [r for r in CS.read_queue(root) if r['analysis'] == a.analysis]
    rep = []
    try:
        for r in q:
            d = os.path.join(root, a.runs, a.analysis, r['label']); stem = r['xml_stem']
            if not os.path.isdir(d):
                continue
            parts = [d]; k = 1
            while os.path.isdir(os.path.join(d, f'continuation_{k}')):
                parts.append(os.path.join(d, f'continuation_{k}')); k += 1
            last = parts[-1]
            ecf = os.path.join(last, 'exit_code.txt')
            if not (os.path.exists(ecf) and open(ecf).read().strip() not in ('', '129', '130', '137', '143')):
                continue                                   # not finished
            t = os.path.join(tmp, r['label']); os.makedirs(t)
            src = os.path.join(d, stem + '.log'); dst = os.path.join(t, stem + '.log')
            if len(parts) == 1:
                shutil.copy2(src, dst); assert sha(src) == sha(dst)
                head, cols, rows, cut = read_log(dst)
                n = len(rows); nb = int(n * a.burnin)
                rep.append(dict(analysis=a.analysis, chain=r['label'], seed=r['seed'], continued='no', saved_state_from_which_the_chain_was_continued='', state_at_which_the_first_part_was_ended='',
                                logged_states_of_the_first_part=n, logged_states_of_the_first_part_kept=n, logged_states_of_the_first_part_left_out=0, logged_states_of_the_continuation='', last_state_of_the_cut_chain=rows[-1][0],
                                logged_states_discarded_as_burn_in=nb, states_used_in_the_cut_summary=n - nb, sha256_of_the_log_of_the_cut_summary=sha(dst), log_of_the_cut_summary_is='a copy of the log of the chain'))
                print(f"{r['label']}: not continued, log copied unchanged ({n} logged states)")
                continue
            ls = [f for f in os.listdir(parts[1]) if re.fullmatch(r'loaded_state_(\d+)\.state', f)]
            assert len(ls) == 1, (parts[1], ls)
            N = int(re.fullmatch(r'loaded_state_(\d+)\.state', ls[0]).group(1))
            assert CS.saved_state(os.path.join(parts[1], ls[0])) == N, (parts[1], N)
            head, cols, rows, cut = read_log(src)
            keep = [x for x in rows if x[0] <= N]
            want = list(range(0, N + 1, a.log_every)); got = [x[0] for x in keep]
            if got != want or N % a.log_every != 0:
                raise SystemExit(f"STOP: the first part of {r['label']} does not hold every state once from 0 to {N} in steps of {a.log_every}: {len(got)} rows")
            with open(dst, 'w', newline='') as fh:
                fh.write(''.join(head)); fh.write(''.join(x[1] for x in keep))
            # the cut log must be the beginning of the log of the first part, byte for byte
            b_cut = open(dst, 'rb').read(); b_src = open(src, 'rb').read()
            assert b_src.startswith(b_cut) and b_cut.endswith(b'\n'), 'the cut log is not the beginning of the log of the first part'
            ncont = sum(len(read_log(os.path.join(p, stem + '.log'))[2]) for p in parts[1:])
            n = len(keep); nb = int(n * a.burnin)
            rep.append(dict(analysis=a.analysis, chain=r['label'], seed=r['seed'], continued=f'yes, {len(parts) - 1} continuation(s)', saved_state_from_which_the_chain_was_continued=N,
                            state_at_which_the_first_part_was_ended=rows[-1][0], logged_states_of_the_first_part=len(rows), logged_states_of_the_first_part_kept=n, logged_states_of_the_first_part_left_out=len(rows) - n,
                            logged_states_of_the_continuation=ncont, last_state_of_the_cut_chain=keep[-1][0], logged_states_discarded_as_burn_in=nb, states_used_in_the_cut_summary=n - nb,
                            sha256_of_the_log_of_the_cut_summary=sha(dst), log_of_the_cut_summary_is=f'the log of the first part up to state {N}, byte for byte its beginning'))
            print(f"{r['label']}: continued from state {N}; first part logged {len(rows)} states (to state {rows[-1][0]}), kept {n}, left out {len(rows) - n}; continuation logged {ncont} states, none used")
        with open(os.path.join(tmp, 'cut_report.csv'), 'w', newline='') as fh:
            if rep:
                w = csv.DictWriter(fh, fieldnames=list(rep[0].keys())); w.writeheader(); w.writerows(rep)
            else:
                fh.write('no finished chain\n')
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    shutil.rmtree(out, ignore_errors=True); os.replace(tmp, out)


if __name__ == '__main__':
    main()
