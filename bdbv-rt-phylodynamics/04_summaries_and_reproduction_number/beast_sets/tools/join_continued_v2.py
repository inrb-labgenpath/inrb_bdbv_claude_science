#!/usr/bin/env python3
"""join_continued_v2.py - input of the analysis for chains that were CONTINUED from a saved state. Reads the chain directories;
changes nothing in them.
VERSION 2 of 4 October 2026 (tools/join_continued.py is kept unchanged beside it: it wrote the first version of the joined input of beast_A_C00).
What is changed: the joined TREE FILE. Version 1 left out of the head every line 'End;', and so the line that closes the block 'Begin taxa;'.
Version 2 writes the head as it stands in the first part up to the
first tree, and after the last tree the lines that the program wrote after the last tree of the LAST part; no line of the tree file is
written by this tool itself. It STOPS if the heads of the parts differ (the tree lines refer to the numbers of the block 'Translate') or if
the program did not close the tree file of the last part. The joined logs, join_report.csv and checkpoint_comparison.csv are written as by
version 1. New: tree_files_report.csv, one row per chain.
usage: python tools/join_continued_v2.py <root> <analysis> <outdir> [--steps 40000000 --log-every 20000 --tree-every 400000]
For every FINISHED chain of the analysis it writes <outdir>/<analysis>/<label>/<stem>.log and <stem>.trees:
  chain without a continuation: copies of the log and of the tree file of the chain directory, unchanged (SHA-256 compared);
  continued chain: COPIES of the log and tree file of every part, cut at the state of the checkpoint that the next part loaded, and joined:
      part 1 holds the states 0 .. N1, continuation 1 the states above N1 (.. N2), and so on. The state of a checkpoint is taken from the
      part BEFORE the checkpoint (the chain as it was before it was ended).
  The script STOPS (exit code 1, no file for the analysis is left) unless the joined log holds every state once, in steps of
  --log-every from 0 to --steps, and the joined tree file every state once in steps of --tree-every.
<outdir>/<analysis>/join_report.csv   one row per chain: parts, states of the checkpoints, rows taken from each part
<outdir>/<analysis>/checkpoint_comparison.csv   where the part before and the continuation both log the state of a checkpoint: one row
                                                per column with the two values and their difference
prints one line per chain."""
import argparse, csv, hashlib, os, re, shutil, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chain_status as CS


RX_LEN = r'chain length of the XML of the continuation (\d+)'
RX_LNL = r'difference of the log-likelihoods (\S+)'


def read_log(path):
    head, cols, rows = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                head.append(line); continue
            if cols is None:
                cols = line.rstrip('\n').split('\t'); head.append(line); continue
            if not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != len(cols) or not f[0].isdigit():
                continue                                  # a last line that was cut when the chain was ended
            rows.append((int(f[0]), line if line.endswith('\n') else line + '\n'))
    return head, cols, rows


RX_TREE = re.compile(r'\s*tree STATE_(\d+)\b')
RX_END = re.compile(r'\s*End;\s*$', re.I)


def read_trees(path):
    """head: every line before the first tree, as it stands; trees: the complete tree lines; tail: every line after the last tree line,
    as it stands; cut: number of tree lines that were cut when the chain was ended (left out)."""
    head, trees, tail, cut, seen = [], [], [], 0, False
    with open(path, newline='') as fh:
        for line in fh:
            m = RX_TREE.match(line)
            if m:
                seen = True
                assert not tail, f'{path}: a tree line stands after a line that is not a tree'
                if line.rstrip().endswith(';') and line.endswith('\n'):
                    trees.append((int(m.group(1)), line))
                else:
                    cut += 1
            elif not seen:
                head.append(line)
            else:
                tail.append(line)
    return head, trees, tail, cut


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def parts_of(d, label):
    parts = [d]
    k = 1
    while os.path.isdir(os.path.join(d, f'continuation_{k}')):
        parts.append(os.path.join(d, f'continuation_{k}')); k += 1
    cks = []
    for p in parts[1:]:
        ls = [f for f in os.listdir(p) if re.fullmatch(r'loaded_state_(\d+)\.state', f)]
        assert len(ls) == 1, (p, ls)
        n = int(re.fullmatch(r'loaded_state_(\d+)\.state', ls[0]).group(1))
        assert CS.saved_state(os.path.join(p, ls[0])) == n, (p, n)
        assert not os.path.exists(os.path.join(p, 'continuation_rejected.txt')), p
        cks.append(n)
    return parts, cks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root'); ap.add_argument('analysis'); ap.add_argument('outdir')
    ap.add_argument('--steps', type=int, default=40_000_000); ap.add_argument('--log-every', type=int, default=20_000); ap.add_argument('--tree-every', type=int, default=400_000)
    ap.add_argument('--runs', default='runs')
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    out = os.path.join(os.path.abspath(a.outdir), a.analysis)
    tmp = out + '.tmp'
    shutil.rmtree(tmp, ignore_errors=True); os.makedirs(tmp)
    q = [r for r in CS.read_queue(root) if r['analysis'] == a.analysis]
    rep, cmp_rows, trep = [], [], []
    try:
        for r in q:
            d = os.path.join(root, a.runs, a.analysis, r['label']); stem = r['xml_stem']
            if not os.path.isdir(d):
                continue
            parts, cks = parts_of(d, r['label'])
            last = parts[-1]
            if not (os.path.exists(os.path.join(last, 'exit_code.txt')) and open(os.path.join(last, 'exit_code.txt')).read().strip() not in ('', '129', '130', '137', '143')):
                continue                                  # not finished
            t = os.path.join(tmp, r['label']); os.makedirs(t)
            if len(parts) == 1:
                for ext in ('.log', '.trees'):
                    shutil.copy2(os.path.join(d, stem + ext), os.path.join(t, stem + ext))
                    assert sha(os.path.join(d, stem + ext)) == sha(os.path.join(t, stem + ext))
                rep.append(dict(chain=r['label'], seed=r['seed'], parts=1, checkpoints='', rows_from_each_part='all', trees_from_each_part='all', continued='no',
                                log_states=len(read_log(os.path.join(t, stem + '.log'))[2]), log_last_state=read_log(os.path.join(t, stem + '.log'))[2][-1][0]))
                hj, tj, lj, cj = read_trees(os.path.join(t, stem + '.trees'))
                trep.append(dict(chain=r['label'], seed=r['seed'], parts=1, lines_of_the_head=len(hj), head_identical_in_all_parts='one part',
                                 lines_End_in_the_head=sum(1 for x in hj if RX_END.match(x)), trees=len(tj), trees_from_each_part='all',
                                 tree_lines_cut_when_a_part_was_ended_and_left_out=str(cj), lines_after_the_last_tree=len(lj),
                                 lines_End_after_the_last_tree=sum(1 for x in lj if RX_END.match(x)), sha256_of_the_joined_tree_file=sha(os.path.join(t, stem + '.trees')),
                                 sha256_of_the_joined_log=sha(os.path.join(t, stem + '.log'))))
                print(f"{r['label']}: 1 part, copied unchanged")
                continue
            assert all(cks[i] < cks[i + 1] for i in range(len(cks) - 1)), cks
            logs = [read_log(os.path.join(p, stem + '.log')) for p in parts]
            trs = [read_trees(os.path.join(p, stem + '.trees')) for p in parts]
            assert all(l[1] == logs[0][1] for l in logs), 'columns of the logs differ'
            if not all(t_[0] == trs[0][0] for t_ in trs):
                raise SystemExit(f"STOP: the heads of the tree files of the parts of {r['label']} differ")
            if sum(1 for x in trs[-1][2] if RX_END.match(x)) != 1 or any(x.strip() and not RX_END.match(x) for x in trs[-1][2]):
                raise SystemExit(f"STOP: the tree file of the last part of {r['label']} is not closed by one line 'End;' after its last tree")
            bounds = [-1] + cks + [a.steps]              # part i holds the states above bounds[i] up to bounds[i+1]
            jl, jt, nrow, ntree = [], [], [], []
            for i in range(len(parts)):
                lo, hi = bounds[i], bounds[i + 1]
                rr = [x for x in logs[i][2] if lo < x[0] <= hi]
                tt = [x for x in trs[i][1] if lo < x[0] <= hi]
                jl += rr; jt += tt; nrow.append(len(rr)); ntree.append(len(tt))
                if i > 0:                                 # the state of the checkpoint in the part before and in this continuation
                    n = cks[i - 1]
                    before = [x for x in logs[i - 1][2] if x[0] == n]; after = [x for x in logs[i][2] if x[0] == n]
                    if before and after:
                        fb, fa = before[0][1].rstrip('\n').split('\t'), after[0][1].rstrip('\n').split('\t')
                        for c, vb, va in zip(logs[0][1], fb, fa):
                            if c == 'state':
                                continue
                            cmp_rows.append(dict(analysis=a.analysis, chain=r['label'], seed=r['seed'], continuation=i, state_of_the_checkpoint=n, column=c, value_in_the_part_before=vb,
                                                 value_in_the_continuation=va, difference=repr(float(va) - float(vb)), identical_text=(va == vb)))
                    else:
                        cmp_rows.append(dict(analysis=a.analysis, chain=r['label'], seed=r['seed'], continuation=i, state_of_the_checkpoint=n, column='(all)', value_in_the_part_before='' if not before else 'logged',
                                             value_in_the_continuation='' if not after else 'logged', difference='', identical_text='not compared: the state of the checkpoint is not logged in both parts'))
            want = list(range(0, a.steps + 1, a.log_every)); got = [x[0] for x in jl]
            if got != want:
                raise SystemExit(f"STOP: joined log of {r['label']} does not hold every state once from 0 to {a.steps} in steps of {a.log_every}: {len(got)} rows, first difference at "
                                 f"{next((w for w, g_ in zip(want, got + [None] * len(want)) if w != g_), None)}")
            wantt = list(range(0, a.steps + 1, a.tree_every)); gott = [x[0] for x in jt]
            if gott != wantt:
                raise SystemExit(f"STOP: joined tree file of {r['label']} does not hold every state once from 0 to {a.steps} in steps of {a.tree_every}: {len(gott)} trees")
            with open(os.path.join(t, stem + '.log'), 'w') as fh:
                fh.write(''.join(logs[0][0])); fh.write(''.join(x[1] for x in jl))
            with open(os.path.join(t, stem + '.trees'), 'w', newline='') as fh:
                fh.write(''.join(trs[0][0])); fh.write(''.join(x[1] for x in jt)); fh.write(''.join(trs[-1][2]))
            # the joined tree file without its tree lines must be the tree file of the last part without its tree lines
            hj, tj, lj, cj = read_trees(os.path.join(t, stem + '.trees'))
            assert hj == trs[-1][0] and lj == trs[-1][2] and cj == 0 and [x[0] for x in tj] == wantt, 'joined tree file is not as written'
            trep.append(dict(chain=r['label'], seed=r['seed'], parts=len(parts), lines_of_the_head=len(hj), head_identical_in_all_parts='yes',
                             lines_End_in_the_head=sum(1 for x in hj if RX_END.match(x)), trees=len(tj), trees_from_each_part=' '.join(map(str, ntree)),
                             tree_lines_cut_when_a_part_was_ended_and_left_out=' '.join(str(x[3]) for x in trs), lines_after_the_last_tree=len(lj),
                             lines_End_after_the_last_tree=sum(1 for x in lj if RX_END.match(x)), sha256_of_the_joined_tree_file=sha(os.path.join(t, stem + '.trees')),
                             sha256_of_the_joined_log=sha(os.path.join(t, stem + '.log'))))
            if logs[-1][2][-1][0] != a.steps:
                raise SystemExit(f"STOP: the last part of {r['label']} ends at state {logs[-1][2][-1][0]}, not at {a.steps}")
            det = []
            for i, pth in enumerate(parts[1:], 1):
                info = open(os.path.join(pth, 'continuation_info.txt')).read() if os.path.exists(os.path.join(pth, 'continuation_info.txt')) else ''
                g = lambda rx: (re.search(rx, info).group(1) if re.search(rx, info) else 'not recorded')
                xml_len = g(RX_LEN); lnl = g(RX_LNL)
                det.append(f"continued from state {cks[i - 1]} (the part before was ended at state {logs[i - 1][2][-1][0]}; chain length of the XML of the continuation "
                           f"{xml_len}; difference of the log-likelihoods {lnl})")
            rep.append(dict(chain=r['label'], seed=r['seed'], parts=len(parts), checkpoints=' '.join(map(str, cks)), rows_from_each_part=' '.join(map(str, nrow)),
                            trees_from_each_part=' '.join(map(str, ntree)), continued='; '.join(det), log_states=len(jl), log_last_state=jl[-1][0]))
            print(f"{r['label']}: {len(parts)} parts, checkpoints {cks}, rows {nrow}, trees {ntree}: joined")
        for name, rows in (('join_report.csv', rep), ('checkpoint_comparison.csv', cmp_rows), ('tree_files_report.csv', trep)):
            with open(os.path.join(tmp, name), 'w', newline='') as fh:
                if rows:
                    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
                else:
                    fh.write('no row: no chain of this analysis was continued\n' if name.startswith('checkpoint') else 'no finished chain\n')
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    shutil.rmtree(out, ignore_errors=True); os.replace(tmp, out)


if __name__ == '__main__':
    main()
