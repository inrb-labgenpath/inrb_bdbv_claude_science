#!/usr/bin/env python3
"""postprocess_sensitivity_v3.py - derives ONE analysis from its SAVED archive of raw output and checks it.

postprocess_sensitivity_v2.py with ONE change: it calls derive_sensitivity_v3.py (same code as v2, head corrected). v2 derived the first analysis
(delphy_sens_allClasses); v3 derives the analyses after it.
Against postprocess_sensitivity.py: the rows without a criterion are given the verdict 'not assessed' in place of 'yes'; the first version wrote
two checks of a SHA-256 as passed without comparing them with anything, here the archive is compared with the saved copy of the archive (its path
is an argument), the alignment is compared with its row of genome_sets_sensitivity_20261002.csv, and the check of the degenerate state names the
precision as fixed only where the log shows it.

  process(analysis, periods_csv, schedule_csv, stored_archive=<path of the saved archive>)

Reads   run/archives/<analysis>_raw_output_20261002.tar.gz, delivery/commands_sensitivity_20261002.csv,
        run/aln/<set>_20261002.fasta, the table of the periods, the schedule
Writes  delivery/analyses/<analysis>/{chains,estimates,intervals,trajectories,root_date_intervals}_<analysis>_20261002.csv,
        posterior_<analysis>_20261002.csv.gz, checks_<analysis>_20261002.csv, recomputation_<analysis>_20261002.csv, halves_<analysis>_20261002.csv
"""
import csv
import hashlib
import json
import os
import shlex
import shutil
import sys
import tarfile
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_sensitivity as CK          # noqa: E402
import derive_sensitivity_v3 as DS      # noqa: E402

LABEL = '20261002'
PIPE = 'tree/pipeline_round_20261002/work/pipeline_ext'
STEPS = 1500000000
FIXED_END = {'P1': '2026-05-15', 'P2': '2026-06-15', 'P4': '2026-08-15', 'P5': '2026-08-01'}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def options(cmd):
    t = shlex.split(cmd)
    o = {}
    i = 1
    while i < len(t):
        if t[i].startswith('--') and i + 1 < len(t) and not t[i + 1].startswith('--'):
            o[t[i]] = t[i + 1]
            i += 2
        else:
            o[t[i]] = None
            i += 1
    return o


def setting_of(analysis, description):
    if 'smooth' in analysis:
        return 'Skygrid, ' + description
    return 'Skygrid, smoothing fixed (program default); ' + description


def process(analysis, periods_csv, schedule_csv, workdir='.', stored_archive=None):
    t_begin = time.time()
    cmds = pd.read_csv(os.path.join(workdir, 'delivery', f'commands_sensitivity_{LABEL}.csv'), dtype=str, keep_default_na=False)
    mine = cmds[cmds['analysis'] == analysis].sort_values('seed')
    assert len(mine) > 0, analysis
    sched = pd.read_csv(schedule_csv, dtype=str, keep_default_na=False)
    srow = sched[sched['analysis'] == analysis]
    assert list(srow['seed'].sort_values()) == list(mine['seed']), 'the seeds of the commands are not those of the schedule'
    description = srow['description'].iloc[0]
    genome_set = mine['genome_set'].iloc[0]
    fasta = os.path.join(workdir, 'run', mine['alignment'].iloc[0])
    arch = os.path.join(workdir, 'run', 'archives', f'{analysis}_raw_output_{LABEL}.tar.gz')
    src = os.path.join(workdir, 'derive_in', analysis)
    if os.path.isdir(src):
        shutil.rmtree(src)
    os.makedirs(os.path.join(workdir, 'derive_in'), exist_ok=True)
    with tarfile.open(arch) as tf:
        tf.extractall(os.path.join(workdir, 'derive_in'))
    o = options(mine['command_line'].iloc[0])
    K = int(o['--v0-skygrid-num-parameters'])
    first_knot, last_knot = o.get('--v0-skygrid-first-knot-date') or '', o.get('--v0-skygrid-last-knot-date') or ''
    cutoff = float(o['--v0-skygrid-cutoff']) if '--v0-skygrid-cutoff' in o else None
    assert int(o['--v0-steps']) == STEPS
    led = list(csv.DictReader(open(os.path.join(src, f'ledger_rows_{analysis}.csv'), newline='')))
    outdir = os.path.join(workdir, 'delivery', 'analyses', analysis)
    out = DS.derive(analysis, setting_of(analysis, description), genome_set, list(mine['seed']), STEPS, src, fasta, PIPE, periods_csv, outdir, ledger_rows=led,
                    num_parameters=K, cutoff=cutoff if cutoff is not None else 0.8, first_knot=first_knot, last_knot=last_knot)
    DS.write_tables(out, analysis, outdir)
    san = os.path.join(outdir, 'copies_without_the_cut_line')
    ch = out['chains']
    an = out['analysis']
    g = an.grid
    checks = []

    def check(name, expected, observed, passed):
        checks.append(dict(analysis=analysis, check=name, expected=str(expected), observed=str(observed), passed=passed))
    if stored_archive is not None and os.path.exists(stored_archive):
        check('archive of the raw output: SHA-256 of the file from which the tables are derived against the file saved in the artifact store', sha256(stored_archive), sha256(arch),
              'yes' if sha256(stored_archive) == sha256(arch) else 'no')
    else:
        check('archive of the raw output: SHA-256 of the file from which the tables are derived against the file saved in the artifact store', 'not assessed: no path of the saved file was given',
              sha256(arch), 'not assessed')
    state = open(os.path.join(src, 'ARCHIVE_STATE.txt')).read()
    check('archive of the raw output is complete (ARCHIVE_STATE.txt)', 'complete: yes', [x for x in state.split('\n') if x.startswith('complete')][0], 'yes' if 'complete: yes' in state else 'no')
    man = pd.read_csv(os.path.join(src, 'MANIFEST.csv'), dtype=str)
    bad = [m for m, s_ in zip(man['member'], man['sha256']) if sha256(os.path.join(workdir, 'derive_in', m)) != s_]
    check('members of the archive against its MANIFEST.csv (SHA-256)', f'{len(man)} members agree', f'{len(man) - len(bad)} of {len(man)} agree', 'yes' if not bad else 'no')
    check('every chain of the schedule has a row in the table of chains, and no other chain has', ' '.join(srow['seed'].sort_values()), ' '.join(ch['seed'].astype(str)),
          'yes' if list(ch['seed'].astype(str)) == list(srow['seed'].sort_values()) else 'no')
    same_cmd = all(a == b for a, b in zip(ch['command_line'], mine['command_line']))
    check('command line of every chain (file .cmd) is the command of the table of commands', f'{len(mine)} of {len(mine)}', f'{sum(a == b for a, b in zip(ch["command_line"], mine["command_line"]))} of {len(mine)}', 'yes' if same_cmd else 'no')
    seeds_in_log = []
    for s in mine['seed']:
        head = open(os.path.join(src, f'{analysis}_s{s}.log'), errors='replace').read(4000)
        seeds_in_log.append(('# Seed: ' + s) in head or (f'--v0-seed {s}' in head))
    check('seed written by the program into the head of every log is the seed of the schedule', f'{len(mine)} of {len(mine)}', f'{sum(seeds_in_log)} of {len(mine)}', 'yes' if all(seeds_in_log) else 'no')
    check('limit of the address space of every chain (file .env)', '8388608 KiB', ' '.join(sorted(set(ch['address_space_limit_kib'].astype(str)))), 'yes' if set(ch['address_space_limit_kib'].astype(str)) == {'8388608'} else 'no')
    gs = pd.read_csv(os.path.join(workdir, 'delivery', f'genome_sets_sensitivity_{LABEL}.csv'), dtype=str, keep_default_na=False).set_index('set').loc[genome_set]
    exp_aln = f"{gs['sha256']}; {gs['genomes']} genomes; {gs['first_collection_date']} to {gs['most_recent_collection_date']}"
    obs_aln = f'{sha256(fasta)}; {an.n_tips} genomes; {an.first_tip.date()} to {an.last_tip.date()}'
    check('alignment of the analysis against its row of genome_sets_sensitivity_20261002.csv: SHA-256, genomes, first and most recent collection date', exp_aln, obs_aln,
          'yes' if exp_aln == obs_aln else 'no')
    in_cmd = sorted(set(options(c)['--v0-in-fasta'] for c in ch['command_line']))
    check('alignment named in the command line of every chain is the alignment of the analysis', mine['alignment'].iloc[0], ' '.join(in_cmd), 'yes' if in_cmd == [mine['alignment'].iloc[0]] else 'no')
    n_run = len(ch)
    counts = dict(chains_run=n_run, ended_by_the_program=int((ch['end_of_the_last_attempt_by_the_ledger'] == 'ended by the program').sum()),
                  started_again=int((ch['started_again'] == 'yes').sum()), did_not_reach_the_last_step=int((ch['reached_last_step'] == 'no').sum()),
                  root_date_rule_met=int((ch['root_date_rule_met'] == 'yes').sum()), density_criterion_met=int((ch['density_criterion_met'] == 'yes').sum()),
                  discarded_for_the_tree_of_step_0_alone=int((ch['discarded_for_the_tree_of_step_0_alone'] == 'yes').sum()), retained=int((ch['retained'] == 'yes').sum()))
    check('chains run, ended by the program, discarded by each rule, retained', 'no criterion is set: the numbers are reported and not acted on', json.dumps(counts), 'not assessed')
    deg = ch[(ch['retained'] == 'yes') & (ch['entered_degenerate_state_by_the_definition_of_the_brief'] == 'yes')]
    fixed_all = all(str(a_) == str(b_) and str(a_) not in ('', 'nan') for a_, b_ in zip(ch['smallest_precision'], ch['largest_precision']))
    check('precision of the Skygrid in the logs of the chains (smallest and largest logged value of every chain)', 'the same value in every logged state of every chain (the smoothing is fixed in every analysis of this track)',
          'fixed in every chain at ' + ' '.join(sorted(set(ch['smallest_precision'].astype(str)))) if fixed_all else 'not fixed in every chain: the running median of the precision is NOT assessed by this script', 'yes' if fixed_all else 'not assessed')
    check('a retained chain entered the degenerate state (root-date rule met in any logged state'
          + ('; the precision is fixed in every chain, so that its running median does not apply)' if fixed_all else '; running median of the precision not assessed)'),
          '0 retained chains', f'{len(deg)} retained chains', ('yes' if len(deg) == 0 else 'no') if fixed_all else 'not assessed')
    rc = hv = None
    if out.get('posterior') is not None:
        est = pd.read_csv(os.path.join(outdir, f'estimates_{analysis}_{LABEL}.csv'), dtype=str, keep_default_na=False)
        post = pd.read_csv(os.path.join(outdir, f'posterior_{analysis}_{LABEL}.csv.gz'))
        per = {}
        for P, d in out['periods'].items():
            per[P] = (d['start'], (g.t_anchor - 0.5 * g.D) if P in ('P3', 'P6') else CK.dec_year(FIXED_END[P]), d['state'] == 'evaluated')
        rc = CK.recompute(post, est, str(an.last_tip.date()), g.t_anchor, g.D, g.K, per)
        rc.insert(0, 'analysis', analysis)
        rc.to_csv(os.path.join(outdir, f'recomputation_{analysis}_{LABEL}.csv'), index=False)
        n_tie = int(rc['agrees'].astype(str).str.startswith('tie').sum())
        n_no = int((rc['agrees'] == 'no').sum())
        check('estimates recomputed from the posterior file by code that does not call the stored functions', f'{len(rc)} values agree (ties named)',
              f"{int((rc['agrees'] == 'yes').sum())} agree, {n_tie} ties, {n_no} do not agree", 'yes' if n_no == 0 else 'no')
        hv = CK.halves(post, est, str(an.last_tip.date()))
        hv.insert(0, 'analysis', analysis)
        hv.to_csv(os.path.join(outdir, f'halves_{analysis}_{LABEL}.csv'), index=False)
        col = 'difference_of_the_medians_in_per_cent_of_the_width_of_the_95_interval'
        check('the two halves of the retained chains: difference of the medians in per cent of the width of the 95 % interval (R of P1, P2, P3; rate; root date)', 'no criterion is set: the differences are reported',
              '; '.join(f"{q.replace('reproduction number, ', 'R ').replace('date of the most recent common ancestor (decimal year)', 'root date')}: "
                        + (f'{v:.2f}' if isinstance(v, float) else str(v)) for q, v in zip(hv['quantity'], hv[col])), 'not assessed')
        rows_post = len(post)
        check('rows of the posterior file = states after the burn-in of the retained chains', int(ch.loc[ch['retained'] == 'yes', 'n_states_after_burn_in'].astype(int).sum()), rows_post,
              'yes' if rows_post == int(ch.loc[ch['retained'] == 'yes', 'n_states_after_burn_in'].astype(int).sum()) else 'no')
    else:
        check('estimates recomputed from the posterior file by code that does not call the stored functions', 'values agree', 'not assessed: no chain retained, no posterior file', 'not assessed')
        check('the two halves of the retained chains', 'reported', 'not assessed: no chain retained', 'not assessed')
    C = pd.DataFrame(checks)
    C.to_csv(os.path.join(outdir, f'checks_{analysis}_{LABEL}.csv'), index=False)
    if os.path.isdir(san) and not os.listdir(san):
        os.rmdir(san)
    e = out['estimates'].set_index('quantity')
    summary = dict(analysis=analysis, seconds=round(time.time() - t_begin), **counts, checks_not_passed=int((C['passed'] == 'no').sum()),
                   estimates_not_meeting_the_criterion=int((out['estimates']['meets_convergence_criterion'] == 'no').sum()))
    return summary, out


if __name__ == '__main__':
    s, _ = process(sys.argv[1], sys.argv[2], sys.argv[3], stored_archive=(sys.argv[4] if len(sys.argv) > 4 else None))
    print(json.dumps(s))
