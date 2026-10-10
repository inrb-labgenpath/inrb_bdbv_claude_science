#!/usr/bin/env python3
"""build_sets_sensitivity.py (round of October 2026)

Builds the genome sets of the sensitivity analyses and of the removal sets from ONE snapshot. The same script is run on the
snapshot of 24 September 2026 (reproduction test) and on the snapshot of 2 October 2026.

Steps (the commands of steps 0 to 2 are those of round_20261002/regress_curation.sh of the working tree, with the names of the
files of the snapshot given as arguments; the scripts called are the frozen screen.py and build_alignments.py):
  0  restricted-use genomes of groups other than the main data provider (rule of run_curation_20260924.sh, step 0)
  1  first-stage screen (rules v1) and final screen (rules v2)
  2  build_alignments.py with the set specification given by --set-spec
  3  truncated sets: genomes of the reference set collected up to a date (--upto NAME=DATE)
  4  area of every genome 'by the recorded first administrative level' (rule of the code that wrote regional_set_20260924.csv,
     read from the recorded code of that table), after the removal of blanks at the ends of place names
  5  removal sets and controls (rule of the code that wrote the sets of 24 September, read from the recorded code of each)

Nothing is drawn or decided here beyond what the arguments say: the seeds of the draws, the reference period, the blocks and
the name of the reference set are arguments.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

ITURI_MAIN = ['Bunia', 'Rwampara', 'Nizi', 'Mongbwalu', 'Lita', 'Mangala', 'Bambu']      # as in the code of regional_set_20260924.csv


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, cwd, log):
    with open(log, 'a') as fh:
        fh.write('$ ' + ' '.join(cmd) + '\n')
        fh.flush()
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
        rc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode
    if rc != 0:
        raise SystemExit(f'command failed with exit code {rc}: see {log}')


def read_alignment(path):
    """records in the order of the file: list of (header, sequence)"""
    out = []
    for rec in open(path).read().split('>')[1:]:
        h, s = rec.split('\n', 1)
        out.append((h.strip(), s.replace('\n', '')))
    return out


def write_alignment(path, records, keep_accessions):
    keep = set(keep_accessions)
    n = 0
    with open(path, 'w') as fh:
        for h, s in records:
            if h.split('|')[0] in keep:
                fh.write(f'>{h}\n{s}\n')
                n += 1
    assert n == len(keep), (path, n, len(keep))
    return n


def strip_blanks(x):
    return x.strip() if isinstance(x, str) else x


def area_by_recorded_admin1(country, admin1, zone):
    """rule of the code that wrote regional_set_20260924.csv (function area_admin1), blanks at the ends removed before"""
    country, admin1, zone = strip_blanks(country), strip_blanks(admin1), strip_blanks(zone)
    if country == 'Uganda':
        return 'Uganda'
    if country != 'Democratic Republic of the Congo':
        return None
    if admin1 == 'Nord-Kivu':
        return 'Nord-Kivu'
    if admin1 == 'Sud-Kivu':
        return 'Sud-Kivu'
    if admin1 == 'Ituri':
        return zone if zone in ITURI_MAIN else 'other Ituri'
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--label', required=True, help='suffix of the output files, e.g. 20260924 or 20261002')
    ap.add_argument('--pipeline', required=True, help='directory with the frozen screen.py and build_alignments.py')
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--alignment', required=True)
    ap.add_argument('--validation', required=True)
    ap.add_argument('--duplicate-pairs', required=True)
    ap.add_argument('--amplicons', required=True)
    ap.add_argument('--previous-screen', required=True)
    ap.add_argument('--exclusion-list', required=True)
    ap.add_argument('--expected-n', type=int, required=True)
    ap.add_argument('--set-spec', required=True)
    ap.add_argument('--reference-set', required=True, help='name of the set of the specification that the truncated sets and the removal sets start from')
    ap.add_argument('--upto', action='append', default=[], help='NAME=YYYY-MM-DD: genomes of the reference set collected up to the date')
    ap.add_argument('--area', default='Bunia')
    ap.add_argument('--reference-period', nargs=2, required=True, metavar=('FIRST_DAY', 'DAY_AFTER_THE_LAST'))
    ap.add_argument('--blocks', nargs='+', required=True, help='FIRST_DAY:DAY_AFTER_THE_LAST for every block')
    ap.add_argument('--recent-from', required=True, help='first day of the genomes that the removal sets and controls concern')
    ap.add_argument('--thinned', action='append', default=[], help='NAME=SEED: draw of the thinned set')
    ap.add_argument('--without-recent', required=True, help='NAME of the set without the genomes of the area collected from --recent-from')
    ap.add_argument('--control-a', action='append', default=[], help='NAME=SEED: control of the thinned set (numbers per interval of the first thinned draw)')
    ap.add_argument('--control-b', action='append', default=[], help='NAME=SEED: control of the set without recent genomes of the area')
    ap.add_argument('--num-parameters', type=int, default=20)
    ap.add_argument('--cutoff', type=float, default=0.8)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--python', default=sys.executable)
    a = ap.parse_args()
    ab = os.path.abspath
    out = ab(a.outdir)
    for d in ('run_final', 'run_v2', 'run_align', 'sets', 'logs'):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    sfx = a.label
    log = os.path.join(out, 'logs', 'pipeline.log')
    open(log, 'w').close()
    A24 = ['--metadata', ab(a.metadata), '--alignment', ab(a.alignment)]
    AMP = ['--amplicons', ab(a.amplicons)]
    PREV = ['--previous-screen', ab(a.previous_screen)]
    V1SW = ['--undated-comparators', 'exclude', '--threshold-precision', 'exact', '--batch-rule', 'holm-exact', '--loo-scan', '--extended-columns']

    # ---- 0. restricted-use genomes (rule of run_curation_20260924.sh, step 0) ----
    m = pd.read_csv(a.metadata, low_memory=False)
    o = m[m.outbreak2026 == True]                                                            # noqa: E712 (as in the earlier script)
    main_group = o.groupName.value_counts().index[0]
    d7 = o[(o.dataUseTerms == 'RESTRICTED') & (o.groupName != main_group)]
    assert len(d7) == 2 and d7.groupName.nunique() == 2, 'expected two groups with one restricted-use genome each'
    D7 = ','.join(d7.accessionVersion)

    # ---- 1. screens ----
    scr = ab(os.path.join(a.pipeline, 'screen.py'))
    run([a.python, scr] + A24 + AMP + ['--outdir', '.', '--out-suffix', sfx, '--expected-n', str(a.expected_n)] + PREV + V1SW
        + ['--public-table', '--validation', ab(a.validation), '--withhold', D7, '--withhold-name', 'minusD7'], os.path.join(out, 'run_final'), log)
    run([a.python, scr] + A24 + AMP + ['--outdir', '.', '--out-suffix', sfx, '--expected-n', str(a.expected_n)] + PREV + ['--rules', 'v2',
        '--attach-class', f'rules_v1=../run_final/public_quality_screen_{sfx}.csv',
        '--public-table', '--validation', ab(a.validation), '--withhold', D7, '--withhold-name', 'minusD7'], os.path.join(out, 'run_v2'), log)

    # ---- 2. alignments of the specification ----
    spec = json.load(open(a.set_spec))
    spec_copy = os.path.join(out, 'run_align', os.path.basename(a.set_spec))
    with open(spec_copy, 'wb') as fh:
        fh.write(open(a.set_spec, 'rb').read())
    run([a.python, ab(os.path.join(a.pipeline, 'build_alignments.py'))] + A24 + ['--screen', f'../run_v2/public_quality_screen_{sfx}.csv',
        '--out-suffix', sfx, '--outdir', '.', '--set-spec', os.path.basename(a.set_spec), '--completeness-precision', 'exact',
        '--duplicate-pairs', ab(a.duplicate_pairs), '--data-use-exclude', D7, '--batch-table', f'../run_v2/screen_batches_{sfx}.csv',
        '--validation', ab(a.validation), '--exclusion-list', ab(a.exclusion_list), '--tail-days', '42'], os.path.join(out, 'run_align'), log)
    rows = []
    for name in spec['sets']:
        src = os.path.join(out, 'run_align', f'{name}_{sfx}.fasta')
        dst = os.path.join(out, 'sets', f'{name}_{sfx}.fasta')
        with open(dst, 'wb') as fh:
            fh.write(open(src, 'rb').read())
        rows.append(dict(set=name, built_by='build_alignments.py (frozen), set specification', file=os.path.basename(dst)))

    # ---- reference set ----
    ref_file = os.path.join(out, 'sets', f'{a.reference_set}_{sfx}.fasta')
    REC = read_alignment(ref_file)
    ST = pd.read_csv(os.path.join(out, 'run_align', f'alignment_sets_{sfx}.csv'), low_memory=False)
    prim = ST[ST[f'in_{a.reference_set}'].astype(str).str.lower() == 'true'][['accessionVersion', 'collection_date']].copy()
    assert set(prim.accessionVersion) == set(h.split('|')[0] for h, _ in REC) and len(prim) == len(REC)
    prim['d'] = pd.to_datetime(prim.collection_date)

    # ---- 3. truncated sets ----
    for item in a.upto:
        name, date = item.split('=')
        keep = list(prim.accessionVersion[prim.d <= pd.Timestamp(date)])
        write_alignment(os.path.join(out, 'sets', f'{name}_{sfx}.fasta'), REC, keep)
        rows.append(dict(set=name, built_by=f'genomes of the set {a.reference_set} collected up to {date}', file=f'{name}_{sfx}.fasta'))

    # ---- 4. area by the recorded first administrative level ----
    V = pd.read_csv(a.validation, low_memory=False).set_index('accessionVersion')
    AR = pd.DataFrame({'accessionVersion': ST.accessionVersion.values})
    AR['country'] = V['geoLocCountry'].reindex(AR.accessionVersion).values
    AR['recorded_admin1'] = V['geoLocAdmin1'].reindex(AR.accessionVersion).values
    AR['health_zone'] = V['zone'].reindex(AR.accessionVersion).values
    AR['area_by_recorded_admin1'] = [area_by_recorded_admin1(c, a1, z) for c, a1, z in zip(AR.country, AR.recorded_admin1, AR.health_zone)]
    AR['place_name_with_blank_at_an_end'] = [any(isinstance(x, str) and x != x.strip() for x in t) for t in zip(AR.country, AR.recorded_admin1, AR.health_zone)]
    AR[f'in_{a.reference_set}'] = AR.accessionVersion.isin(set(prim.accessionVersion)).values
    AR.to_csv(os.path.join(out, f'area_by_recorded_admin1_{sfx}.csv'), index=False)

    # ---- 5. removal sets and controls (code of the sets of 24 September, with the dates and seeds as arguments) ----
    pr = prim.merge(AR[['accessionVersion', 'area_by_recorded_admin1']], on='accessionVersion', how='left').reset_index(drop=True)
    pr['is_area'] = pr.area_by_recorded_admin1.eq(a.area)
    prx = pr.set_index('accessionVersion')
    r0, r1 = pd.Timestamp(a.reference_period[0]), pd.Timestamp(a.reference_period[1])
    ref = prx[(prx.d >= r0) & (prx.d < r1)]
    share = float(ref.is_area.mean())
    blocks = [tuple(b.split(':')) for b in a.blocks]
    recent_from = pd.Timestamp(a.recent_from)
    anchor = prx.d.max()
    iv = a.cutoff * 365 / (a.num_parameters - 1)
    pm = pr[['accessionVersion', 'collection_date', 'area_by_recorded_admin1', 'd']].copy()
    pm['interval'] = np.floor((anchor - pm.d).dt.days / iv).astype(int) + 1
    interval_of = dict(zip(pm.accessionVersion, pm.interval))
    memb = pd.DataFrame({'accessionVersion': list(prx.index), 'collection_date': prx.collection_date.values,
                         'area_by_recorded_admin1': prx.area_by_recorded_admin1.values, 'skygrid_interval': [interval_of[x] for x in prx.index]})
    draws = []
    plan_first = None
    for item in a.thinned:
        name, seed = item.split('=')
        rng = np.random.default_rng(int(seed))
        drop = []
        for b0, b1 in blocks:
            blk = prx[(prx.d >= pd.Timestamp(b0)) & (prx.d < pd.Timestamp(b1))]
            nb = int(blk.is_area.sum())
            no = int((~blk.is_area).sum())
            target = int(round(no * share / (1 - share)))
            bun = sorted(blk.index[blk.is_area])
            keepb = set(rng.choice(bun, size=min(target, nb), replace=False))
            drop += [x for x in bun if x not in keepb]
            draws.append(dict(set=name, seed_of_the_draw=int(seed), block_first_day=b0, block_day_after_the_last=b1, genomes_of_the_area_in_the_block=nb,
                              other_genomes_in_the_block=no, share_of_the_area_in_the_reference_period=share, genomes_of_the_area_kept=min(target, nb),
                              genomes_of_the_area_removed=nb - min(target, nb)))
        keep = [x for x in prx.index if x not in set(drop)]
        write_alignment(os.path.join(out, 'sets', f'{name}_{sfx}.fasta'), REC, keep)
        memb[f'in_{name}'] = memb.accessionVersion.isin(set(keep)).values
        rows.append(dict(set=name, built_by=f'thinned set, seed of the draw {seed}', file=f'{name}_{sfx}.fasta'))
        if plan_first is None:
            cnt = pd.Series([interval_of[x] for x in drop]).value_counts()
            plan_first = {int(k): int(cnt[k]) for k in sorted(cnt.index)}                # ascending interval, as the literal of the earlier code ({1: .., 2: ..})
    keep = list(prx.index[~(prx.is_area & (prx.d >= recent_from))])
    write_alignment(os.path.join(out, 'sets', f'{a.without_recent}_{sfx}.fasta'), REC, keep)
    memb[f'in_{a.without_recent}'] = memb.accessionVersion.isin(set(keep)).values
    rows.append(dict(set=a.without_recent, built_by=f'without the genomes of the area collected from {a.recent_from}', file=f'{a.without_recent}_{sfx}.fasta'))
    # controls
    recent = pm[pm.d >= recent_from].copy()
    last = pm[pm.d == pm.d.max()]
    protect = set(last.accessionVersion)
    MEM = memb[['accessionVersion', f'in_{a.without_recent}']].merge(pm[['accessionVersion', 'interval']], on='accessionVersion')
    rem_wo = MEM[MEM[f'in_{a.without_recent}'].astype(str).str.lower() != 'true']
    vc = rem_wo.interval.value_counts()
    rem_wo_int = {int(k): int(v) for k, v in vc.items()}                                  # order of value_counts, as in the earlier code
    ties_in_plan_b = bool(vc.duplicated(keep=False).any())

    def draw(seed, plan):
        rng_ = np.random.default_rng(seed)
        dropped = []
        for k, nrem in plan.items():
            cand = sorted(set(recent[recent.interval == k].accessionVersion) - protect)
            dropped += list(rng_.choice(cand, size=nrem, replace=False))
        return set(pm.accessionVersion) - set(dropped), dropped

    for group, items, plan in (('A', a.control_a, plan_first), ('B', a.control_b, rem_wo_int)):
        for item in items:
            name, seed = item.split('=')
            keep, dropped = draw(int(seed), plan)
            write_alignment(os.path.join(out, 'sets', f'{name}_{sfx}.fasta'), REC, keep)
            memb[f'in_{name}'] = memb.accessionVersion.isin(keep).values
            rows.append(dict(set=name, built_by=f'control {group}, seed of the draw {seed}, numbers removed per Skygrid interval {json.dumps(plan)}', file=f'{name}_{sfx}.fasta'))
            for k, nrem in plan.items():
                draws.append(dict(set=name, seed_of_the_draw=int(seed), skygrid_interval=k, genomes_removed=nrem,
                                  candidates_in_the_interval=len(set(recent[recent.interval == k].accessionVersion) - protect)))
    memb.to_csv(os.path.join(out, f'removal_sets_membership_{sfx}.csv'), index=False)
    pd.DataFrame(draws).to_csv(os.path.join(out, f'removal_sets_draws_{sfx}.csv'), index=False)

    # ---- summary ----
    S = []
    for r in rows:
        p = os.path.join(out, 'sets', r['file'])
        rec = read_alignment(p)
        d = pd.to_datetime(pd.Series([h.split('|')[-1] for h, _ in rec]))
        S.append(dict(r, genomes=len(rec), positions=len(rec[0][1]), first_collection_date=str(d.min().date()), most_recent_collection_date=str(d.max().date()),
                      genomes_on_the_most_recent_collection_date=int((d == d.max()).sum()), sha256=sha256(p)))
    pd.DataFrame(S).to_csv(os.path.join(out, f'sets_built_{sfx}.csv'), index=False)
    json.dump(dict(label=sfx, restricted_use_genomes_left_out=len(d7), reference_set=a.reference_set, area=a.area, reference_period=a.reference_period,
                   share_of_the_area_in_the_reference_period=share, genomes_in_the_reference_period=int(len(ref)),
                   genomes_of_the_area_in_the_reference_period=int(ref.is_area.sum()), blocks=blocks, recent_from=a.recent_from,
                   anchor_of_the_grid=str(anchor.date()), interval_days=iv, genomes_on_the_most_recent_collection_date=int(len(last)),
                   numbers_removed_per_interval_control_A=plan_first, numbers_removed_per_interval_control_B=rem_wo_int,
                   equal_numbers_in_the_plan_of_control_B=ties_in_plan_b,
                   arguments=vars(a)), open(os.path.join(out, f'sets_parameters_{sfx}.json'), 'w'), indent=1, default=str)
    print(f'{sfx}: {len(S)} sets written')
    return 0


if __name__ == '__main__':
    sys.exit(main())
