#!/usr/bin/env python3
"""set_tables_20261007.py (snapshot of 7 October 2026)

Writes into <out dir>, from the files that build_sets_sensitivity_v2.py wrote into sets07/ and from the directories of the tests:
  genome_sets_sensitivity_20261007.csv              one row per genome set (form of genome_sets_sensitivity_20261002.csv)
  genome_sets_membership_sensitivity_20261007.csv   one row per genome of the files that the program of the quality checks read
  removal_sets_parameters_sensitivity_20261007.json, removal_sets_draws_sensitivity_20261007.csv   as the builder wrote them
  alignments_sensitivity_20261007.tar.gz            the 19 alignments
  INTERNAL_reproduction_test_sensitivity_20261007.csv   reproduction tests and checks of the sets
It writes no file twice. Every number is read from a file at the time of the run.

usage: python3 code/set_tables_20261007.py --sets-dir sets07 --queue run/queue.csv --set-spec <json> --stored-alignments <tar.gz of 2 October>
       --stored-sets-table <csv of 2 October> --delivered-S <fasta> --quality-screen-table <csv> --genome-sets-table <csv>
       --test NAME=DIR=CODE=CALL (three times) --build-call specs/build_command_20261007.txt --out delivery
"""
import argparse, gzip, hashlib, io, json, os, tarfile, time
import pandas as pd

L = '20261007'
ap = argparse.ArgumentParser()
for o in ('--sets-dir', '--queue', '--set-spec', '--stored-alignments', '--stored-sets-table', '--delivered-S', '--quality-screen-table', '--genome-sets-table', '--build-call', '--build-code', '--out'):
    ap.add_argument(o, required=True)
ap.add_argument('--test', action='append', default=[])
a = ap.parse_args()
now = lambda: time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def sha(b):
    return hashlib.sha256(b).hexdigest()


def sha_file(p):
    return sha(open(p, 'rb').read())


def heads(b):
    return [l[1:].strip() for l in b.decode().split('\n') if l.startswith('>')]


def acc(b):
    return [h.split('|')[0] for h in heads(b)]


for f in ('genome_sets_sensitivity', 'genome_sets_membership_sensitivity', 'removal_sets_draws_sensitivity', 'INTERNAL_reproduction_test_sensitivity'):
    assert not os.path.exists(os.path.join(a.out, f'{f}_{L}.csv')), f
S = a.sets_dir
B = pd.read_csv(f'{S}/sets_built_{L}.csv', dtype=str, keep_default_na=False)
P = json.load(open(f'{S}/sets_parameters_{L}.json'))
Q = pd.read_csv(a.queue, dtype=str, keep_default_na=False)
spec = json.load(open(a.set_spec))
# ---- table of sets ----
GS = B.copy()
for r in GS.itertuples():
    assert sha_file(f'{S}/sets/{r.file}') == r.sha256
GS['alignment'] = 'aln/' + GS['file']
GS['analyses_that_run_on_it'] = [' '.join(sorted(Q[Q.genome_set == x].analysis.unique())) for x in GS['set']]
GS['chains'] = [int((Q.genome_set == x).sum()) for x in GS['set']]
GS['definition'] = [json.dumps(spec['sets'][x]) if x in spec['sets'] else b for x, b in zip(GS['set'], GS['built_by'])]
GS = GS[['set', 'genomes', 'positions', 'first_collection_date', 'most_recent_collection_date', 'genomes_on_the_most_recent_collection_date', 'sha256',
         'alignment', 'built_by', 'definition', 'analyses_that_run_on_it', 'chains']]
assert GS.chains.sum() == len(Q) == 222
GS.to_csv(os.path.join(a.out, f'genome_sets_sensitivity_{L}.csv'), index=False)
# ---- membership ----
AS = pd.read_csv(f'{S}/run_align/INTERNAL_alignment_sets_{L}.csv', dtype=str, keep_default_na=False)
RM = pd.read_csv(f'{S}/removal_sets_membership_{L}.csv', dtype=str, keep_default_na=False)
AR = pd.read_csv(f'{S}/area_by_recorded_admin1_{L}.csv', dtype=str, keep_default_na=False)
M = AS[['accessionVersion', 'collection_date', 'called_fraction_exact', 'screen_class', 'in_C00', 'in_allClasses', 'in_minusExcluded1046', 'in_retained1046',
        'in_dupNone', 'in_dupKeepBest', 'in_drcOnly']].copy()
for c in [c for c in M.columns if c.startswith('in_')]:
    M[c] = M[c].str.lower().map({'true': 'True', 'false': 'False'})
    assert M[c].notna().all()
members = {}
for r in B.itertuples():
    members[r.set] = acc(open(f'{S}/sets/{r.file}', 'rb').read())
    assert len(members[r.set]) == int(r.genomes) == len(set(members[r.set]))
for s_ in ('C00', 'allClasses', 'minusExcluded1046', 'retained1046', 'dupNone', 'dupKeepBest', 'drcOnly'):
    assert set(M[M[f'in_{s_}'] == 'True'].accessionVersion) == set(members[s_]), s_
for s_ in ('upTo0831', 'upTo0817'):
    M[f'in_{s_}'] = ['True' if x in set(members[s_]) else 'False' for x in M.accessionVersion]
M = M.merge(AR[['accessionVersion', 'area_by_recorded_admin1']], on='accessionVersion', how='left', validate='one_to_one')
M = M.merge(RM[['accessionVersion', 'skygrid_interval']].rename(columns={'skygrid_interval': 'skygrid_interval_of_S'}), on='accessionVersion', how='left', validate='one_to_one')
M['skygrid_interval_of_S'] = M['skygrid_interval_of_S'].fillna('')
for s_ in [c[3:] for c in RM.columns if c.startswith('in_')]:
    M[f'in_{s_}'] = ['True' if x in set(members[s_]) else 'False' for x in M.accessionVersion]
    assert set(RM[RM[f'in_{s_}'].str.lower() == 'true'].accessionVersion) == set(members[s_]), s_
assert len(M) == len(AS) and M.accessionVersion.is_unique
M.to_csv(os.path.join(a.out, f'genome_sets_membership_sensitivity_{L}.csv'), index=False)
# ---- removal sets: parameters and draws as written by the builder ----
open(os.path.join(a.out, f'removal_sets_parameters_sensitivity_{L}.json'), 'wb').write(open(f'{S}/sets_parameters_{L}.json', 'rb').read())
open(os.path.join(a.out, f'removal_sets_draws_sensitivity_{L}.csv'), 'wb').write(open(f'{S}/removal_sets_draws_{L}.csv', 'rb').read())
# ---- archive of the alignments ----
tp = os.path.join(a.out, f'alignments_sensitivity_{L}.tar.gz')
assert not os.path.exists(tp)
with open(tp, 'wb') as raw:
    with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode='w') as tf:
            for r in sorted(B.file):
                ti = tf.gettarinfo(f'{S}/sets/{r}', arcname=f'aln/{r}')
                ti.mtime = 0; ti.uid = ti.gid = 0; ti.uname = ti.gname = ''
                tf.addfile(ti, open(f'{S}/sets/{r}', 'rb'))
# ---- reproduction tests ----
stored = {}
with tarfile.open(a.stored_alignments) as t:
    for m_ in t.getmembers():
        if m_.isfile() and m_.name.endswith('.fasta'):
            stored[os.path.basename(m_.name)] = t.extractfile(m_).read()
T02 = pd.read_csv(a.stored_sets_table, dtype=str, keep_default_na=False).set_index('set')
rows = []
for item in a.test:
    name, d, code, call = item.split('=', 3)
    times = open(os.path.join(d, 'TEST_times.txt')).read().split('\n')
    for fn in sorted(stored):
        s_ = fn.replace('_20261002.fasta', '')
        p = os.path.join(d, 'sets02', 'sets', 'TEST_' + fn)
        rb = open(p, 'rb').read() if os.path.exists(p) else None
        a_s, a_r = acc(stored[fn]), (acc(rb) if rb is not None else [])
        only_s = sorted(set(a_s) - set(a_r)); only_r = sorted(set(a_r) - set(a_s))
        rows.append(dict(kind='reproduction of the sets of 2 October 2026', test=name, set=s_, genomes_stored_or_expected=len(a_s), genomes_built=len(a_r),
                         accessions_only_in_the_stored_set=len(only_s), accessions_only_in_the_built_set=len(only_r),
                         accessions_that_differ=' '.join(only_s + only_r), same_order=('yes' if a_s == a_r else 'no'),
                         bytes_stored_or_expected=len(stored[fn]), bytes_built=(len(rb) if rb is not None else ''), sha256_stored_or_expected=sha(stored[fn]),
                         sha256_built=(sha(rb) if rb is not None else ''), sha256_in_the_stored_table_of_sets=T02.loc[s_, 'sha256'],
                         byte_for_byte_identical=('yes' if rb == stored[fn] else 'no'), reproduced=('yes' if rb == stored[fn] else 'no'),
                         code=os.path.basename(code), sha256_of_the_code=sha_file(code), call=open(call).read().strip(), directory_of_the_test=d,
                         run_utc=' '.join(x for x in times if x), compared_utc=now()))
# ---- checks of the sets ----
bS = open(a.delivered_S, 'rb').read(); bC = open(f'{S}/sets/C00_{L}.fasta', 'rb').read()
dS = sorted(h.split('|')[1] for h in heads(bS))
common = dict(code=os.path.basename(a.build_code), sha256_of_the_code=sha_file(a.build_code), call=open(a.build_call).read().strip(), directory_of_the_test=S, run_utc='', compared_utc=now())
rows.append(dict(kind='sets of the update', test='reference set S built by this run against the alignment of S given in the call (C00_20261007.fasta)', set='C00',
                 genomes_stored_or_expected=len(acc(bS)), genomes_built=len(acc(bC)), accessions_only_in_the_stored_set=len(set(acc(bS)) - set(acc(bC))),
                 accessions_only_in_the_built_set=len(set(acc(bC)) - set(acc(bS))), accessions_that_differ='', same_order=('yes' if acc(bS) == acc(bC) else 'no'),
                 bytes_stored_or_expected=len(bS), bytes_built=len(bC), sha256_stored_or_expected=sha(bS), sha256_built=sha(bC), sha256_in_the_stored_table_of_sets='',
                 byte_for_byte_identical=('yes' if bS == bC else 'no'), reproduced=('yes' if bS == bC else 'no'), **common))
qs = open(a.quality_screen_table, 'rb').read(); ms = open(f'{S}/run_v2/INTERNAL_public_quality_screen_{L}.csv', 'rb').read()
rows.append(dict(kind='sets of the update', test='table of the program of the quality checks written by this run against the table of the quality checks given in the call', set='',
                 genomes_stored_or_expected=qs.count(b'\n') - 1, genomes_built=ms.count(b'\n') - 1, accessions_only_in_the_stored_set='', accessions_only_in_the_built_set='',
                 accessions_that_differ='', same_order='', bytes_stored_or_expected=len(qs), bytes_built=len(ms), sha256_stored_or_expected=sha(qs), sha256_built=sha(ms),
                 sha256_in_the_stored_table_of_sets='', byte_for_byte_identical=('yes' if qs == ms else 'no'), reproduced=('yes' if qs == ms else 'no'), **common))
G = pd.read_csv(a.genome_sets_table, dtype=str, keep_default_na=False)
g = set(G[G['in_C00'].str.lower().isin(['true', 'yes', '1'])].iloc[:, 0])
gc = [c for c in G.columns if 'accession' in c.lower()][0]
g = set(G[G['in_C00'].str.lower().isin(['true', 'yes', '1'])][gc])
rows.append(dict(kind='sets of the update', test='members of S by the table of genome sets of the update (genome_sets_by_called_fraction_20261007.csv) against the set built by this run', set='C00',
                 genomes_stored_or_expected=len(g), genomes_built=len(acc(bC)), accessions_only_in_the_stored_set=len(g - set(acc(bC))), accessions_only_in_the_built_set=len(set(acc(bC)) - g),
                 accessions_that_differ='', same_order='', bytes_stored_or_expected='', bytes_built='', sha256_stored_or_expected='', sha256_built='', sha256_in_the_stored_table_of_sets='',
                 byte_for_byte_identical='', reproduced=('yes' if g == set(acc(bC)) else 'no'), **common))
# numbers that the plan gives (section 8, eighth entry) against the sets built
MS = M[M.in_C00 == 'True'].copy()
bun = MS.area_by_recorded_admin1 == 'Bunia'
d_ = MS.collection_date
plan = [('genomes of S collected up to 31 August 2026 (upTo0831)', 727, len(members['upTo0831'])), ('genomes of S collected up to 17 August 2026 (upTo0817)', 666, len(members['upTo0817'])),
        ('genomes of S collected from 5 August 2026', 204, int((d_ >= '2026-08-05').sum())), ('genomes of Bunia among them', 110, int(((d_ >= '2026-08-05') & bun).sum())),
        ('genomes of Bunia in the block 5 to 19 August', 47, int(((d_ >= '2026-08-05') & (d_ < '2026-08-20') & bun).sum())),
        ('genomes of Bunia in the block 20 August to 12 September', 63, int(((d_ >= '2026-08-20') & (d_ < '2026-09-13') & bun).sum())),
        ('genomes of Bunia in the reference period 4 June to 4 August', 140, int(((d_ >= '2026-06-04') & (d_ < '2026-08-05') & bun).sum())),
        ('other genomes in the reference period 4 June to 4 August', 322, int(((d_ >= '2026-06-04') & (d_ < '2026-08-05') & ~bun).sum())),
        ('genomes of S', 775, len(MS)), ('most recent collection date of S', '2026-09-12', d_.max()), ('first collection date of S', '2026-05-02', d_.min())]
for what, exp, obs in plan:
    rows.append(dict(kind='sets of the update', test='expected number (plan of the update, section 8, eighth entry, or the specification of the update) against the sets built: ' + what, set='C00',
                     genomes_stored_or_expected=exp, genomes_built=obs, accessions_only_in_the_stored_set='', accessions_only_in_the_built_set='', accessions_that_differ='', same_order='',
                     bytes_stored_or_expected='', bytes_built='', sha256_stored_or_expected='', sha256_built='', sha256_in_the_stored_table_of_sets='', byte_for_byte_identical='',
                     reproduced=('yes' if str(exp) == str(obs) else 'no'), **common))
# the sets beside the stored sets of 2 October 2026 (description)
for r in B.itertuples():
    o = stored.get(f'{r.set}_20261002.fasta')
    nb = open(f'{S}/sets/{r.file}', 'rb').read()
    rows.append(dict(kind='description: set of the update beside the set of 2 October 2026 (no test)', test='members of the set of the update and of the set of 2 October', set=r.set,
                     genomes_stored_or_expected=len(acc(o)), genomes_built=len(acc(nb)), accessions_only_in_the_stored_set=len(set(acc(o)) - set(acc(nb))),
                     accessions_only_in_the_built_set=len(set(acc(nb)) - set(acc(o))), accessions_that_differ='', same_order='', bytes_stored_or_expected=len(o), bytes_built=len(nb),
                     sha256_stored_or_expected=sha(o), sha256_built=sha(nb), sha256_in_the_stored_table_of_sets=T02.loc[r.set, 'sha256'],
                     byte_for_byte_identical=('yes' if o == nb else 'no'), reproduced='not applicable', **common))
R = pd.DataFrame(rows)
R.to_csv(os.path.join(a.out, f'INTERNAL_reproduction_test_sensitivity_{L}.csv'), index=False)
out = R[R.kind != 'description: set of the update beside the set of 2 October 2026 (no test)'].groupby(['kind', 'test'], sort=False).reproduced.agg(lambda x: f"{(x == 'yes').sum()} of {len(x)}")
for (k, t_), v in out.items():
    print(v, '|', t_[:150])
same = R[(R.kind.str.startswith('description')) & (R.byte_for_byte_identical == 'yes')].set.tolist()
print('sets of the update with the same bytes as on 2 October:', same)
