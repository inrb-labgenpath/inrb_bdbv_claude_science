"""spatial_sets_lib_20261002.py - the regional set by the rule of the earlier round.

The statements are those of the code that wrote regional_set_20260924.csv and the three alignments of the
time trees (code kept with the archive of the regional runs of the earlier round; set rule, ADAR masking,
writer of the alignments, area sets and subsamples of prepare_area_sets.py).  Nothing is estimated here.
What is passed in instead of being read from fixed paths: the table of the genome sets, the validation
table, the aligned genomes, the names of the columns that say whether a genome is in the analysed set and
why not, and the three numbers of the rule (threshold of the base rule, lower and upper bound of the
exception).  Every place where a constant of the earlier round stood is marked 'CONSTANT'.
"""
import gzip
import hashlib
import os

import numpy as np
import pandas as pd

L = 18900
GL = 18940
ITURI_MAIN = ['Bunia', 'Rwampara', 'Nizi', 'Mongbwalu', 'Lita', 'Mangala', 'Bambu']
NONIT = ['Nord-Kivu', 'Sud-Kivu']
VAL_COLS = ['geoLocCountry', 'geoLocAdmin1', 'geoLocAdmin2', 'zone', 'zone_rule', 'zone_ref_province',
            'flag_province_conflict', 'flag_zone_not_in_reference', 'flag_zone_interpreted',
            'flag_zone_missing', 'collection_date_precision', 'collection_date_iso8601',
            'sampleCollectionDate']


def read_fasta(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    seqs, name, buf = {}, None, []
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None: seqs[name] = ''.join(buf)
                name = line[1:]; buf = []
            else: buf.append(line.strip())
    if name is not None: seqs[name] = ''.join(buf)
    return seqs


class Snapshot:
    """Aligned genomes of a snapshot: called fraction (positions 1 to 18,900), founder alleles, ADAR-type
    clusters and the masked sequence, as in the earlier code."""

    def __init__(self, aligned_fasta):
        snap = read_fasta(aligned_fasta)
        self.snap = snap
        names = list(snap)
        self.names = names
        b2i = np.full(256, -1, dtype=np.int8)
        for i, c in enumerate('ACGT'): b2i[ord(c)] = i
        M = np.stack([b2i[np.frombuffer(snap[n].upper().encode(), dtype=np.uint8)] for n in names])
        self.M = M
        cf_raw = (M[:, :L] >= 0).mean(1)
        self.cf = pd.Series(cf_raw, index=names)
        pol = np.round(cf_raw, 3) >= 0.95
        counts = np.stack([(M[pol][:, :L] == i).sum(0) for i in range(4)])
        self.founder = counts.argmax(0)
        self.nocov = counts.sum(0) == 0
        self.idx = {n: k for k, n in enumerate(names)}
        self.clusters = {n: self.adar_masks(M[self.idx[n]]) for n in names}

    def adar_masks(self, row, window=300, kmin=3):
        founder, nocov = self.founder, self.nocov
        dpos = np.where((row[:L] >= 0) & (row[:L] != founder) & ~nocov)[0]
        cl = []
        for ref, alt in ((3, 1), (0, 2)):
            pos = sorted(int(p) + 1 for p in dpos if founder[p] == ref and row[p] == alt)
            g = []
            for p in pos:
                if g and p - g[-1][-1] <= window: g[-1].append(p)
                else: g.append([p])
            cl += [c for c in g if len(c) >= kmin]
        return cl

    def masked_seq(self, n):
        s = list(self.snap[n].upper())
        for c in self.clusters[n]:
            for p in c[1:]: s[p - 1] = 'N'
        for p in range(L, GL): s[p] = 'N'
        return ''.join(s)


def area_admin1(r):
    a1 = r['geoLocAdmin1']
    if r['geoLocCountry'] == 'Uganda': return 'Uganda'
    if r['geoLocCountry'] != 'Democratic Republic of the Congo': return None
    if a1 == 'Nord-Kivu': return 'Nord-Kivu'
    if a1 == 'Sud-Kivu': return 'Sud-Kivu'
    if a1 == 'Ituri': return r['zone'] if r['zone'] in ITURI_MAIN else 'other Ituri'
    return None


def area_zoneprov(r):
    if r['flag_province_conflict']:
        zp = r['zone_ref_province']
        if zp in ('Nord-Kivu', 'Sud-Kivu'): return zp
        if zp == 'Ituri': return r['zone'] if r['zone'] in ITURI_MAIN else 'other Ituri'
        return None
    return area_admin1(r)


def build_T(sets, val, S, base_min_cf=0.90, exc_min_cf=0.70, exc_max_cf=0.90):
    """sets: table of the genome sets with the column names of alignment_sets_20260924.csv
    (accessionVersion, collection_date, screen_class, release, shares_sample_identifier,
    restricted_use_agreement_pending, masked_positions_1based, in_primary, reason_not_in_primary);
    val: validation table of the snapshot; S: Snapshot.
    CONSTANT base_min_cf: 0.90 in the earlier code ('base = ... & (cf_exact >= 0.90)'); None = no bound.
    CONSTANT exc_min_cf, exc_max_cf: 0.70 and 0.90 in the earlier code."""
    v2 = val.set_index('accessionVersion')
    s2 = sets.set_index('accessionVersion')
    T = s2.join(v2[VAL_COLS])
    T['date_parsed'] = pd.to_datetime(T['sampleCollectionDate'], format='ISO8601', errors='coerce')
    T['day_precision'] = T['sampleCollectionDate'].astype(str).str.fullmatch(r'\d{4}-\d{2}-\d{2}') & T['date_parsed'].notna()
    T['area_by_admin1'] = T.apply(area_admin1, axis=1)
    T['area_by_zone_province'] = T.apply(area_zoneprov, axis=1)
    T['cf_exact'] = S.cf.reindex(T.index)
    cls_ok = T['screen_class'].isin(['A', 'A*'])
    base = T['in_primary'] & cls_ok & T['day_precision']
    if base_min_cf is not None:
        base = base & (T['cf_exact'] >= base_min_cf)
    other_excl = T['shares_sample_identifier'] | T['restricted_use_agreement_pending']
    exc_pool = (~T['in_primary']) & cls_ok & T['day_precision'] & (T['cf_exact'] >= exc_min_cf) & (T['cf_exact'] < exc_max_cf) \
        & (T['reason_not_in_primary'] == 'called fraction below threshold') & ~other_excl
    for v in ['admin1', 'zone_province']:
        a = T[f'area_by_{v}']
        T[f'exception_{v}'] = exc_pool & a.isin(NONIT)
        T[f'in_regional_{v}'] = (base & a.notna()) | T[f'exception_{v}']
    T['in_regional_union'] = T.in_regional_admin1 | T.in_regional_zone_province
    T['adar_clusters'] = pd.Series({n: len(S.clusters[n]) for n in S.names}, name='adar_clusters').reindex(T.index)
    return T


def reason(r, set_name='the primary set', base_min_cf=0.90, exc_min_cf=0.70):
    if r['in_regional_union']: return ''
    if r['screen_class'] not in ('A', 'A*'): return f"screen class {r['screen_class']} (regional set keeps A and A* only)"
    if not r['day_precision']: return 'collection date not given to the day'
    if r['area_by_admin1'] is None and r['area_by_zone_province'] is None and r['in_primary']: return 'no area: recorded outside the Democratic Republic of the Congo and Uganda'
    if not r['in_primary']:
        rs = r['reason_not_in_primary']
        if rs == 'called fraction below threshold':
            a = r['area_by_admin1']
            if r['cf_exact'] < exc_min_cf: return 'called fraction below 0.70'
            return f"called fraction below 0.90 (area {a}; the 0.70 exception applies to provinces of the DRC other than Ituri)"
        return f'not in {set_name}: {rs}'
    return 'not classified'


def regional_set_table(T, set_name='the primary set', in_set_column='in_primary_set'):
    RS = pd.DataFrame({
        'accessionVersion': T.index, 'collection_date': T['collection_date'].values, 'collection_date_to_the_day': T['day_precision'].values,
        'called_fraction_exact': T['cf_exact'].round(6).values, 'screen_class_rules_v2': T['screen_class'].values, 'release': T['release'].values,
        'country': T['geoLocCountry'].values, 'recorded_admin1': T['geoLocAdmin1'].values, 'health_zone_harmonised': T['zone'].values, 'zone_harmonisation_rule': T['zone_rule'].values,
        'province_of_health_zone': T['zone_ref_province'].values, 'flag_province_conflict': T['flag_province_conflict'].values, in_set_column: T['in_primary'].values,
        'area_by_recorded_admin1': T['area_by_admin1'].values, 'area_by_province_of_health_zone': T['area_by_zone_province'].values,
        'in_regional_set_by_recorded_admin1': T['in_regional_admin1'].values, 'in_regional_set_by_province_of_health_zone': T['in_regional_zone_province'].values,
        'entered_under_completeness_exception_by_recorded_admin1': T['exception_admin1'].values, 'entered_under_completeness_exception_by_province_of_health_zone': T['exception_zone_province'].values,
        'in_time_tree_alignment': T['in_regional_union'].values, 'adar_type_clusters': T['adar_clusters'].values, 'masked_positions_1based': T['masked_positions_1based'].values,
        'reason_not_in_regional_set': [reason(r, set_name) for _, r in T.iterrows()]})
    return RS


def write_fasta(fn, ids, T, S):
    ids = list(ids)
    with open(fn, 'w') as fh:
        for n in ids:
            fh.write(f">{n}|{T.loc[n, 'collection_date']}\n{S.masked_seq(n)}\n")
    return len(ids)


def write_time_tree_alignments(T, S, outdir, suffix):
    """the three alignments of the earlier code; order of the genomes = order of the aligned snapshot"""
    os.makedirs(outdir, exist_ok=True)
    order = [n for n in S.names if n in T.index]
    union_ids = [n for n in order if T.loc[n, 'in_regional_union']]
    base_ids = [n for n in order if T.loc[n, 'in_regional_union'] and not (T.loc[n, 'exception_admin1'] or T.loc[n, 'exception_zone_province'])]
    adm_ids = [n for n in order if T.loc[n, 'in_regional_admin1']]
    out = {}
    for key, ids in [('union', union_ids), ('noexception', base_ids), ('admin1set', adm_ids)]:
        fn = os.path.join(outdir, f'regional_{key}_{suffix}.fasta')
        write_fasta(fn, ids, T, S)
        out[key] = (fn, ids)
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def fasta_ids(path):
    return [l[1:].strip() for l in open(path) if l.startswith('>')]


def prepare_area_sets(T, union_fasta, outdir_areas, outdir_ss, subsample_seed=20260924, min_n=25,
                      sizes=(20, 30, 45, 70, 100, 130, 180), nrep=8, fasta_prefix='regional'):
    """prepare_area_sets.py of the earlier round with the paths passed in.
    CONSTANT subsample_seed: 20260924 is fixed in the earlier code (generator of the subsamples)."""
    seqs = {}
    name = None
    for line in open(union_fasta):
        if line.startswith('>'):
            name = line[1:].strip()
        else:
            seqs[name.split('|')[0]] = (name, line.strip())
    os.makedirs(outdir_areas, exist_ok=True); os.makedirs(outdir_ss, exist_ok=True)
    order = list(seqs)

    def write(fn, ids):
        with open(fn, 'w') as fh:
            for a in ids:
                fh.write('>%s\n%s\n' % seqs[a])
    rows = []
    sets = {}
    for assign in ['admin1', 'zone_province']:
        memb = T['in_regional_' + assign]
        area = T['area_by_' + assign]
        for ar, g in T[memb].groupby(area[memb]):
            ids = [a for a in order if a in set(g.index)]
            sets.setdefault(ar, {})[assign] = ids
        hub = [a for a in order if memb.get(a, False) and area[a] in ('Bunia', 'Rwampara')]
        sets.setdefault('Bunia+Rwampara', {})[assign] = hub
    for ar, d in sets.items():
        same = d.get('admin1') == d.get('zone_province')
        for assign, ids in d.items():
            if same and assign == 'zone_province':
                continue
            tag = ar.replace(' ', '_').replace('+', '_') + ('' if same else '_' + assign)
            n = len(ids)
            dd = pd.to_datetime(T.loc[ids, 'collection_date'], format='ISO8601')
            fit = n >= min_n
            fn = os.path.join(outdir_areas, f'{tag}.fasta')
            if fit:
                write(fn, ids)
            rows.append(dict(tag=tag, area=ar, assignment='both (identical genome set)' if same else assign, n=n, first=str(dd.min().date()), last=str(dd.max().date()),
                             fitted=fit, n_exception_genomes=int((T.loc[ids, 'exception_admin1'] | T.loc[ids, 'exception_zone_province']).sum()),
                             fasta=f'{fasta_prefix}/areas/{tag}.fasta' if fit else '', sha256=hashlib.sha256(open(fn, 'rb').read()).hexdigest() if fit else ''))
    A = pd.DataFrame(rows).sort_values('n', ascending=False)
    A.to_csv(os.path.join(outdir_areas, 'area_sets.csv'), index=False)
    best = 'Bunia'
    src = sets[best]['admin1']
    assert sets[best]['admin1'] == sets[best]['zone_province']
    rng = np.random.default_rng(subsample_seed)
    jobs = []
    for rep in range(nrep):
        for n_s in sizes:
            ids = sorted(rng.choice(src, n_s, replace=False).tolist(), key=order.index)
            tag = f'bunia_{n_s}_{rep}'
            write(os.path.join(outdir_ss, f'{tag}.fasta'), ids)
            dd = pd.to_datetime(T.loc[ids, 'collection_date'], format='ISO8601')
            jobs.append(dict(tag=tag, n=n_s, rep=rep, first=str(dd.min().date()), last=str(dd.max().date()), accessions=';'.join(ids)))
    J = pd.DataFrame(jobs)
    J.to_csv(os.path.join(outdir_ss, 'subsamples.csv'), index=False)
    return A, J, sets
