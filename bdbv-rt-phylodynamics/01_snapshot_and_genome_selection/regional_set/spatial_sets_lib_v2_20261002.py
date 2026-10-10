"""spatial_sets_lib_v2_20261002.py - second version of spatial_sets_lib_20261002.py (that file is kept unchanged).

What differs: the text of the reason why a genome is not in the regional set
takes the two thresholds from the parameters of the rule instead of holding 0.90 and 0.70 as text
(CONSTANT of the earlier round: 0.90 and 0.70 in the two sentences of reason()), and it has one branch more, for a
genome that is in the analysed set and lies below the threshold of the base rule (it cannot occur where the analysed
set has that threshold itself, as in the earlier round).  With the parameters 0.90 and 0.70 the table of the
earlier round is written byte for byte (tested).  The rule of the set itself is that of spatial_sets_lib_20261002.py.
"""
import pandas as pd

from spatial_sets_lib_20261002 import *          # noqa: F401,F403
import spatial_sets_lib_20261002 as _v1


def reason_v2(r, set_name='the primary set', base_min_cf=0.90, exc_min_cf=0.70):
    if r['in_regional_union']: return ''
    if r['screen_class'] not in ('A', 'A*'): return f"screen class {r['screen_class']} (regional set keeps A and A* only)"
    if not r['day_precision']: return 'collection date not given to the day'
    if r['area_by_admin1'] is None and r['area_by_zone_province'] is None and r['in_primary']: return 'no area: recorded outside the Democratic Republic of the Congo and Uganda'
    below = (lambda a: f"called fraction below {base_min_cf:.2f} (area {a}; the {exc_min_cf:.2f} exception applies to provinces of the DRC other than Ituri)")
    if not r['in_primary']:
        rs = r['reason_not_in_primary']
        if rs == 'called fraction below threshold':
            a = r['area_by_admin1']
            if r['cf_exact'] < exc_min_cf: return f'called fraction below {exc_min_cf:.2f}'
            return below(a)
        return f'not in {set_name}: {rs}'
    if base_min_cf is not None and r['cf_exact'] < base_min_cf:
        return below(r['area_by_admin1'])
    return 'not classified'


def regional_set_table_v2(T, set_name='the primary set', in_set_column='in_primary_set', base_min_cf=0.90, exc_min_cf=0.70):
    RS = _v1.regional_set_table(T, set_name, in_set_column)
    RS['reason_not_in_regional_set'] = [reason_v2(r, set_name, base_min_cf, exc_min_cf) for _, r in T.iterrows()]
    return RS
