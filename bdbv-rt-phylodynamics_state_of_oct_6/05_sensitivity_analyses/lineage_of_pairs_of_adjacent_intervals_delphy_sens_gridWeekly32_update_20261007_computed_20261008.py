import importlib
import sys
import os
import tarfile

import numpy as np
import pandas as pd


class Stop(Exception):
    pass


def compute(paths, pipeline_dir, out_path, expect):
    if pipeline_dir not in sys.path:
        sys.path.insert(0, pipeline_dir)
    L = importlib.import_module("rt_lib")
    GEN = L.GenTime(15.3, 9.3)

    def pairs(pp, g, tau):
        K = g.K
        seeds = list(dict.fromkeys(pp["seed"]))
        cols = ["ln_Ne_tau_%d" % i for i in range(1, K + 1)]
        lp = [pp.loc[pp["seed"] == s, cols].values.astype(float) for s in seeds]
        rows = []
        for i in range(1, K):
            c = L.pair_contrast(g, i)
            rs = [x @ c for x in lp]
            allr = np.concatenate(rs)
            Rr = GEN.r2R(allr)
            lo, hi = L.hpd(Rr)
            s = L.summarise(rs, "pair %d" % i)
            sd = float(np.std(allr)) / L.YEAR_DAYS
            psd = L.prior_sd_contrast(c, tau) / L.YEAR_DAYS
            rows.append(dict(
                pair="%d to %d" % (i + 1, i),
                newer_interval=i,
                older_interval=i + 1,
                date_on_which_the_two_intervals_meet=str(L.dec2date(g.t_anchor - i * g.D)),
                mid_point_of_the_older_interval=str(L.dec2date(g.mids[i])) if i + 1 < K else "",
                mid_point_of_the_newer_interval=str(L.dec2date(g.mids[i - 1])),
                R_median=float(np.median(Rr)),
                R_hpd95_lower=lo,
                R_hpd95_upper=hi,
                probability_that_R_exceeds_1=float((allr > 0).mean()),
                growth_rate_per_day_median=float(np.median(allr)) / L.YEAR_DAYS,
                posterior_sd_of_the_growth_rate_per_day=sd,
                prior_sd_of_the_growth_rate_per_day=float(psd),
                sd_ratio=sd / psd,
                sd_ratio_of_0_90_or_more="yes" if sd / psd >= 0.90 else "no",
                states_with_R_set_to_0=GEN.n_undefined(allr),
                ess_pooled=float(s["ess_pooled"]),
                rhat_split=float(s["rhat_split"]),
                meets_convergence_criterion="yes" if (s["ess_pooled"] >= 200 and s["rhat_split"] <= 1.05) else "no",
                chains_retained=len(seeds),
                states_used=int(len(allr))
            ))
        return pd.DataFrame(rows)

    facts = {}
    names = [l[1:].strip() for l in open(paths["alignment"], encoding="utf-8") if l.startswith(">")]
    dS = pd.to_datetime(pd.Series([n.split("|")[-1] for n in names]))
    last = str(dS.max().date())
    if len(names) != expect["genomes"] or last != expect["last"]:
        raise Stop("the alignment holds %d genomes, most recent %s" % (len(names), last))
    facts.update(genomes=len(names), most_recent_collection_date=last)

    pp0 = pd.read_csv(paths["posterior_primary"], low_memory=False)
    K0, cut0 = expect["primary_grid"]
    g0 = L.Grid.from_cutoff(K0, cut0, last)
    T0 = pairs(pp0, g0, expect["primary_precision"])
    P = pd.read_csv(paths["primary_pairs"], dtype=str, keep_default_na=False)
    if len(T0) != len(P):
        raise Stop("control 1: %d pairs computed, %d in the stored table" % (len(T0), len(P)))
    numeric = ("R_median", "R_hpd95_lower", "R_hpd95_upper", "probability_that_R_exceeds_1", "growth_rate_per_day_median",
               "posterior_sd_of_the_growth_rate_per_day", "prior_sd_of_the_growth_rate_per_day", "sd_ratio", "ess_pooled",
               "rhat_split")
    dev = {c: float(np.max(np.abs(T0[c].values.astype(float) - P[c].astype(float).values))) for c in numeric}
    texts = ("pair", "date_on_which_the_two_intervals_meet", "mid_point_of_the_older_interval",
             "mid_point_of_the_newer_interval", "sd_ratio_of_0_90_or_more", "meets_convergence_criterion")
    same = {c: bool(list(T0[c].astype(str)) == list(P[c])) for c in texts}
    counts = {c: bool(list(T0[c].astype(int)) == list(P[c].astype(int))) for c in ("states_used", "chains_retained",
                                                                                 "states_with_R_set_to_0")}
    facts["control_1"] = dict(pairs=len(T0), largest_differences=dev, text_columns_equal=same, count_columns_equal=counts,
                              states=int(len(pp0)))
    if max(dev.values()) >= 1e-9 or not all(same.values()) or not all(counts.values()):
        raise Stop("control 1 is not met: %s %s %s" % (dev, same, counts))

    Ch = pd.read_csv(paths["chains_weekly"], dtype=str, keep_default_na=False)
    if set(Ch.precision_is_estimated) != {"no"} or len(set(Ch.smallest_precision) | set(Ch.largest_precision)) != 1:
        raise Stop("the chains table of the weekly analysis does not show one fixed precision")
    tau_w = float(Ch.smallest_precision.iloc[0])
    facts["precision_of_the_weekly_analysis"] = Ch.smallest_precision.iloc[0]

    Iw = pd.read_csv(paths["intervals_weekly"], dtype=str, keep_default_na=False)
    Kw, cutw = expect["weekly_grid"]
    if len(Iw) != Kw:
        raise Stop("control 2: the intervals table holds %d rows" % len(Iw))
    m0, m1 = float(Iw.mid_point_decimal[0]), float(Iw.mid_point_decimal[1])
    Dw = m0 - m1
    anchor = m0 + Dw / 2
    gw = L.Grid.from_cutoff(Kw, cutw, last)
    tb = gw.table()
    eq = bool((tb.newer_edge.values == Iw.newer_edge.values).all()
              and (tb.mid_point.values[:Kw - 1] == Iw.mid_point.values[:Kw - 1]).all()
              and (tb.older_edge.values[:Kw - 1] == Iw.older_edge.values[:Kw - 1]).all())
    facts["control_2"] = dict(difference_of_the_interval_length=abs(gw.D - Dw), difference_of_the_anchor=abs(gw.t_anchor - anchor),
                              edges_and_mid_points_equal=eq, interval_days=gw.D * L.YEAR_DAYS)
    if abs(gw.D - Dw) >= 1e-9 or abs(gw.t_anchor - anchor) >= 1e-9 or not eq:
        raise Stop("control 2 is not met: %s" % facts["control_2"])

    ppw = pd.read_csv(paths["posterior_weekly"], low_memory=False)
    med = np.median(ppw[["ln_Ne_tau_%d" % i for i in range(1, Kw + 1)]].values, axis=0)
    dmed = float(np.max(np.abs(med - Iw.ln_Ne_tau_median.astype(float).values)))
    seedsw = list(dict.fromkeys(ppw["seed"]))
    dess = []
    drh = []
    for i in (1, 5, 14, 17, 25):
        s = L.summarise([ppw.loc[ppw["seed"] == sd_, "ln_Ne_tau_%d" % i].values.astype(float) for sd_ in seedsw], "x")
        dess.append(abs(s["ess_pooled"] - float(Iw.ess_pooled[i - 1])))
        drh.append(abs(s["rhat_split"] - float(Iw.rhat_split[i - 1])))
    facts["control_3"] = dict(largest_difference_of_the_medians=dmed, largest_difference_of_ess=max(dess),
                              largest_difference_of_rhat=max(drh), states=int(len(ppw)), chains=len(seedsw),
                              states_as_the_table=bool({str(len(ppw))} == set(Iw.states_used)),
                              chains_as_the_table=bool({str(len(seedsw))} == set(Iw.chains_retained)))
    c3 = facts["control_3"]
    if dmed >= 1e-9 or max(dess) >= 1e-6 or max(drh) >= 1e-9 or not c3["states_as_the_table"] or not c3["chains_as_the_table"]:
        raise Stop("control 3 is not met: %s" % c3)

    Tw = pairs(ppw, gw, tau_w)
    Tw.insert(0, "analysis", expect["analysis"])
    Tw.insert(1, "program", expect["program"])
    Tw.insert(2, "setting", expect["setting"])
    Tw.insert(3, "genome_set", expect["genome_set"])
    Tw.insert(4, "genomes", len(names))
    Tw.insert(5, "most_recent_collection_date", last)
    Tw["interval_days"] = gw.D * L.YEAR_DAYS
    Tw["skygrid_precision_fixed_at"] = tau_w
    form = list(pd.read_csv(paths["form"], dtype=str, keep_default_na=False, nrows=1).columns)
    if list(Tw.columns) != form or len(Tw) != Kw - 1:
        raise Stop("the table has %d rows and columns %s; expected %d rows and %s" % (len(Tw), list(Tw.columns), Kw - 1, form))
    Tw.to_csv(out_path, index=False)
    back = pd.read_csv(out_path, dtype=str, keep_default_na=False)
    if len(back) != len(Tw) or list(back.columns) != form:
        raise Stop("the table read back is not the table written")
    facts["table"] = dict(rows=len(Tw), columns=len(form), open_symbols=int((Tw.sd_ratio_of_0_90_or_more == "yes").sum()),
                          not_converged=int((Tw.meets_convergence_criterion == "no").sum()),
                          smallest_sd_ratio=float(Tw.sd_ratio.min()), largest_sd_ratio=float(Tw.sd_ratio.max()),
                          largest_upper_bound_of_R=float(Tw.R_hpd95_upper.max()))
    return facts


# Extract pipeline
os.makedirs('work/zoom/tree', exist_ok=True)
with tarfile.open('inputs/pipeline_round_20261002.tar.gz') as tf:
    tf.extractall('work/zoom/tree')

pdir = 'work/zoom/tree/pipeline_round_20261002/work/pipeline_ext'

# Load without_chains result for precision
import json
W = json.load(open("inputs/without_chains_update_result.json"))

# Load form columns from the first computed table
form_path = 'inputs/INTERNAL_pairs_of_adjacent_intervals_delphy_sens_gridWeekly32_computed_20261007.csv'

# Load chains to get program/setting info
Ch_path = 'inputs/chains_delphy_sens_gridWeekly32_20261007.csv'
Ch = pd.read_csv(Ch_path, dtype=str, keep_default_na=False)

form1 = pd.read_csv(form_path, dtype=str, keep_default_na=False, nrows=1).iloc[0]

paths = dict(
    alignment='inputs/C00_20261007.fasta',
    posterior_primary='inputs/posterior_delphy_fixed_C00_20261007.csv.gz',
    primary_pairs='inputs/primary_pairs_sensitivity_20261007.csv',
    chains_weekly=Ch_path,
    intervals_weekly='inputs/intervals_delphy_sens_gridWeekly32_20261007.csv',
    posterior_weekly='inputs/posterior_delphy_sens_gridWeekly32_20261007.csv.gz',
    form=form_path,
)

expect = dict(
    genomes=775,
    last="2026-09-12",
    primary_grid=(20, 0.8),
    primary_precision=float(W["facts"]["precision"]),
    weekly_grid=(32, 0.613699),
    analysis="delphy_sens_gridWeekly32",
    program=form1["program"],
    setting=form1["setting"],
    genome_set="C00",
)

out = "pairs_of_adjacent_intervals_delphy_sens_gridWeekly32_update_20261007_computed_20261008.csv"
facts = compute(paths, pdir, out, expect)
print("done:", facts["table"])