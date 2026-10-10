#!/usr/bin/env python3
"""analyses_config.py - the analyses of the settings C and D and of the Delphy analysis with exponential growth (round of October 2026).

Names and command templates only.  Seeds are read from the schedule (phase5_schedule_20261003.csv); the commands of
the counterparts are read from the run manifests of the earlier round.  No estimate and no number of genomes is written here.
"""
import csv

SUFFIX = "20261002"
TRACK = "beast_full"

# BEAST X: launch options of the commands of the run manifest of BEAST X of 24 September; the number of threads is
# one thread per chain.
BEAST_TEMPLATE = ("beast -threads {threads} -beagle_CPU -beagle_SSE -beagle_threads {threads} -seed {seed} "
                  "-save_every 500000 -save_state {chain}.state -overwrite ../../../xml/{stem}.xml")
# environment that the launch script of BEAST X of the earlier round sets (work/beast/run_beast_chains.sh)
BEAST_PRELUDE = "export JAVA_TOOL_OPTIONS='-XX:ParallelGCThreads=2 -XX:ConcGCThreads=1'"
# Delphy: the command of the counterpart 'exp_c8k_primary' of the earlier round with the alignment, the seed and the output names replaced
DELPHY_EXP_TEMPLATE = ("programs/delphy_1.4.1/bin/delphy --v0-in-fasta {alignment} --v0-steps 1500000000 --v0-seed {seed} "
                       "--v0-threads 2 --v0-pop-model exponential --v0-log-every 300000 --v0-tree-every 30000000 "
                       "--v0-out-log-file runs/{analysis}/{chain}.log --v0-out-trees-file runs/{analysis}/{chain}.trees "
                       "--v0-target-coal-prior-cells 8000")
# the primary command (README_round_20261002.md, section 4), for the comparison that rule 7 asks for
DELPHY_PRIMARY_TEMPLATE = ("programs/delphy_1.4.1/bin/delphy --v0-in-fasta {alignment} --v0-steps 1500000000 --v0-seed {seed} "
                           "--v0-threads 2 --v0-pop-model skygrid --v0-skygrid-num-parameters 20 --v0-skygrid-cutoff 0.8 "
                           "--v0-log-every 300000 --v0-tree-every 30000000 "
                           "--v0-out-log-file runs/{analysis}/{chain}.log --v0-out-trees-file runs/{analysis}/{chain}.trees "
                           "--v0-target-coal-prior-cells 8000")

ANALYSES = [
    dict(analysis="beast_C", program="BEAST X 10.5.0", kind="beast", genome_set="C00",
         setting="C: Skygrid, 32 parameters, cut-off 0.613699 years, precision estimated (Gamma(0.001, 1000) prior), GTR, strict clock",
         alignment="aln/C00_20261002.fasta", stem="beastC_weekly32_C00_20261002", threads=1, steps=40000000,
         log_every=20000, tree_every=2000000,
         counterpart_manifest="beast", counterpart_analysis="C_inrblike_primary", counterpart_stem="beastC_inrblike_primary_20260924"),
    dict(analysis="beast_D", program="BEAST X 10.5.0", kind="beast", genome_set="C00",
         setting="D: constant exponential growth, HKY, strict clock",
         alignment="aln/C00_20261002.fasta", stem="beastD_exp_C00_20261002", threads=1, steps=40000000,
         log_every=20000, tree_every=2000000,
         counterpart_manifest="beast", counterpart_analysis="D_exp_primary", counterpart_stem="beastD_exp_primary_20260924"),
    dict(analysis="delphy_exponential", program="Delphy 1.4.1", kind="delphy", genome_set="C00",
         setting="constant exponential growth, 8,000 cells",
         alignment="aln/C00_20261002.fasta", stem="", threads=2, steps=1500000000, log_every=300000, tree_every=30000000,
         counterpart_manifest="delphy", counterpart_analysis="exp_c8k_primary", counterpart_stem=""),
    dict(analysis="beast_C_retained1046", program="BEAST X 10.5.0", kind="beast", genome_set="retained1046",
         setting="C: Skygrid, 32 parameters, cut-off 0.613699 years, precision estimated (Gamma(0.001, 1000) prior), GTR, strict clock",
         alignment="aln/retained1046_20261002.fasta", stem="beastC_weekly32_retained1046_20261002", threads=1, steps=40000000,
         log_every=20000, tree_every=2000000,
         counterpart_manifest="beast", counterpart_analysis="C_inrblike_inrbRetained",
         counterpart_stem="beastC_inrblike_inrbRetained_20260924"),
]


def seeds_of(schedule_csv, analysis):
    with open(schedule_csv, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["analysis"] == analysis]
    rows.sort(key=lambda r: int(r["chain"]))
    n = {int(r["chains_of_the_analysis"]) for r in rows}
    assert rows and n == {len(rows)}, (analysis, len(rows), n)
    return [int(r["seed"]) for r in rows]


def chain_name(analysis, seed):
    return f"{analysis}_s{seed}"


def command_of(a, seed):
    chain = chain_name(a["analysis"], seed)
    if a["kind"] == "beast":
        return BEAST_TEMPLATE.format(threads=a["threads"], seed=seed, chain=chain, stem=a["stem"])
    return DELPHY_EXP_TEMPLATE.format(alignment=a["alignment"], seed=seed, analysis=a["analysis"], chain=chain)
