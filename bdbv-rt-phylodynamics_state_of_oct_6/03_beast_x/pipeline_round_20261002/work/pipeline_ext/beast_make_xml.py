#!/usr/bin/env python3
"""beast_make_xml.py - BEAST X (v10.5) XML for the dated, masked alignment.

Recovered from make_xml.py of beast_reanalysis_482.tar.gz and from the fixed-precision template
(skygridfix1.xml) that was re-used on 12 Sep 2026, and parameterised.

Tip-date text (--tip-date-text): 'fixed6' writes six decimals (2026.397260), as the XML builder of
12 Sep did; 'shortest' writes the same number without trailing zeros (2026.39726), as the builder of
the bundle did.  Comparison with the bundle (inputs of the bundle, 482 genomes; regression report,
section 7.1): with 'fixed6' the generated exp1.xml, skygrid1.xml and skygridfix1.xml differ from the
bundle's files in 56 of 482 taxon lines, only in the date text.  The 12 Sep XML files themselves
were not saved, so identity with them cannot be checked directly.

Models (substitution and clock model as in Delphy: HKY, strict clock, tip-dated)
  exp                exponential-growth coalescent; 1/x prior on the population size, Laplace prior
                     (mean 0, scale 100 per year) on the growth rate
  skygrid-fixed      Skygrid with the GMRF precision FIXED at --precision (no prior, no operator on
                     the precision); the log population sizes are moved by two random-walk operators
                     (window 1.0 and 0.3, weight 15 each).  12 Sep value: 4.06, i.e. Delphy's fixed
                     precision 4.06335 (double-half time 30 d, 20 parameters, cutoff 0.8 yr) rounded.
  skygrid-estimated  Skygrid with ESTIMATED precision: Gamma(shape 0.001, scale 1000) prior, GMRF
                     block-update operator (scale 2.0, weight 2) and scale operator on the precision
                     (BEAST default parameterisation)

12 Sep run settings: chain length 12,000,000; parameters logged every 5,000, trees every 60,000,
screen every 50,000 states;  beast -threads 4 -seed <S> -overwrite <file>.xml  run inside the output
directory; seeds 111, 112 (exponential) and 113, 114 (fixed-precision Skygrid); 30 % burn-in;
Summary tree on 12 Sep: treeannotator -burnin 3600000 -heights median <trees> <out.nexus>.
TreeAnnotator v10.5.0 builds a HIPSTR tree unless -type is given, so the tree saved on 12 Sep under
the name *.mcc.nexus is a HIPSTR tree.  A maximum-clade-credibility tree needs  -type mcc.
(The bundle runs of 1 Sep were started with  -seed <S> -threads 4 -beagle_CPU -beagle_SSE -overwrite,
seeds 101, 102 for the exponential and estimated-precision models and 401, 402 for fixed precision.)

Added for the 24 Sep 2026 round (defaults unchanged, so the 12 Sep XML files are still reproduced)
  --subst-model gtr       GTR instead of HKY: six relative rates (BEAST X parameterisation, Dirichlet(1)
                          prior on rates summing to 6, delta-exchange operator), frequencies as for HKY
  --num-parameters / --cutoff   any Skygrid grid (BEAUti convention: number of parameters = number of
                          transition points + 1; cutoff = time of the last transition point, years)
  --root-height-min X     sensitivity analysis only: uniform prior with lower bound X (years before the
                          last tip) on the root height, i.e. an upper bound on the root date
  --logpop-lower X        sensitivity analysis only: lower bound X on every log population size
  --screen-every N        screen log frequency (default 10 x log-every)
Checkpointing is a command-line matter (beast -save_every N -save_state FILE, -load_state FILE).

Usage
  beast_make_xml.py --fasta cleaned_<suffix>.fasta --model skygrid-fixed --tag 1 \
                    --out-suffix 20260924 --outdir beastruns [--precision 4.06] [--chain-length N]
writes <outdir>/<name>.xml with name = <model-name><tag>_<suffix> (or --name), and
<outdir>/beast_xml_parameters_<name>.json
"""
import argparse
import gzip
import hashlib
import json
import os
import sys

import pandas as pd

MODEL_NAMES = {"exp": "exp", "skygrid-fixed": "skygridfix", "skygrid-estimated": "skygrid"}


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def decimal_year(ts):
    ts = pd.Timestamp(ts)
    start = pd.Timestamp(ts.year, 1, 1)
    end = pd.Timestamp(ts.year + 1, 1, 1)
    return ts.year + (ts - start).days / (end - start).days


def fmt_precision(x):
    s = repr(float(x))
    return s


def xml(model, taxa_labels, decimals, SEQ, chain=12_000_000, log_every=5_000, tree_every=60_000,
        prefix="run", precision=4.06, num_parameters=20, cutoff=0.8, subst_model="hky",
        root_height_min=None, logpop_lower=None, screen_every=None):
    taxa = "\n".join(
        f'    <taxon id="{t}"><date value="{d}" direction="forwards" units="years"/></taxon>'
        for t, d in zip(taxa_labels, decimals))
    seqs = "\n".join(f'    <sequence><taxon idref="{t}"/>{SEQ[t]}</sequence>' for t in taxa_labels)
    grid_points = f"{float(num_parameters - 1)}"

    if model == "exp":
        demo = """
  <exponentialGrowth id="exponential" units="years">
    <populationSize><parameter id="exponential.popSize" value="0.1" lower="0.0"/></populationSize>
    <growthRate><parameter id="exponential.growthRate" value="10.0"/></growthRate>
  </exponentialGrowth>
  <coalescentLikelihood id="coalescent">
    <model><exponentialGrowth idref="exponential"/></model>
    <populationTree><treeModel idref="treeModel"/></populationTree>
  </coalescentLikelihood>"""
        demo_ops = """
    <scaleOperator scaleFactor="0.75" weight="3"><parameter idref="exponential.popSize"/></scaleOperator>
    <randomWalkOperator windowSize="1.0" weight="3"><parameter idref="exponential.growthRate"/></randomWalkOperator>"""
        demo_priors = """
      <oneOnXPrior><parameter idref="exponential.popSize"/></oneOnXPrior>
      <laplacePrior mean="0.0" scale="100.0"><parameter idref="exponential.growthRate"/></laplacePrior>
      <coalescentLikelihood idref="coalescent"/>"""
        demo_log = """
    <parameter idref="exponential.popSize"/>
    <parameter idref="exponential.growthRate"/>
    <coalescentLikelihood idref="coalescent"/>"""
    else:
        start_precision = "0.1" if model == "skygrid-estimated" else fmt_precision(precision)
        demo = f"""
  <gmrfSkyGridLikelihood id="skygrid">
    <populationSizes><parameter id="skygrid.logPopSize" dimension="{num_parameters}" value="1.0"/></populationSizes>
    <precisionParameter><parameter id="skygrid.precision" value="{start_precision}" lower="0.0"/></precisionParameter>
    <numGridPoints><parameter id="skygrid.numGridPoints" value="{grid_points}"/></numGridPoints>
    <cutOff><parameter id="skygrid.cutOff" value="{cutoff}"/></cutOff>
    <populationTree><treeModel idref="treeModel"/></populationTree>
  </gmrfSkyGridLikelihood>"""
        if model == "skygrid-estimated":
            demo_ops = """
    <gmrfGridBlockUpdateOperator scaleFactor="2.0" weight="2">
      <gmrfSkyGridLikelihood idref="skygrid"/>
    </gmrfGridBlockUpdateOperator>
    <scaleOperator scaleFactor="0.75" weight="1"><parameter idref="skygrid.precision"/></scaleOperator>"""
            demo_priors = """
      <gammaPrior shape="0.001" scale="1000.0" offset="0.0"><parameter idref="skygrid.precision"/></gammaPrior>
      <gmrfSkyGridLikelihood idref="skygrid"/>"""
        else:
            demo_ops = """
    <randomWalkOperator windowSize="1.0" weight="15"><parameter idref="skygrid.logPopSize"/></randomWalkOperator>
    <randomWalkOperator windowSize="0.3" weight="15"><parameter idref="skygrid.logPopSize"/></randomWalkOperator>"""
            demo_priors = """
      
      <gmrfSkyGridLikelihood idref="skygrid"/>"""
        demo_log = """
    <parameter idref="skygrid.precision"/>
    <parameter idref="skygrid.logPopSize"/>
    <gmrfSkyGridLikelihood idref="skygrid"/>"""

    if subst_model == "gtr":
        subst_block = """  <gtrModel id="gtr">
    <frequencies>
      <frequencyModel dataType="nucleotide">
        <frequencies><parameter id="frequencies" value="0.25 0.25 0.25 0.25"/></frequencies>
      </frequencyModel>
    </frequencies>
    <rates><parameter id="gtr.rates" dimension="6" value="1.0" lower="0.0"/></rates>
  </gtrModel>
  <siteModel id="siteModel">
    <substitutionModel><gtrModel idref="gtr"/></substitutionModel>
  </siteModel>"""
        subst_ops = '    <deltaExchange delta="0.01" weight="1"><parameter idref="gtr.rates"/></deltaExchange>'
        subst_priors = '        <dirichletPrior alpha="1.0" sumsTo="6.0"><parameter idref="gtr.rates"/></dirichletPrior>'
        subst_log = '      <parameter idref="gtr.rates"/>'
    else:
        subst_block = """  <HKYModel id="hky">
    <frequencies>
      <frequencyModel dataType="nucleotide">
        <frequencies><parameter id="frequencies" value="0.25 0.25 0.25 0.25"/></frequencies>
      </frequencyModel>
    </frequencies>
    <kappa><parameter id="kappa" value="2.0" lower="0.0"/></kappa>
  </HKYModel>
  <siteModel id="siteModel">
    <substitutionModel><HKYModel idref="hky"/></substitutionModel>
  </siteModel>"""
        subst_ops = '    <scaleOperator scaleFactor="0.75" weight="1"><parameter idref="kappa"/></scaleOperator>'
        subst_priors = '        <logNormalPrior mu="1.0" sigma="1.25" offset="0.0"><parameter idref="kappa"/></logNormalPrior>'
        subst_log = '      <parameter idref="kappa"/>'
    extra_priors = ""
    if root_height_min is not None:
        extra_priors += (f'\n        <uniformPrior lower="{float(root_height_min)}" upper="1.0E100">'
                         '<parameter idref="treeModel.rootHeight"/></uniformPrior>')
    if logpop_lower is not None and model != "exp":
        extra_priors += (f'\n        <uniformPrior lower="{float(logpop_lower)}" upper="1.0E100">'
                         '<parameter idref="skygrid.logPopSize"/></uniformPrior>')

    return f"""<?xml version="1.0" standalone="yes"?>
<beast version="1.10.4">
  <taxa id="taxa">
{taxa}
  </taxa>
  <alignment id="alignment" dataType="nucleotide">
{seqs}
  </alignment>
  <patterns id="patterns" from="1" strip="false"><alignment idref="alignment"/></patterns>

  <constantSize id="initialDemo" units="years">
    <populationSize><parameter id="initialDemo.popSize" value="0.1"/></populationSize>
  </constantSize>
  <coalescentSimulator id="startingTree">
    <taxa idref="taxa"/>
    <constantSize idref="initialDemo"/>
  </coalescentSimulator>
  <treeModel id="treeModel">
    <coalescentTree idref="startingTree"/>
    <rootHeight><parameter id="treeModel.rootHeight"/></rootHeight>
    <nodeHeights internalNodes="true"><parameter id="treeModel.internalNodeHeights"/></nodeHeights>
    <nodeHeights internalNodes="true" rootNode="true"><parameter id="treeModel.allInternalNodeHeights"/></nodeHeights>
  </treeModel>
  <treeLengthStatistic id="treeLength"><treeModel idref="treeModel"/></treeLengthStatistic>
{demo}

  <strictClockBranchRates id="branchRates">
    <rate><parameter id="clock.rate" value="1.0E-3" lower="0.0"/></rate>
  </strictClockBranchRates>
  <rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">
    <treeModel idref="treeModel"/><strictClockBranchRates idref="branchRates"/>
  </rateStatistic>

{subst_block}

  <treeDataLikelihood id="treeLikelihood" useAmbiguities="false">
    <partition>
      <patterns idref="patterns"/>
      <siteModel idref="siteModel"/>
    </partition>
    <treeModel idref="treeModel"/>
    <strictClockBranchRates idref="branchRates"/>
  </treeDataLikelihood>

  <operators id="operators" optimizationSchedule="default">
{subst_ops}
    <deltaExchange delta="0.01" weight="1"><parameter idref="frequencies"/></deltaExchange>
    <scaleOperator scaleFactor="0.75" weight="3"><parameter idref="clock.rate"/></scaleOperator>
    <subtreeSlide size="1.0E-3" gaussian="true" weight="30"><treeModel idref="treeModel"/></subtreeSlide>
    <narrowExchange weight="30"><treeModel idref="treeModel"/></narrowExchange>
    <wideExchange weight="3"><treeModel idref="treeModel"/></wideExchange>
    <wilsonBalding weight="3"><treeModel idref="treeModel"/></wilsonBalding>
    <scaleOperator scaleFactor="0.75" weight="3"><parameter idref="treeModel.rootHeight"/></scaleOperator>
    <uniformOperator weight="30"><parameter idref="treeModel.internalNodeHeights"/></uniformOperator>
    <upDownOperator scaleFactor="0.75" weight="3">
      <up><parameter idref="clock.rate"/></up>
      <down><parameter idref="treeModel.allInternalNodeHeights"/></down>
    </upDownOperator>{demo_ops}
  </operators>

  <mcmc id="mcmc" chainLength="{chain}" autoOptimize="true" operatorAnalysis="{prefix}.ops">
    <joint id="joint">
      <prior id="prior">
{subst_priors}
        <uniformPrior lower="0.0" upper="1.0"><parameter idref="frequencies"/></uniformPrior>
        <ctmcScalePrior><ctmcScale><parameter idref="clock.rate"/></ctmcScale><treeModel idref="treeModel"/></ctmcScalePrior>{demo_priors}{extra_priors}
      </prior>
      <likelihood id="likelihood">
        <treeDataLikelihood idref="treeLikelihood"/>
      </likelihood>
    </joint>
    <operators idref="operators"/>
    <log id="screenLog" logEvery="{screen_every or log_every * 10}">
      <column label="Joint" dp="4" width="12"><joint idref="joint"/></column>
      <column label="rootHeight" sf="6" width="12"><parameter idref="treeModel.rootHeight"/></column>
      <column label="clock.rate" sf="6" width="12"><parameter idref="clock.rate"/></column>
    </log>
    <log id="fileLog" logEvery="{log_every}" fileName="{prefix}.log" overwrite="true">
      <joint idref="joint"/><prior idref="prior"/><likelihood idref="likelihood"/>
      <parameter idref="treeModel.rootHeight"/>
      <treeLengthStatistic idref="treeLength"/>
      <parameter idref="clock.rate"/>
      <rateStatistic idref="meanRate"/>
{subst_log}
      <parameter idref="frequencies"/>{demo_log}
    </log>
    <logTree id="treeFileLog" logEvery="{tree_every}" nexusFormat="true" fileName="{prefix}.trees" sortTranslationTable="true">
      <treeModel idref="treeModel"/>
      <trait name="rate" tag="rate"><strictClockBranchRates idref="branchRates"/></trait>
      <joint idref="joint"/>
    </logTree>
  </mcmc>
  <report>
    <property name="timer"><mcmc idref="mcmc"/></property>
  </report>
</beast>
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="dated alignment, headers accession|YYYY-MM-DD")
    ap.add_argument("--model", required=True, choices=list(MODEL_NAMES))
    ap.add_argument("--tag", default="1", help="chain tag appended to the model name")
    ap.add_argument("--out-suffix", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--name", default=None, help="file stem; default <model-name><tag>_<suffix>")
    ap.add_argument("--file-prefix", default=None,
                    help="prefix written into the XML for .log/.trees/.ops (default: the file stem, "
                         "i.e. relative to the directory BEAST is started in)")
    ap.add_argument("--chain-length", type=int, default=12_000_000)
    ap.add_argument("--log-every", type=int, default=5_000)
    ap.add_argument("--tree-every", type=int, default=60_000)
    ap.add_argument("--precision", type=float, default=4.06, help="fixed GMRF precision (skygrid-fixed)")
    ap.add_argument("--num-parameters", type=int, default=20)
    ap.add_argument("--cutoff", type=float, default=0.8)
    ap.add_argument("--subst-model", choices=["hky", "gtr"], default="hky")
    ap.add_argument("--root-height-min", type=float, default=None,
                    help="sensitivity: lower bound on the root height (years before the last tip)")
    ap.add_argument("--logpop-lower", type=float, default=None,
                    help="sensitivity: lower bound on every Skygrid log population size")
    ap.add_argument("--screen-every", type=int, default=None)
    ap.add_argument("--tip-date-text", choices=["fixed6", "shortest"], default="fixed6",
                    help="fixed6: 2026.397260 (12 Sep); shortest: 2026.39726 (bundle of 1 Sep)")
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)

    SEQ = read_fasta(args.fasta)
    labels = list(SEQ)
    dates = pd.to_datetime(pd.Series([t.split("|")[-1] for t in labels]), format="ISO8601", errors="coerce")
    assert dates.notna().all(), "every header must end with |YYYY-MM-DD"
    assert len(set(len(s) for s in SEQ.values())) == 1
    decimals = [f"{decimal_year(d):.6f}" for d in dates]
    if args.tip_date_text == "shortest":
        decimals = [str(float(x)) for x in decimals]
    name = args.name or f"{MODEL_NAMES[args.model]}{args.tag}_{args.out_suffix}"
    prefix = args.file_prefix or name
    txt = xml(args.model, labels, decimals, SEQ, chain=args.chain_length, log_every=args.log_every,
              tree_every=args.tree_every, prefix=prefix, precision=args.precision,
              num_parameters=args.num_parameters, cutoff=args.cutoff, subst_model=args.subst_model,
              root_height_min=args.root_height_min, logpop_lower=args.logpop_lower,
              screen_every=args.screen_every)
    fn = os.path.join(args.outdir, name + ".xml")
    with open(fn, "w") as fh:
        fh.write(txt)
    rec = dict(file=os.path.basename(fn), sha256=hashlib.sha256(txt.encode()).hexdigest(), model=args.model,
               tip_date_text=args.tip_date_text, n_taxa=len(labels), first_tip=str(dates.min().date()), last_tip=str(dates.max().date()),
               last_tip_decimal=decimal_year(dates.max()), chain_length=args.chain_length,
               log_every=args.log_every, tree_every=args.tree_every, screen_every=(args.screen_every or args.log_every * 10),
               burnin_states_30pct=int(args.chain_length * 0.3),
               skygrid_precision=("estimated; prior Gamma(shape 0.001, scale 1000); start value 0.1"
                                  if args.model == "skygrid-estimated"
                                  else (f"fixed at {args.precision}" if args.model == "skygrid-fixed" else None)),
               num_parameters=args.num_parameters if args.model != "exp" else None,
               cutoff=args.cutoff if args.model != "exp" else None,
               subst_model=args.subst_model, root_height_min=args.root_height_min,
               logpop_lower=args.logpop_lower,
               run_command=f"beast -threads 4 -seed <S> -overwrite {name}.xml",
               summary_tree_command_12sep=(f"treeannotator -burnin {int(args.chain_length * 0.3)} -heights median "
                                           f"{prefix}.trees {prefix}.hipstr.nexus   (HIPSTR: default type of "
                                           "TreeAnnotator 10.5)"),
               mcc_command=(f"treeannotator -type mcc -burnin {int(args.chain_length * 0.3)} -heights median "
                            f"{prefix}.trees {prefix}.mcc.nexus"))
    with open(os.path.join(args.outdir, f"beast_xml_parameters_{name}.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    print(f"[beast_make_xml] {fn}: {args.model}, {len(labels)} taxa, last tip {rec['last_tip']}, "
          f"precision {rec['skygrid_precision']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
