# Code and parameter settings: the reproduction number of the Bundibugyo virus disease outbreak in Ituri, estimated from genomes

## What this repository is

This repository holds the code and the parameter settings of the analyses that the short version of the article
reports: a post of the INRB on the reproduction number of the outbreak of Bundibugyo virus disease in Ituri,
estimated from genomes with Delphy and compared with BEAST X. It holds

- the programs that selected the genomes, ran the chains, summarised them and drew the figures;
- the settings of every analysis: the command line of every Delphy chain, the XML file of every BEAST X analysis
  (without the sequences), the schedules of the chains, the periods, and the rules by which chains were kept;
- the tables from which the numbers and the figures of the short version (Figure 1 and Figure 2) are read
  (`results/`).

It holds no genome sequence, no alignment, no output of a chain and no tree, and it does not hold Delphy or
BEAST X themselves. The estimates are in the tables of `results/` and in the article; this text states none.

Licence: Apache License, Version 2.0; its text is in `LICENSE`. It covers what this repository holds. The genome
sequences and the case counts are not in the repository and stay under the terms of their sources (see 'What is
not included, and why' and 'The case counts').

## What differs from the code as it ran

Every program and settings file is, byte for byte, the file that ran or that the analysis stored, with these
exceptions: paths of the machines of the analysis were made relative or replaced by a placeholder; look-ups of
stored files and identifiers of stored files were replaced by relative paths or by the names of the files; the
sequences were removed from the XML files of BEAST X, where the accession of every genome stays in its place;
tables whose stored names began with `INTERNAL_` (`figure1b_response_events_20261006.csv`, and
`main_figures_tokens_20261004.csv`, which is not in the repository) carry their names without it, also in the
lines of the programs and in the cells of the tables that name them; and the table with the outcome of the
quality checks of the consensus
sequences is named `public_quality_screen_<date>.csv` in the lines of the programs that name it. No name of a
machine or of a user had to be replaced. In tables of `results/` the cells that held identifiers of stored files
are empty (columns `source_version` and `source_version_id`).

Comments and documentation strings were edited in 40 of the programs: sentences about the organisation of the work
and about the computers on which it ran were removed or reworded, and version histories were shortened to the
technical differences between the versions. These edits touch no line of
code: for every edited program the syntax tree of the code without its documentation strings and comments is the
same before and after the edit. Lines of code were changed in these places: in
`04_summaries_and_reproduction_number/beast_full/code/make_config_v4.py` the import of the module of the
configuration, which was renamed `analyses_config.py`; in
`07_figures/code_supplementary_figures_20261006/supp_fig_common_20261006.py` the path from which the common module
of the figures is imported (the copy of that module in the same folder was removed); and in
`07_figures/fig_ne_R_cases_short_version_20261006.py`, `fig_ne_R_cases_short_version_event_legend_20261006.py` and
`fig_R_robustness_focus_short_version_periods_20261006.py` the texts of five error messages.
Further lines of code were changed in `07_figures/code_main_figures_20261006/main_fig_common_20261006.py`,
`fig_ne_R_cases_20261006.py`, `code_supplementary_figures_20261006/supp_fig_common_20261006.py`,
`fig_program_trajectories_20261006.py`,
`04_summaries_and_reproduction_number/beast_full/code/compare_programs_v3.py`,
`beast_sets/code/derive_beast_sets.py`, `pipeline_round_20261002/work/pipeline_ext/beast_first_delivery.py` and
`01_snapshot_and_genome_selection/membership/lineage_of_membership_of_every_genome_in_every_analysed_set_20261003.py`:
the texts of notes and messages that the programs write, the names of the two header forms of the table of the
response events (`first form`, `second form`) and the name of one column of the table of inputs of the figures
(`state_of_the_file`); and in `07_figures/add_case_counts_to_figure_1_values.py` the SHA-256 against which the file
of values of Figure 1 is compared. The notes that these programs wrote were shortened in the same way in the cells
that hold them: three cells of the column `note` of
`results/fig_ne_R_cases_20261006_values_without_case_counts.csv`, 41 cells of the column `drawn_as` of
`results/fig_program_trajectories_20261006_values.csv` and 189 cells of the column `note_on_the_pair` of
`results/program_comparison_trajectories_20261002.csv`; no value of a table changed.
In the shell scripts
`02_delphy/run_chain.sh`, `03_beast_x/beast_sets/run_chain_v2.sh`, `03_beast_x/beast_sets/config.sh` and
`03_beast_x/beast_full/code/build_sets_round.sh` comments were reworded in the same way. The column
`code_that_ran_the_chains` of `settings/analyses.csv` names the command files and wrappers of this repository in
place of the queue runners, and its column `code_that_summarised_the_chains` names the version of each program that
is in the repository (see 'What is not included, and why').

The table with the outcome of the quality checks is written by the program of the quality checks, which is not
in the repository. As it ran, that program writes the table under another name, so that the three programs that
call it and then read the table (`03_beast_x/beast_full/code/build_sets_round.sh`,
`05_sensitivity_analyses/code/build_sets_sensitivity.py` and
`05_sensitivity_analyses/lineage_of_genome_sets_sensitivity_20261002.py`) run as they stand only with a program
of the quality checks that writes the table under the name `public_quality_screen_<date>.csv`.

## Software

The versions are those of `settings/software_versions.csv`, which is a part of the version record of the analysis
(the rows of the programs, libraries, interpreters and imported Python packages of its environments `beast` and
`python`, and the rows of matplotlib and pillow of the environment `python`).

| software | version | row of `settings/software_versions.csv` (columns category, item, environment) |
|---|---|---|
| Delphy | 1.4.1 (build 2056, commit 1e8d4a9) | program, Delphy (programs/delphy_1.4.1/bin/delphy) |
| BEAST X | v10.5.0 | program, BEAST X (lib/beast.jar), beast |
| BEAGLE | v4.0.1 (PRE-RELEASE) | library, BEAGLE (libhmsbeagle.so.1), beast |
| Java runtime | 25.0.2-internal | program, Java runtime (bin/java), beast |
| Python | 3.11.15 | interpreter, Python, python (and beast) |
| numpy | 2.4.6 | python package (import), numpy, python (and beast) |
| pandas | 2.3.3 | python package (import), pandas, python (in the environment beast: 3.0.3) |
| scipy | 1.17.1 | python package (import), scipy, python (and beast) |
| matplotlib | 3.11.0 | package (conda), matplotlib, python |
| pillow | 12.3.0 | package (conda), pillow, python |

The file also lists the other programs of the environment `beast`: the tools `logcombiner` and `treeannotator` of
BEAST X, which `beast_analysis.py` calls for the trees, and IQ-TREE and MAFFT, which no program of this repository
calls.
The figures set the font Liberation Sans (`07_figures/code_main_figures_20261006/main_fig_common_20261006.py`,
line 610); with another font or another version of matplotlib the pictures differ in their pixels.

## What each folder holds

| folder | holds |
|---|---|
| `data/` | `genomes_by_analysed_set.csv`: every genome of the outbreak in the snapshot by accession and version, and for every analysed set whether the genome is in it; `snapshot_fetch_log_20261002.csv`: the queries with which the snapshot was downloaded from Pathoplexus, with the version of the data and the SHA-256 of every answer |
| `01_snapshot_and_genome_selection/` | `snapshot/`: the programs that froze the download and built the metadata and the alignment of the outbreak; `pipeline_round_20261002/`: the builder of the alignments and the specifications of the genome sets; `curation/`: the programs that wrote the tables of the sets; `regional_set/`: the functions that wrote the table of the regional set, which is not in the repository; `membership/`: the code that wrote the table of the membership of every genome in every set |
| `02_delphy/` | `run_chain.sh`, the wrapper that started one Delphy chain of the genome sets `C95`, `C90`, `C80` and `C00`: the command line, the limits and the environment of the process |
| `03_beast_x/` | `pipeline_round_20261002/work/pipeline_ext/beast_make_xml.py`, the generator of the XML files; `beast_sets/`: the wrapper that ran one part of one chain of the settings A and B, its configuration and the tool that wrote the XML of a continuation; `beast_full/code/`: the command file that built the genome sets of the settings C and D and the program that prepared the continuation of a chain from a saved state |
| `04_summaries_and_reproduction_number/` | `pipeline_round_20261002/work/pipeline_ext/`: the functions that read the chains, apply the rules by which chains are kept, and compute the reproduction number (`analysis.py`, `rt_lib.py`, and for BEAST X `beast_analysis.py`, `beast_first_delivery.py`); the programs that call them and write the tables of estimates and of chains |
| `05_sensitivity_analyses/` | the programs that built the genome sets of the sensitivity analyses and summarised their chains; `specs/`: the specifications of these sets and the command that built them |
| `06_case_counts_and_response_events/` | the program that built the weekly case counts, with its description (`README_case_series_20261002.md`); the code that wrote the table of the response events of Figure 1; `regenerate_case_tables.py`, written for this repository, which regenerates the case tables from the public source |
| `07_figures/` | the programs that draw Figure 1 and Figure 2, the programs that write the files of values from which they are drawn, and the program that adds the case counts to the file of values of Figure 1 |
| `settings/` | the parameter settings of every analysis in one place (see below) |
| `results/` | the tables from which the numbers and the figures of the short version are read, and the files of values of the figures; that of Figure 1 is held without its rows of case counts (see 'The case counts') |
| `MANIFEST.csv` | every file of the repository with its size in bytes and its SHA-256 |
| `LICENSE` | the text of the Apache License, Version 2.0 |

Folders named `pipeline_round_20261002/` (in `01_`, `03_` and `04_`) are parts of one working tree: the programs
find `work/pipeline_ext/`, `pipeline_20260924/` and `set_specs/` below one root, so copy the three parts into one
folder before a program is run that takes this tree as an argument. In the same way `beast_sets/` and `beast_full/`
are each one working tree that is divided between `03_` and `04_`.

## The settings of every analysis

`settings/analyses.csv` has one row for each of the 27 analyses: program and version, genome set (and its column
in `data/genomes_by_analysed_set.csv`), number of genomes, chains and steps, the files of `settings/` that hold
its commands and its XML, the code that ran its chains, the code that summarised them, and its table of estimates.

| analysis | program | genome set | genomes | chains | steps per chain |
|---|---|---|---|---|---|
| beast_A_C00 | BEAST X 10.5.0 | C00 | 680 | 8 | 40000000 |
| beast_A_C80 | BEAST X 10.5.0 | C80 | 674 | 8 | 40000000 |
| beast_B_C00 | BEAST X 10.5.0 | C00 | 680 | 8 | 40000000 |
| beast_C | BEAST X 10.5.0 | C00 | 680 | 8 | 40000000 |
| beast_C_retained1046 | BEAST X 10.5.0 | retained1046 | 515 | 8 | 40000000 |
| beast_D | BEAST X 10.5.0 | C00 | 680 | 8 | 40000000 |
| delphy_estimated_C00 | Delphy 1.4.1 | C00 | 680 | 16 | 1500000000 |
| delphy_exponential | Delphy 1.4.1 | C00 | 680 | 6 | 1500000000 |
| delphy_fixed_C00 | Delphy 1.4.1 | C00 | 680 | 24 | 1500000000 |
| delphy_fixed_C80 | Delphy 1.4.1 | C80 | 674 | 24 | 1500000000 |
| delphy_fixed_C90 | Delphy 1.4.1 | C90 | 555 | 24 | 1500000000 |
| delphy_fixed_C95 | Delphy 1.4.1 | C95 | 392 | 24 | 1500000000 |
| delphy_sens_allClasses | Delphy 1.4.1 | allClasses | 729 | 6 | 1500000000 |
| delphy_sens_cells1000 | Delphy 1.4.1 | C00 | 680 | 24 | 1500000000 |
| delphy_sens_cells16000 | Delphy 1.4.1 | C00 | 680 | 12 | 1500000000 |
| delphy_sens_drcOnly | Delphy 1.4.1 | drcOnly | 659 | 6 | 1500000000 |
| delphy_sens_dupKeepBest | Delphy 1.4.1 | dupKeepBest | 711 | 6 | 1500000000 |
| delphy_sens_dupNone | Delphy 1.4.1 | dupNone | 739 | 6 | 1500000000 |
| delphy_sens_gridCalendar20 | Delphy 1.4.1 | C00 | 680 | 6 | 1500000000 |
| delphy_sens_gridCoarse10 | Delphy 1.4.1 | C00 | 680 | 6 | 1500000000 |
| delphy_sens_gridWeekly32 | Delphy 1.4.1 | C00 | 680 | 6 | 1500000000 |
| delphy_sens_minusExcluded1046 | Delphy 1.4.1 | minusExcluded1046 | 616 | 6 | 1500000000 |
| delphy_sens_retained1046 | Delphy 1.4.1 | retained1046 | 515 | 6 | 1500000000 |
| delphy_sens_smooth15d | Delphy 1.4.1 | C00 | 680 | 6 | 1500000000 |
| delphy_sens_smooth60d | Delphy 1.4.1 | C00 | 680 | 6 | 1500000000 |
| delphy_sens_upTo0817 | Delphy 1.4.1 | upTo0817 | 615 | 6 | 1500000000 |
| delphy_sens_upTo0831 | Delphy 1.4.1 | upTo0831 | 668 | 6 | 1500000000 |

The numbers of chains are the numbers of rows of the analysis in its schedule; the steps of a chain are read from
the file that the column `steps_per_chain_read_from` of `settings/analyses.csv` names.

| in `settings/` | holds |
|---|---|
| `schedules/chain_schedule_round_20261002.csv` | the chains of the analyses of the genome sets `C95`, `C90`, `C80` and `C00`: program, setting, genome set, seed, steps, logging intervals, threads and burn-in share of every chain |
| `schedules/phase5_schedule_20261003.csv` | the chains of the other analyses: what is varied, and the seed of every chain |
| `schedules/lineage_of_*.py` | the code that wrote these schedules (see 'Code that was not stored as a program') |
| `delphy/<analysis>/<chain>.cmd` | the command line of the chain as it ran, with every option and the seed; `<chain>.env`: the limits and the environment settings of its process |
| `beast_x/<analysis>/xml/` | the XML file of the analysis without the sequences, the record of its parameters (`beast_xml_parameters_*.json`) and, for the settings C and D, the command that wrote it |
| `beast_x/<analysis>/<chain>/command.txt` | the command line of the chain with its seed; where a chain was continued from a saved state, a folder below it holds the command, the XML file and the record of the continuation |
| `beast_x/beast_full_control/queue.json` | the chains of the settings C and D and of the Delphy analysis with exponential growth: setting, genome set, seed, steps, logging intervals, command, working directory and environment of every chain |
| `commands_sensitivity_20261002.csv` | the command lines of the sensitivity analyses as they were planned, option by option, beside the command of the primary analysis |
| `periods_round_20261002.csv` | the periods for which reproduction numbers are given |
| `rules_and_constants_in_the_code.csv` | where the rules and constants stand that are fixed in the code, each with file, line and the text of the line |
| `software_versions.csv` | the versions of the software |

Settings that are fixed in the code (`settings/rules_and_constants_in_the_code.csv` quotes the lines):

- burn-in: the first 0.30 of the logged states of every chain are discarded (`analysis.py`, line 27; the column
  `burn_in_share` of `settings/schedules/chain_schedule_round_20261002.csv` holds the same share, 0.3);
- a chain is kept unless a logged state has a root date on or after the day before the first collection date
  (`analysis.py`, lines 27, 68 and 69), and, for Delphy, unless the density of the coalescent prior that the
  program computes exceeds the exact density of a sampled tree by more than 100 log units (`analysis.py`, lines
  607 and 633); the tables `chains_beast_sets_20261002.csv`, `chains_beast_full_20261002.csv` and
  `chains_sensitivity_20261002.csv` of `results/` say for every chain of their analyses whether it was kept, and
  why not; for the analyses `delphy_fixed_*` and `delphy_estimated_C00` the tables of estimates give the numbers
  of chains run and kept (columns `chains_run` and `chains_retained`);
- an estimate counts as converged with a pooled effective sample size of at least 200 and a split R-hat of at
  most 1.05 (`analysis.py`, line 495);
- generation time: gamma distribution with a mean of 15.3 days and a standard deviation of 9.3 days (`rt_lib.py`,
  line 194); the other generation times of the sensitivity analysis are in
  `results/generation_times_sensitivity_20261002.csv`.

`analysis.py` and `rt_lib.py` are in `04_summaries_and_reproduction_number/pipeline_round_20261002/work/pipeline_ext/`.
The genome sets are defined with the code that builds them: the sets `C95`, `C90`, `C80` and `C00` in
`01_snapshot_and_genome_selection/pipeline_round_20261002/set_specs/` (`README_set_specs.md` there gives their
thresholds of the called fraction), the sets of the sensitivity analyses in `05_sensitivity_analyses/specs/`;
`data/genomes_by_analysed_set.csv` gives the members of every set.
The SHA-256 values that settings files quote for XML files are those of the XML files with their sequences.

## The steps in their order, with the command of each

A program is started as its head says; the commands below are quoted from these heads. `TREE` is the working tree
`pipeline_round_20261002/` (see above). The chains were run from a run root, a folder that holds `programs/`
(Delphy), `aln/` (one alignment for each genome set, `aln/<set>_20261002.fasta`) and `runs/`; the command files
of `settings/` are written relative to it.

**0. The genomes.** Fetch the genomes of `data/genomes_by_analysed_set.csv` from Pathoplexus by accession and
version; `data/snapshot_fetch_log_20261002.csv` gives the queries of the download of 2 October 2026 (version of
the data: 1790610472). A later download can differ from the snapshot.

**1. Snapshot and genome selection** (`01_snapshot_and_genome_selection/`).

    python code/01_freeze.py                 # snapshot/: checks the download and freezes it
    python code/03_build_20261002.py         # snapshot/: metadata, alignment of the outbreak, validation table
    python write_restricted_list.py STORED_TABLE NEW_METADATA OUT.txt            # curation/code/
    python round_20261002/make_completeness_specs.py --frozen pipeline_20260924/alignment_set_spec_20260924.json --outdir set_specs
    python build_completeness_sets.py --pipeline TREE/pipeline_20260924 --specs TREE/set_specs --metadata META.csv
    python sets_independent.py --metadata META --alignment ALN --screen SCREEN --pairs PAIRS --restricted A,B --suffix SFX --outdir OUT
    python make_set_tables.py --sets-dir RUN/sets --independent RUN/independent --metadata META --alignment ALN --validation VAL

`01_freeze.py` and `03_build_20261002.py` read `code/inputs.json`, which is not in the repository because it held
paths of the machine and identifiers of stored files. It is a JSON object with the keys `paths` and `version_ids`;
the programs read `paths.zip_20261002` (the zip with the download), `version_ids.zip_20261002` (a label of that
zip, which is written into a record) and `paths.old_metadata`, `paths.old_fasta` and `paths.old_zone_map` (files of
the earlier snapshot). The program of the quality checks of the consensus sequences is not in the repository, and
neither is the table that it writes with its outcome; the programs that build the genome sets call the program
or read that table, which their lines name `public_quality_screen_<date>.csv` (see 'What differs from the code
as it ran'). The genomes of every analysed set are given in `data/genomes_by_analysed_set.csv`, and every later
step starts from them. `build_alignments.py` of `pipeline_20260924/` is called by `build_completeness_sets.py`,
whose head gives its further options. The functions of `regional_set/` wrote the table of the regional set,
which is not in the repository; they were called interactively, and the calling lines were not stored.

**2. Delphy chains of the genome sets `C95`, `C90`, `C80` and `C00`** (`02_delphy/`), from the run root, with the
alignments that step 1 writes copied to `aln/`. One chain of the schedule `chain_schedule_round_20261002.csv`:

    bash run_chain.sh <analysis> <genome_set> <seed> <fixed|estimated> <attempt>

The command line that the wrapper builds is the one of `settings/delphy/<analysis>/<chain>.cmd`, and the limits it
sets are those of `<chain>.env`; a chain can also be run by that command line alone. The program that went through
the schedule and started the wrapper once per chain is not in the repository; a chain that was ended from outside
was started again from state 0 with the same seed.

**3. BEAST X** (`03_beast_x/`).

    beast_make_xml.py --fasta cleaned_<suffix>.fasta --model skygrid-fixed --tag 1      # pipeline_round_20261002/work/pipeline_ext/; further options in its head
    bash run_chain_v2.sh <root> <label> <analysis> <genome_set> <setting> <seed> <threads> <xml-stem> fresh 1    # beast_sets/: one chain of the settings A and B
    python3 tools/make_continuation_xml.py <XML of the chain> <XML of the continuation> <length of the chain> <state of the checkpoint>
    prepare_continuation.py --root ROOT --control control --chain CHAIN [--check-seconds 600]    # beast_full/code/: settings C and D

The options with which each XML file was written are in its `beast_xml_parameters_*.json` (and, for the settings C
and D, in `*_make_xml_command.txt`) in `settings/beast_x/<analysis>/xml/`. The command line of every chain is in
`settings/beast_x/<analysis>/<chain>/command.txt`, to be run in the directory of the chain; the settings C and D
and the Delphy analysis with exponential growth have, in addition, `settings/beast_x/beast_full_control/queue.json`
with the working directory and the environment of every chain. `beast_sets/config.sh` and `queue.json` name the
folder of the launcher of BEAST X as `PATH_OF_THE_BEAST_X_ENVIRONMENT/bin`; put your own path there. The programs
that went through the schedules and started the chains are not in the repository. A chain that was ended from
outside was continued from its last saved state: for the settings A and B by `run_chain_v2.sh` in its mode
`continue`, with the XML of `tools/make_continuation_xml.py`; for the settings C and D by
`prepare_continuation.py`. A saved state was used only if the difference between the saved and the recomputed
log-likelihood that BEAST X reports was below 1e-6; `settings/beast_x/<analysis>/<chain>/` holds the command, the
XML and the record of every continuation.

**4. Summaries and reproduction number** (`04_summaries_and_reproduction_number/`).

    python derive_round3.py --analysis delphy_fixed_C00 --run-root RUN_ROOT --pipeline TREE/work/pipeline_ext --outdir OUT
    python tools/join_continued_v2.py <root> <analysis> <outdir> [--steps 40000000 --log-every 20000 --tree-every 400000]
    python tools/cut_continued.py <root> <analysis> <outdir> [--steps 40000000 --log-every 20000] [--runs runs]
    python code/derive_beast_sets.py --root <working directory> --tree <root of the working tree> --out <dir> --analyses a,b,... [--test]
    make_config_v4.py --analysis NAME --extracted DIR --alignment A.fasta --periods P.csv --pipeline work/pipeline_ext --out config.json
    derive_analysis_v5.py --config analysis.json --out DIR
    continued_chains_table.py --extracted DIR --analysis NAME --out continued_chains_NAME_20261002.csv
    compare_programs_v3.py --pairs pairs.json --periods periods.csv --out-values V.csv --out-trajectories T.csv

`delphy_sets/derivation/derive_round3.py` wrote the tables of the set `C00`; the tables of the sets `C95`, `C90`
and `C80` were written by earlier versions of the same program, which computed the same estimates and differed in
bookkeeping columns of the table of the chains (see 'What is not included, and why'). In the same way
`beast_full/code/derive_analysis_v5.py` with `make_config_v4.py` summarised the analyses of the settings C and D
and the Delphy analysis with exponential growth; the analysis `beast_C` was summarised by the version before it.
The file `pairs.json` that `compare_programs_v3.py` reads names, for every pair of analyses that is compared, the
tables of estimates and of intervals of the two analyses; it is not in the repository because it held identifiers
of stored files. `derive_beast_sets.py` is run from the root of the working tree `pipeline_round_20261002/`.

**5. Sensitivity analyses** (`05_sensitivity_analyses/`). The genome sets are built with the command of
`specs/build_command_20261002.txt`; the chains are run with the command lines of
`settings/commands_sensitivity_20261002.csv` (column `command_line`), from the run root and with the limits of
the `.env` file of the chain; every analysis is summarised by `process(analysis, periods_csv, schedule_csv, ...)` of
`code/postprocess_sensitivity_v3.py`, which calls `code/derive_sensitivity_v3.py`; then

    python combine_sensitivity.py SCHEDULE.csv   (from the working directory of the sensitivity analyses)
    python quantities_without_chains.py POSTERIOR.csv.gz ALIGNMENT.fasta PERIODS.csv PIPELINE_DIR OUTDIR

and `main(out_of_derive_sensitivity, periods_csv, pipeline, outdir)` of `code/primary_derived_quantities.py`.

**6. Case counts and response events** (`06_case_counts_and_response_events/`): the repository holds no table of
case counts. They are regenerated from the public source (see 'The case counts' below), from the top folder of the
repository:

    git clone https://github.com/INRB-UMIE/BDBV2026-Data
    git -C BDBV2026-Data checkout 8eb57154cf967e7dfdb7e690230f61f2a298d2b0
    python 06_case_counts_and_response_events/regenerate_case_tables.py BDBV2026-Data

`regenerate_case_tables.py` was written for this repository. It runs `build_case_series_20261002.py`, the program
of the analysis, with the arguments of the analysis and writes its tables to `results/case_counts_regenerated/`; it
copies the weekly table (`weekly_confirmed_cases_by_province_20261002.csv`) to `results/`; and it calls
`07_figures/add_case_counts_to_figure_1_values.py`, which writes the file of values of Figure 1
(`results/fig_ne_R_cases_20261006_values.csv`). It compares the weekly table and the file of values with the SHA-256
of the files of the analysis and prints whether they are the same. The cumulative counts by province that the article
prints are the last row of `results/case_counts_regenerated/cumulative_confirmed_cases_by_province_20261002.csv`.
The files that this step writes are not listed in `MANIFEST.csv`.

The command file of the analysis is not in the repository. It ran `build_case_series_20261002.py` with the
arguments that `regenerate_case_tables.py` passes and, in addition, with `--source-note` and a text that names the
source; that text fills the column `source` of `case_growth_seven_zones_treatments_20261002.csv` and the entry
`source_note` of `case_series_checks_20261002.json` and enters no other file. Its two further commands called a
program that compared the numbers of the article with the case counts and a program that wrote reports; neither
is in the repository.

**7. Figures** (`07_figures/`), from the top folder of the repository:

    MPLBACKEND=Agg PYTHONHASHSEED=0 python 07_figures/draw_figure_1.py OUTDIR
    MPLBACKEND=Agg PYTHONHASHSEED=0 python 07_figures/draw_figure_2.py OUTDIR

The first command draws Figure 1 (`OUTDIR/fig_ne_R_cases_short_version_event_legend_20261006.png` and `.pdf`) from
`results/fig_ne_R_cases_20261006_values.csv`, which step 6 writes, and
`results/fig_program_trajectories_20261006_values.csv`. The second draws Figure 2
(`OUTDIR/fig_R_robustness_focus_short_version_periods_20261006.png` and `.pdf`) from
`results/fig_R_robustness_focus_20261006_values.csv`. `draw_figure_1.py` and `draw_figure_2.py` were written for this
repository and only call the programs of the two figures:
`fig_ne_R_cases_short_version_event_legend_20261006.py`, which draws through
`fig_ne_R_cases_short_version_20261006.py`, `fig_ne_R_cases_lanes_20261006.py` and
`code_main_figures_20261006/fig_ne_R_cases_20261006.py`; and
`fig_R_robustness_focus_short_version_periods_20261006.py`, which draws through
`code_main_figures_20261006/fig_R_robustness_focus_20261006.py`. Drawn from this repository with the versions of the
software given above, each picture was byte for byte the picture of the article. After the edits described under 'What differs from the
code as it ran', step 6 wrote the weekly table with the SHA-256 of the analysis, and Figure 2 was drawn with the
same bytes as from the unedited repository; Figure 1 was not drawn again, because the font was not available on the
machine used.

The files of values are written by `code_main_figures_20261006/fig_ne_R_cases_20261006.py` and
`fig_R_robustness_focus_20261006.py` (Figures 1 and 2) and by
`code_supplementary_figures_20261006/fig_program_trajectories_20261006.py` (the lines of BEAST X in Figure 1) from
the tables of `results/`:

    python fig_ne_R_cases_20261006.py INPUTS.json OUTDIR --interventions NAME.csv [--sources KEYS.csv]   the figure of the article
    python fig_R_robustness_focus_20261006.py INPUTS.json [OUTDIR]          builds the file of values, then draws
    python fig_program_trajectories_20261006.py INPUTS_SUPPLEMENT.json INPUTS_THIRD_TASK.json [OUTDIR]

`INPUTS.json` maps the name of every input file to its path. This step cannot be run from the repository alone:
beside the tables of `results/` the programs open the table of their inputs, the specification of the numbers of
the article, the table of the tokens of the figures (`main_figures_tokens_20261004.csv`), the document of the
conventions of the figures and a module of an earlier round; these are working files of the article and are not
in the repository.

## Code that was not stored as a program

Five tables were written by code that was typed interactively and not stored as a program. For
these the repository holds the code as it was typed, in files named `lineage_of_<table>.py`:
the schedules (`settings/schedules/`), the table of the membership of the genomes in the sets
(`01_snapshot_and_genome_selection/membership/`), the table of the genome sets of
the sensitivity analyses (`05_sensitivity_analyses/`) and the table of the response events
(`06_case_counts_and_response_events/`). Each was run
again when the repository was made, in the form in which it stands here and with its input files under `inputs/`,
and wrote its table byte for byte; for the code of the genome sets of the sensitivity analyses this was done with
a copy of the program of the quality checks that writes its table under the name
`public_quality_screen_<date>.csv`. They read their input files from a folder
`inputs/`; these are files of the analysis, most of which are not in the repository, and
`lineage_of_phase5_schedule_20261003.py` and `lineage_of_figure1b_response_events_20261006.py` start from an
earlier version of their own table. For `results/periods_round_20261002.csv` no such code could be included.

## What is not included, and why

- **Genome sequences and alignments.** They are fetched from Pathoplexus by accession. Pathoplexus has terms of
  use for open and for restricted-use data (their addresses are in `data/snapshot_fetch_log_20261002.csv`), which a
  user of the genomes has to respect.
- **The program of the quality checks of the consensus sequences, and the tables of its outcome.** The program
  is not in the repository, and neither is the table that it writes with its outcome, which the programs that
  build the genome sets read and which their lines name `public_quality_screen_<date>.csv` (see 'What differs
  from the code as it ran'). The tables `genome_sets_by_called_fraction_20261002.csv` and
  `regional_set_20261002.csv`, which give the outcome of the quality checks genome by genome, are not in the
  repository either. The genomes of every analysed set are given in `data/genomes_by_analysed_set.csv`, and
  every later step starts from them.
- **The output of the chains and the trees.** The chains need many processor hours: the table above gives the
  chains of every analysis and the steps of each chain. No table of the repository holds a run time for the Delphy
  chains of the genome sets `C95`, `C90`, `C80` and `C00` (the analyses `delphy_fixed_*` and
  `delphy_estimated_C00`). For the other analyses the column `wall_time_s` of
  `results/chains_beast_sets_20261002.csv`, `results/chains_beast_full_20261002.csv` and
  `results/chains_sensitivity_20261002.csv` gives the wall time of a chain; for a chain that was continued from a
  saved state it is empty in the first of these tables and is the time of the continuation alone in the second
  (its column `start_end_wall_time_exit_code_and_command_line_are_those_of` says which).
- **Delphy and BEAST X.** Their versions are given above; the launcher and the binaries are not part of the code.
- **The programs that went through the schedules and started the chains** (queue runners, launchers and the
  keepers of saved states), and the queues they read. The command line of every chain is in `settings/`; the
  wrappers that set the environment of a process and the rules for continuing a chain from a saved state are in
  `02_delphy/` and `03_beast_x/`. A chain that was ended from outside was started again from state 0 (Delphy) or
  continued from its last saved state (BEAST X).
- **Checks, tests, second computations, comparisons with an earlier round, reports and logs**, and the programs
  that assembled the tables of a work step for the next. A program of this kind is in the repository only where a
  program that is included imports it or reads what it writes (`sets_independent.py`, `check_sensitivity.py`,
  `chain_status.py`).
- **Earlier versions of programs.** Where a program was run in several versions, the repository holds the last.
  The earlier versions of `derive_round3.py` (they wrote the tables of the sets `C95`, `C90` and `C80`) computed the
  same estimates and differed in bookkeeping columns of the table of the chains and in a comparison file that they
  did not write; the version before `derive_analysis_v5.py` (it summarised `beast_C`) lacked the columns
  `start_end_wall_time_exit_code_and_command_line_are_those_of` and the text 'not applicable' of `reached_last_step`
  for a chain cut at its saved state, and `make_config_v4.py` differs from the version before it in the same
  columns; the versions before `derive_sensitivity_v3.py` and `postprocess_sensitivity_v3.py` (they summarised
  `delphy_sens_allClasses`) had the same code apart from comments and the name of the module that is imported.
- **Working files of the article**: its text, its templates and the tables of its numbers, among them the table of
  the tokens of the figures (`main_figures_tokens_20261004.csv`), which two files of values of `results/` name as
  the source of a label, and the program that compared the numbers of the article with the case counts.

Tables of `results/` also hold rows of analyses that the short version does not report (for example
`beast_A_C90`); the settings of those analyses are not part of the repository.

## What runs from the repository alone, and what `results/` is for

From the repository alone, with the software above: the command of step 7 that draws Figure 2 from its file of
values in `results/`. Figure 1 needs, in addition, the case counts of the source below (step 6, then the first
command of step 7). Every other step needs the genomes (steps 1 to 5), the output of the chains (steps 4 and 5) or
working files of the article (the files of values of step 7).

`results/` can be deleted without touching any program of the analysis outside `07_figures/`: `regenerate_case_tables.py`,
written for this repository, writes into it, and no other program of the other folders
reads it. Without it the figures cannot be drawn, and the tables that the steps above write are no longer at
hand for comparison. The schedules, the table of the periods and `commands_sensitivity_20261002.csv` are in
`settings/` as well as in `results/`.

## The case counts

The repository holds no table of case counts. The weekly counts of confirmed cases of Figure 1c come from the
situation reports of the INSP as transcribed in the public repository INRB-UMIE/BDBV2026-Data
(https://github.com/INRB-UMIE/BDBV2026-Data), commit 8eb57154 of 1 October 2026 (files retrieved 2 October 2026 for
the analysis), and step 6 regenerates them from that commit.
`06_case_counts_and_response_events/README_case_series_20261002.md` describes the source and how the weekly series
was built from it, and quotes the licence of that repository and its terms for the situation reports: reuse with
attribution to the INSP and citation of the number and date of the report, and the request to confirm the terms of
distribution with the INSP before the counts are published elsewhere.

## Authors, acknowledgements, funding and competing interests

The statements below are copied from the article; 'this post' in them is the article. Where the article cites the
transcription of the situation reports by the number of a reference, the name and the address of that repository
are written out here.

**Authors.** In alphabetical order; grouped by institution:

- Collins Tanui (Africa CDC)
- Justus Nsio (Africa CDC)
- Yap Boum II (Africa CDC)
- Yenew Kebede (Africa CDC)
- Laura Luebbert (Anthropic, USA)
- David Blazes (Gates Foundation)
- Sofonias Kifle Tessema (Gates Foundation)
- Adrienne Amuri-Aziza (INRB and Basic Sciences Department at the University of Kinshasa, Kinshasa DRC; Department
  of clinical sciences at the Institute of Tropical Medicine Antwerp Belgium and Department of Microbiology,
  Immunology and Transplantation KU Leuven Belgium, Belgium)
- Dav Ebengo (INRB, INOHA, Kinshasa, DRC)
- Emmanuel Lokilo-Lofiko (INRB, Kinshasa, DRC)
- Gradi Luakanda-Ndelemo (INRB, Kinshasa, DRC)
- Prince Akil (INRB, Kinshasa, DRC)
- Sifa Kavira (INRB, Kinshasa, DRC)
- Daniel Mukadi (INRB, University of Kinshasa, Kinshasa, DRC)
- Eddy Kinganda-Lusamaki (INRB, University of Kinshasa, Kinshasa, DRC; TransVIHMI, Université de Montpellier,
  INSERM, IRD, Montpellier, France)
- Pauline Muswamba (INRB, University of Kinshasa, Kinshasa, DRC)
- Placide Mbala-Kingebeni (INRB, University of Kinshasa, Kinshasa, DRC; Africa CDC)
- Steve Ahuka-Mundeke (INRB, University of Kinshasa, Kinshasa, DRC)
- Tania Bishola (INRB, University of Kinshasa, Kinshasa, DRC)
- Tony Wawina-Bokalanga (INRB, University of Kinshasa, Kinshasa, DRC; Department of Clinical Sciences, Institute of
  Tropical Medicine, Antwerp, Belgium)
- Dieudonne Mwamba (INSP, Kinshasa, DRC)
- Pierre Akilimali (INSP, University of Kinshasa, Kinshasa, DRC)
- Amadou Mouctar Diallo (WHO, DRC)
- Olga Ntumba (WHO, DRC)
- Nicksy Gumede (WHO Regional Office for Africa)

**Claude Science.** All analyses were performed with the assistance of Claude Science (Anthropic), which was used
to run the phylodynamic analyses, estimate reproduction numbers, and generate figures and tables; all outputs were
reviewed by the authors, and the code and parameter settings are available at
https://github.com/inrb-labgenpath/inrb_bdbv_claude_science.

**Acknowledgements and funding.** We gratefully acknowledge the laboratories, health authorities, and response
teams that generated the genome sequences and shared them rapidly through Pathoplexus, and the teams of the
Institut National de Santé Publique that compile the situation reports, and the INRB-UMIE team that transcribes
them in INRB-UMIE/BDBV2026-Data (https://github.com/INRB-UMIE/BDBV2026-Data). We thank the Ministry of Public
Health, Hygiene and Social Welfare of the DRC. The authors gratefully acknowledge the ongoing support of the Africa
Centers for Disease Control and Prevention (Africa CDC), the World Health Organization, and partner
non-governmental organizations. We also acknowledge the support provided by the University of Manitoba (Canada),
South African National Bioinformatics Institute (SANBI), Osaka Metropolitan University (OMU), Unité de Gestion du
Programme de Développement du Système de Santé (UG-PDSS), the Institute of Tropical Medicine (ITM) through Belgian
Directorate-general for Development Cooperation and Humanitarian Aid (DGD FA5 project, the Culmen International
LCC, the US CDC Atlanta, and the Agence Française de Développement through the AFROSCREEN project (grant agreement
CZZ3209), coordinated by ANRS Maladies Infectieuses émérgentes in partnership with Institut Pasteur and Institut de
Recherche pour le Développement, the French Ministry of Europe and Foreign Affairs through FEF programme and
support provided by Institut de Recherche pour le Développement. AREBO project funded by ANRS-MIE. A.A.-A is
supported by a DGD sandwich PhD scholarship. P.M-K. acknowledges the support of the Wellcome Trust through the
ARTIC Network (award 313694/Z/24/Z) and the Gates Foundation. This work was supported by Anthropic, which provided
the credits and compute used for the analyses. Claude Science carried out the phylodynamic runs, the estimation of
reproduction numbers, and the production of figures and tables under the authors’ direction, and the analyses
reported here would not have been completed in this time frame without it.

**Competing interests.** L.L. is an employee of Anthropic and holds equity in the company. Anthropic develops
Claude Science, the software used in this analysis. L.L. contributed to the analysis, interpretation, and writing
of this post; scientific outputs were reviewed by all authors. The remaining authors declare no competing
interests.
