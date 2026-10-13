# Code and parameter settings for 'Reproduction number of the 2026 Bundibugyo virus disease outbreak over time: phylodynamic estimates from genomes collected up to 12 September 2026'

This repository holds the code and the parameter settings of the analyses that the article reports: a post of the INRB
on the reproduction number of the outbreak of Bundibugyo virus disease in Ituri, estimated from genomes with Delphy and
compared with BEAST X. The analyses rest on the genomes of the download of 7 October 2026 from Pathoplexus. The
repository holds

- the programs that selected the genomes, summarised the chains and drew the figures, and the wrapper that started
  one Delphy chain;
- the settings of every analysis: the command line of every Delphy chain, the XML file of every BEAST X analysis
  (without the sequences), the table of the chains of every analysis as they ran, the periods, and the rules by
  which chains were kept;
- the tables from which the numbers and the figures of the article (Figure 1 with its four panels and Figure 2) are
  read (`results/`).

It holds no genome sequence, no alignment, no output of a chain and no tree, and it does not hold Delphy or BEAST X
themselves. The estimates are in the tables of `results/` and in the article; this text states none.

Licence: Apache License, Version 2.0; its text is in `LICENSE`. Copyright 2026 Institut National de Recherche
Biomédicale (INRB). The licence covers what this repository holds. The genome sequences and the case counts are not in
the repository and stay under the terms of their sources (see 'What is not included, and why' and 'The case counts').

## What this repository is for

The purpose of this repository is to show the parameters and the decisions that matter, so that the analyses can be
understood, audited, and reproduced: which genomes were analysed and why others were left out, the model and the
command line of every chain, the rules by which chains were kept, the periods, the generation time, and the tables
from which every number and every figure of the article was read.

This repository is not a pipeline: the genomes, the output of the chains, the control files that issued the calls,
and working files of the article are not included (see 'What is not included, and why'). Paths of the machines on
which the analyses ran were made relative or replaced by placeholders (such as `PATH_OF_THE_BEAST_X_ENVIRONMENT`),
the sequences were removed from the XML files of BEAST X, and comments in the programs were shortened; no line of
code that computes a value was changed. In the files of values of the figures and in the two tables of the comparison
of the programs, the column that named a stored version of a table is empty. 'What runs from the repository alone'
says which commands run as they stand.

The analyses were carried out with Claude Science (Anthropic), which wrote and ran the code under the direction of
the authors. We recommend reproducing them the same way, with Claude Science pointed at this repository: it reads the
settings and the programs here, can fetch the genomes from Pathoplexus by accession, run Delphy and BEAST X with the
command lines of `settings/`, rebuild the tables and the figures, and answer questions about any file of the
repository.

## Reproducing and evaluating analyses with Claude Science

Claude Science is a desktop application for scientific research. It is included in the Pro, Max, Team and Enterprise
plans of Claude (not in the Free plan; on a Team or Enterprise plan an owner of the organisation enables it first).
The documentation is at <https://claude.com/docs/claude-science>.

**Install.**

- macOS (13 or later): download the installer from <https://claude.com/product/claude-science> and double-click it.
  On first launch the app sets up its runtime and its Python and R environments, which takes a few minutes, then
  opens in a browser tab.
- Linux (x64): install `curl`, `bubblewrap` (0.8.0 or later) and `socat` with the package manager of your
  distribution (for example `sudo apt-get install -y curl bubblewrap socat`), then

      curl -fsSL https://claude.ai/install-claude-science.sh | bash
      claude-science serve

- Windows (11, x64): download the installer from the same page.

Sign in with your Claude account; no API key is needed.

**Point it at this repository and ask a question.**

1. Connect the GitHub connector, or get the repository onto your computer:
   `git clone https://github.com/inrb-labgenpath/inrb_bdbv_claude_science`, or download the ZIP from GitHub and
   unpack it.
2. In Claude Science, create a project and start a conversation. Name the folder in your message (type its path, or
   insert it with `@`) and allow the folder-access card when it appears; allow the code-execution card when Claude
   proposes to run code.
3. Ask. For example:

   > I have cloned the repository `inrb_bdbv_claude_science` to `~/inrb_bdbv_claude_science`. It holds the code and
   > the parameter settings of the phylodynamic estimates of the reproduction number of the 2026 Bundibugyo virus
   > disease outbreak; `README.md` describes it. Please read `README.md` and `settings/analyses.csv`, then:
   > (1) list, for the primary analysis `delphy_fixed_C00`, the genome set, the Skygrid grid, the number of cells
   > of the coalescent prior, the smoothing, the number of chains and of steps per chain, the burn-in, and the rules
   > by which chains were kept, each with the file and line where it stands; (2) confirm from
   > `results/estimates_delphy_fixed_C00_20261007.csv` that the reproduction numbers for the periods P1, P2 and P3
   > are 1.97 (95 % HPD 1.61–2.34), 1.26 (1.01–1.53) and 0.83 (0.70–0.97); and (3) draw Figure 2 with
   > `python 07_figures/draw_figure_2.py out` and show me the picture.

   To rerun an analysis, ask it to fetch the genomes of `data/genomes_by_analysed_set.csv` from Pathoplexus by
   accession (respecting the terms of use linked in `data/snapshot_fetch_log_20261007.csv`), to build the alignment
   as step 1 below describes, and to run the chains with the command lines of `settings/`. The chains need many
   processor hours (column `wall_time_s` of the tables of chains in `settings/schedules/`).

## Software

The versions are those of `settings/software_versions.csv`. Every row of that file names the file of the analysis in
which the version stands (column `version_source`): the tables of chains and the files `environment.txt` of
`settings/`, a log of the versions of Python and of its packages, and the records of the runs that drew the figures;
the last two are not in the repository.

| software | version | row of `settings/software_versions.csv` (columns category, item, environment) |
|---|---|---|
| Delphy | 1.4.1 | program, Delphy |
| BEAST X | 10.5.0 | program, BEAST X (lib/beast.jar), beast |
| BEAGLE | 4.0.1 | library, BEAGLE (libhmsbeagle.so.1), beast |
| Java runtime | 25.0.2-internal | program, Java runtime, beast |
| Python | 3.11.15 | interpreter, Python, python |
| numpy | 2.4.6 | python package, numpy, python |
| pandas | 2.3.3 | python package, pandas, python |
| scipy | 1.17.1 | python package, scipy, python |
| matplotlib | 3.11.0 | python package, matplotlib, python |
| pillow | 12.3.0 | python package, pillow, python |

The file gives the SHA-256 of the program of Delphy and of the files of BEAST X and BEAGLE that the chains ran. The
figures set the font Liberation Sans (`07_figures/code_main_figures_20261007/main_fig_common_20261007.py`, line
646); with another font or another version of matplotlib the pictures differ in their pixels.

## What each folder holds

| folder | holds |
|---|---|
| `data/` | `genomes_by_analysed_set.csv`: every genome that is in at least one analysed set, by accession and version, and for every analysed set whether the genome is in it; `BDBV_post_update_seqset_accessions_20261008.csv`: the accessions of the set of sequences that the article cites, with the role of each (focal: the analysed set; background); `snapshot_fetch_log_20261007.csv`: the queries with which the snapshot was downloaded from Pathoplexus, with the version of the data and the SHA-256 of every answer |
| `01_snapshot_and_genome_selection/` | `snapshot/`: the programs that froze the download, built the metadata and the alignment of the outbreak, and wrote the list of the restricted-use records of other groups; `pipeline_round_20261002/`: the builder of the alignments and the specifications of the genome sets; `curation/`: the programs that prepared the input of the quality checks and wrote the tables of the sets |
| `02_delphy/` | `run_chain_r2.sh`, the wrapper that started one Delphy chain of the genome sets `C95`, `C90`, `C80` and `C00`: the command line, the limits and the environment of the process |
| `03_beast_x/` | `pipeline_round_20261002/work/pipeline_ext/beast_make_xml.py`, the generator of the XML files |
| `04_summaries_and_reproduction_number/` | `pipeline_round_20261002/work/pipeline_ext/`: the functions that read the chains, apply the rules by which chains are kept, and compute the reproduction number (`analysis.py`, `rt_lib.py`, and for BEAST X `beast_analysis.py`, `beast_first_delivery.py`); `delphy_sets/code/derive_round_r2.py` and `beast_update/derive/derive_beast_update.py`, the programs that call them and write the tables of estimates and of chains, and `beast_update/tools/build_derive_root.py`, which builds the working root of the second; `beast_full/code/`: the program that compares the two programs, and `compare_programs_update.py`, written for this repository, which writes its input files and runs it |
| `05_sensitivity_analyses/` | the programs that built the genome sets of the further analyses and summarised their chains; `specs/`: the specification of these sets, the command that built them and the periods; the code that computed the table of the pairs of adjacent intervals of the analysis `delphy_sens_gridWeekly32` (see 'Code that was not stored as a program') |
| `06_case_counts_and_response_events/` | the programs that built the weekly case counts; the code that wrote the table of the response events of Figure 1; `regenerate_case_tables.py`, written for this repository, which regenerates the case tables from the public source |
| `07_figures/` | the programs that draw Figure 1 and Figure 2, the programs that write the files of values from which they are drawn, and the program that adds the case counts to the file of values of Figure 1 |
| `settings/` | the parameter settings of every analysis in one place (see below) |
| `results/` | the tables from which the numbers and the figures of the article are read, and the files of values of the figures; that of panels a to c of Figure 1 is held without its rows of case counts (see 'The case counts') |
| `LICENSE` | the text of the Apache License, Version 2.0 |

Folders named `pipeline_round_20261002/` (in `01_`, `03_` and `04_`) are parts of one working tree: the programs
find `work/pipeline_ext/`, `pipeline_20260924/` and `set_specs/` below one root, so copy the three parts into one
folder before a program is run that takes this tree as an argument.

## The settings of every analysis

`settings/analyses.csv` has one row for each of the 33 analyses: program and version, genome set (and its column in
`data/genomes_by_analysed_set.csv`), number of genomes, chains and steps, the table of its chains, the files of
`settings/` from which the chains were taken, the files of `settings/` that hold its commands and its XML, the code
that ran its chains, the code that summarised them, and its table of estimates.

| analysis | program | genome set | genomes | chains | steps per chain |
|---|---|---|---|---|---|
| beast_A_C00 | BEAST X 10.5.0 | C00 | 775 | 8 | 40000000 |
| beast_B_C00 | BEAST X 10.5.0 | C00 | 775 | 8 | 40000000 |
| delphy_estimated_C00 | Delphy 1.4.1 | C00 | 775 | 16 | 1500000000 |
| delphy_fixed_C00 | Delphy 1.4.1 | C00 | 775 | 24 | 1500000000 |
| delphy_fixed_C80 | Delphy 1.4.1 | C80 | 769 | 24 | 1500000000 |
| delphy_fixed_C90 | Delphy 1.4.1 | C90 | 641 | 24 | 1500000000 |
| delphy_fixed_C95 | Delphy 1.4.1 | C95 | 471 | 24 | 1500000000 |
| delphy_rem_buniaThinned_rep1_run2 | Delphy 1.4.1 | buniaThinned_rep1 | 706 | 6 | 1500000000 |
| delphy_rem_buniaThinned_rep2_run2 | Delphy 1.4.1 | buniaThinned_rep2 | 706 | 6 | 1500000000 |
| delphy_rem_buniaThinned_rep3_run2 | Delphy 1.4.1 | buniaThinned_rep3 | 706 | 6 | 1500000000 |
| delphy_rem_randomThinnedA_rep1 | Delphy 1.4.1 | randomThinnedA_rep1 | 706 | 6 | 1500000000 |
| delphy_rem_randomThinnedA_rep2 | Delphy 1.4.1 | randomThinnedA_rep2 | 706 | 6 | 1500000000 |
| delphy_rem_randomThinnedA_rep3 | Delphy 1.4.1 | randomThinnedA_rep3 | 706 | 6 | 1500000000 |
| delphy_rem_randomThinnedB_rep1 | Delphy 1.4.1 | randomThinnedB_rep1 | 665 | 6 | 1500000000 |
| delphy_rem_randomThinnedB_rep2 | Delphy 1.4.1 | randomThinnedB_rep2 | 665 | 6 | 1500000000 |
| delphy_rem_randomThinnedB_rep3 | Delphy 1.4.1 | randomThinnedB_rep3 | 665 | 6 | 1500000000 |
| delphy_rem_withoutRecentBunia_run2 | Delphy 1.4.1 | withoutRecentBunia | 665 | 6 | 1500000000 |
| delphy_sens_allClasses | Delphy 1.4.1 | allClasses | 841 | 6 | 1500000000 |
| delphy_sens_cells1000 | Delphy 1.4.1 | C00 | 775 | 24 | 1500000000 |
| delphy_sens_cells16000 | Delphy 1.4.1 | C00 | 775 | 12 | 1500000000 |
| delphy_sens_cells400 | Delphy 1.4.1 | C00 | 775 | 24 | 1500000000 |
| delphy_sens_drcOnly | Delphy 1.4.1 | drcOnly | 754 | 6 | 1500000000 |
| delphy_sens_dupKeepBest | Delphy 1.4.1 | dupKeepBest | 806 | 6 | 1500000000 |
| delphy_sens_dupNone | Delphy 1.4.1 | dupNone | 834 | 6 | 1500000000 |
| delphy_sens_gridCalendar20 | Delphy 1.4.1 | C00 | 775 | 6 | 1500000000 |
| delphy_sens_gridCoarse10 | Delphy 1.4.1 | C00 | 775 | 6 | 1500000000 |
| delphy_sens_gridWeekly32 | Delphy 1.4.1 | C00 | 775 | 6 | 1500000000 |
| delphy_sens_minusExcluded1046 | Delphy 1.4.1 | minusExcluded1046 | 715 | 6 | 1500000000 |
| delphy_sens_retained1046 | Delphy 1.4.1 | retained1046 | 515 | 6 | 1500000000 |
| delphy_sens_smooth15d | Delphy 1.4.1 | C00 | 775 | 6 | 1500000000 |
| delphy_sens_smooth60d | Delphy 1.4.1 | C00 | 775 | 6 | 1500000000 |
| delphy_sens_upTo0817 | Delphy 1.4.1 | upTo0817 | 666 | 6 | 1500000000 |
| delphy_sens_upTo0831 | Delphy 1.4.1 | upTo0831 | 727 | 6 | 1500000000 |


files; the steps of a chain are read from the file that the column `steps_per_chain_read_from` of
`settings/analyses.csv` names.


|---|---|
| `schedules/chains_<analysis>_20261007.csv` | the table of the chains of the analysis as they ran (`chains_beast_update_20261007.csv` for the two analyses of BEAST X): for every chain the seed, the command line, the last step or state logged, the states used after the part discarded, and whether the chain was kept, and why |
| `delphy/<analysis>/<chain>.cmd` | the command line of the chain as it ran, with every option and the seed; `<chain>.env`: the limits and the environment settings of its process |
| `beast_x/<analysis>/xml/` | the XML file of the analysis without the sequences, the record of its parameters (`beast_xml_parameters_*.json`) and the calls of the generator that wrote the XML files (`calls_that_wrote_the_xml_files_20261007.txt`) |
| `beast_x/<analysis>/<chain>/` | `command.txt`: the command line of the chain with its seed; `environment.txt`: the environment of the chain at its start and the program environment (packages, SHA-256 of the launcher and of the libraries, version of Java) |
| `schedules/chain_schedule_round_20261002.csv` | the schedule from which the chains and the seeds of seven analyses were taken (`delphy_fixed_C95`, `delphy_fixed_C90`, `delphy_fixed_C80`, `delphy_fixed_C00`, `delphy_estimated_C00`, `beast_A_C00`, `beast_B_C00`); the tables of chains say what ran; its other rows are not used here |
| `schedules/phase5_schedule_20261003.csv` | the schedule of the further analyses: what is varied in each analysis, and the seed of every chain. The calls of `process(...)` of `postprocess_sensitivity_20261007_v4.py` gave it this file as `schedule_csv`; the file also holds rows of analyses that are not part of the article |
| `schedules/lineage_of_*.py` | the code that wrote these two schedules (see 'Code that was not stored as a program') |
| `commands_sensitivity_20261007.csv` | the command lines of the further analyses as they were planned, option by option, beside the command of the primary analysis, in the order in which they were queued |
| `periods_of_the_update_20261007.csv`, `ends_by_set_20261007.csv` | the periods for which reproduction numbers are given, and for each of the four genome sets the end of the last period |
| `periods_round_20261002.csv` | the periods with their stored names (`W1` to `W3`), from which the programs of the figures read the names of the periods |
| `rules_and_constants_in_the_code.csv` | where the rules and constants stand that are fixed in the code, each with file, line and the text of the line |
| `software_versions.csv` | the versions of the software |

The queue runners and the control files by which they started the chains are not in the repository (see 'What is not
included, and why'). `commands_sensitivity_20261007.csv` also holds the first runs of four genome sets
(`delphy_rem_buniaThinned_rep1` to `_rep3` and `delphy_rem_withoutRecentBunia`). The analyses of the table above
whose names end in `_run2` are the second runs of these sets with other seeds (the rule by which these analyses
are read needs the difference between two runs of the same set; `settings/schedules/phase5_schedule_20261003.csv`,
column `description`). The article reports the second runs; the settings and tables of the first runs are not in
the repository.

The two XML files hold, for each of the 775 genomes, the accession and the collection date as the analysis used it,
with the sequences removed. The SHA-256 values that settings files quote for XML files are those of the XML files
with their sequences.

Settings that are fixed in the code (`settings/rules_and_constants_in_the_code.csv` quotes the lines):

- burn-in: the first 0.30 of the logged states of every chain are discarded (`analysis.py`, line 27);
- a chain is kept unless a logged state has a root date on or after the day before the first collection date
  (`analysis.py`, lines 27, 68 and 69), and, for Delphy, unless the density of the coalescent prior that the program
  computes exceeds the exact density of a sampled tree by more than 100 log units (`analysis.py`, lines 607 and
  633); the tables of chains in `settings/schedules/` say for every chain whether it was kept, and why (columns
  `retained` and `why`);
- an estimate counts as converged with a pooled effective sample size of at least 200 and a split R-hat of at most
  1.05 (`analysis.py`, line 495);
- generation time: gamma distribution with a mean of 15.3 days and a standard deviation of 9.3 days (`rt_lib.py`,
  line 194); the other generation times are in `results/generation_times_sensitivity_20261007.csv`.

`analysis.py` and `rt_lib.py` are in
`04_summaries_and_reproduction_number/pipeline_round_20261002/work/pipeline_ext/`. The genome sets are defined with
the code that builds them: the sets `C95`, `C90`, `C80` and `C00` in
`01_snapshot_and_genome_selection/pipeline_round_20261002/set_specs/` (`README_set_specs.md` there gives their
thresholds of the called fraction), the sets of the further analyses in `05_sensitivity_analyses/specs/` and in
`results/removal_sets_parameters_sensitivity_20261007.json`; `data/genomes_by_analysed_set.csv` gives the members of
every set, and `results/genome_sets_sensitivity_20261007.csv` and `results/set_steps_by_set_20261007.csv` give the
number of genomes of every set.

## The steps in their order, with the command of each

A program is started as its head says. The commands below are the calls with which the tables of the repository were
written, with every argument; the paths in them are those of the working directories of the analysis, and the control
files that issued them are not in the repository (see 'What is not included, and why'). `TREE` is the working tree
`pipeline_round_20261002/` (see above). The chains were run from a run root, a folder that holds `programs/`
(Delphy), `aln/` (one alignment for each genome set, `aln/<set>_20261007.fasta`) and `runs/`; the command files
of `settings/` are written relative to it. Every chain ran in one attempt (columns `attempts` and `started_again`
of the tables of chains in `settings/schedules/`).

**0. The genomes.** Fetch the genomes of `data/genomes_by_analysed_set.csv` from Pathoplexus by accession and
version; `data/snapshot_fetch_log_20261007.csv` gives the queries of the download of 7 October 2026 (version of the
data: 1791281853). A later download can differ from the snapshot.

**1. Snapshot and genome selection** (`01_snapshot_and_genome_selection/`).

    python code/01_freeze.py                 # snapshot/: checks the download and freezes it
    python code/03_build_20261007.py         # snapshot/: metadata, alignment of the outbreak, validation table
    python code/05_terms_of_use.py           # snapshot/: the list of the restricted-use records of other groups
    python code/make_screen_inputs.py in_new/pp_snapshot_20261007_metadata.csv in_new/pp_snapshot_20261007_outbreak_aligned.fasta.gz in_new/restricted_use_records_of_other_groups_20261007.txt 20261007 screen_in
    python build_completeness_sets.py --pipeline pipeline_copy --specs specs_used --metadata screen_in/INTERNAL_screen_input_metadata_20261007.csv --alignment screen_in/screen_input_outbreak_aligned_20261007.fasta.gz --screen runs/new_snapshot/run_v2/public_quality_screen_20261007.csv --batch-table runs/new_snapshot/run_v2/screen_batches_20261007.csv --duplicate-pairs in_new/duplicate_sample_pairs_20261007.csv --data-use-exclude "" --validation in_new/snapshot_validation_20261007.csv --suffix 20261007 --outdir runs/new_snapshot/sets
    python sets_independent.py --metadata screen_in/INTERNAL_screen_input_metadata_20261007.csv --alignment screen_in/screen_input_outbreak_aligned_20261007.fasta.gz --screen runs/new_snapshot/run_v2/public_quality_screen_20261007.csv --pairs in_new/duplicate_sample_pairs_20261007.csv --restricted "" --suffix 20261007 --outdir runs/new_snapshot/independent
    python code/make_set_tables.py --sets-dir runs/new_snapshot/sets --independent runs/new_snapshot/independent --metadata in_new/pp_snapshot_20261007_metadata.csv --alignment in_new/pp_snapshot_20261007_outbreak_aligned.fasta.gz --validation in_new/snapshot_validation_20261007.csv --zone-map in_new/zone_harmonisation_map_20261007.csv --restricted in_new/restricted_use_records_of_other_groups_20261007.txt --suffix 20261007 --outdir deliver

The first three calls were run in `snapshot/` with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and
`NUMEXPR_NUM_THREADS` set to 4; the others in the working directory of the genome selection, with a limit of the
address space of 16 GiB (`ulimit -v 16777216`), the same four variables set to 1 and `PYTHONDONTWRITEBYTECODE=1`.
There `pipeline_copy` is a copy of `TREE/pipeline_20260924`; `specs_used/` holds the four specifications of the sets
and `alignment_set_spec_completeness_all_20261002.json` of `TREE/set_specs/`; `in_new/` holds the files that the
first three calls write and the table of the duplicate samples; and `runs/new_snapshot/run_v2/` holds the output of
the program of the quality checks. `build_completeness_sets.py` is in `pipeline_round_20261002/round_20261002/` and
`sets_independent.py` in `curation/code/`. After the last call the four alignments
`runs/new_snapshot/sets/<set>/<set>_20261007.fasta` (`C95`, `C90`, `C80`, `C00`) were copied to `deliver/`.
`03_build_20261007.py` calls `build_snapshot.py`, which calls the functions of `snapshot_lib.py`.

`01_freeze.py`, `03_build_20261007.py` and `05_terms_of_use.py` read `code/inputs.json`, which is not in the
repository because it held paths of the machine and identifiers of stored files. It is a JSON object with the keys
`paths` and `version_ids`, each of which maps a label to a path or to a label of a stored file; the programs read
`paths.zip_20261007` (the zip with the download) and four further paths (`paths.zone_map_20261002`,
`paths.metadata_20261002`, `paths.fasta_20261002`, `paths.terms_table_20261002`) of files that are not in the
repository. `03_build_20261007.py` also reads `out/control_facts.json`, the outcome of a control on its inputs, and
stops without it; `05_terms_of_use.py` reads a folder `control_on_20261002/` of the same control; neither is in the
repository. The records of the list that `05_terms_of_use.py` writes are removed from the metadata and from the
alignment by `make_screen_inputs.py` before the genome sets are built. The program of the quality checks of the
consensus sequences is not in the repository, and neither is the table that it writes with its outcome; the programs
that build the genome sets call the program or read that table, which their lines name
`public_quality_screen_<date>.csv`. The genomes of every analysed set are given in
`data/genomes_by_analysed_set.csv`, and every later step starts from them. `build_alignments.py` of
`pipeline_20260924/` is called by `build_completeness_sets.py`, whose head gives its further options; the
specifications of the sets in `set_specs/` were written by `round_20261002/make_completeness_specs.py`.

**2. Delphy chains of the genome sets `C95`, `C90`, `C80` and `C00`** (`02_delphy/`), from the run root, with the
alignments that step 1 writes copied to `aln/`. One chain of the schedule `chain_schedule_round_20261002.csv`:

    bash run_chain_r2.sh <analysis> <genome_set> <seed> <fixed|estimated> <attempt>

The command line that the wrapper builds is the one of `settings/delphy/<analysis>/<chain>.cmd`, and the limits it
sets are those of `<chain>.env`; a chain can also be run by that command line alone. The program that went through
the schedule and started the wrapper once per chain is not in the repository; it checked the program of Delphy
against its SHA-256 (`settings/software_versions.csv`) before it started a chain, and it would have started a chain
that was ended from outside again from state 0 with the same command and seed.

**3. BEAST X** (`03_beast_x/`).

    beast_make_xml.py --fasta cleaned_<suffix>.fasta --model skygrid-fixed --tag 1      # pipeline_round_20261002/work/pipeline_ext/; further options in its head

The calls with which the two XML files were written are in
`settings/beast_x/<analysis>/xml/calls_that_wrote_the_xml_files_20261007.txt`, and the options of each XML file in
its `beast_xml_parameters_*.json`. The command line of every chain is in
`settings/beast_x/<analysis>/<chain>/command.txt`, to be run in the directory of the chain with the XML file copied
into it. `environment.txt` beside it gives what the queue runner set for the chain: `JAVA_TOOL_OPTIONS` (two threads
for the parallel and one for the concurrent garbage collection), the folder of the launcher of BEAST X at the head
of `PATH`, and `OMP_NUM_THREADS` unset. It names the folder of BEAST X as `PATH_OF_THE_BEAST_X_ENVIRONMENT` (and the
folder that this path is a link to as `REAL_PATH_OF_THE_BEAST_X_ENVIRONMENT`); put your own path there. The queue
runner, which is not in the repository, started the 16 chains at once, each from state 0. It never continued a chain
from a saved state; a chain that was ended from outside would have been started again from state 0 with the same
command, and none was (`settings/schedules/chains_beast_update_20261007.csv`, columns `attempts` and
`started_again`).

**4. Summaries and reproduction number** (`04_summaries_and_reproduction_number/`).

    python3 code/derive_round_r2.py --analysis <analysis> --run-root run_root --pipeline tree/pipeline_round_20261002/work/pipeline_ext --outdir derived/<analysis>      # delphy_sets/
    python3 U/tools/build_derive_root.py U W <analysis>                                # beast_update/
    python3 U/derive/derive_beast_update.py --root W --tree W/tree/pipeline_round_20261002 --out W/parts --analyses <analysis> --periods W/checks/periods_of_the_update_20261007.csv --periods-by-set W/checks/ends_by_set_20261007.csv --schedule W/checks/chain_schedule_round_20261002.csv
    python3 U/derive/derive_beast_update.py --root W --tree W/tree/pipeline_round_20261002 --out W/parts --analyses <analyses with a summary, separated by commas> --periods W/checks/periods_of_the_update_20261007.csv --periods-by-set W/checks/ends_by_set_20261007.csv --schedule W/checks/chain_schedule_round_20261002.csv --constants-only
    python 04_summaries_and_reproduction_number/beast_full/code/compare_programs_update.py OUTDIR      # from the top folder of the repository

The first call is that of the Delphy chains of the four sets, once for every analysis, from a working directory that
holds `code/`, `run_root/` and `tree/`, with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and
`NUMEXPR_NUM_THREADS` set to 1.

The next three calls are those of the analyses of BEAST X, once for each of `beast_A_C00` and `beast_B_C00`; `U` is
the working directory of these analyses and `W` is `U/derive_root`. `build_derive_root.py` builds `W` from the
archive of the raw output of the analysis (`U/deliver/<analysis>_raw_output_20261007.tar.gz`), the queue of the
chains (`U/sched/queue.csv`, which is not in the repository), the alignment, the two tables of the periods and the
schedule (`U/stored/`), the records of the XML files (`U/xml/`) and the working tree
(`U/tree/pipeline_round_20261002`); it also copies three tools from `U/code_sets/tools/` (`join_continued_v2.py`,
`cut_continued.py`, `chain_status.py`), which are not in the repository and are needed only for a chain that was
continued from a saved state; none was. The third of these calls, with `--constants-only`, writes the table of the
constants that are passed for the analyses that have a summary. For a chain continued from a saved state
`derive_beast_update.py` would have been called three times in place of the second call, with its options
`--inputs`, `--exclude-chains`, `--cut-chains` and `--tag`; these calls did not run. In a test of the code on short
test chains `--test` was added to the call.

`compare_programs_update.py` was written for this repository. The two tables of the comparison of the two programs
(`results/program_comparison_values_20261007.csv` and `results/program_comparison_trajectories_20261007.csv`) were
written by `compare_programs_v3.py` from two input files that a module wrote which is not in the repository.
`compare_programs_update.py` writes the same two input files (the periods and the three pairs of analyses) from
`settings/periods_of_the_update_20261007.csv`, `settings/ends_by_set_20261007.csv` and the tables of estimates and
of intervals in `results/`, runs `compare_programs_v3.py` and says whether the two tables that it wrote are byte for
byte those of `results/`.

**5. The further analyses** (`05_sensitivity_analyses/`). The genome sets are built with the command of
`specs/build_command_20261007.txt` (`code/build_sets_sensitivity_v2.py`); `code/set_tables_20261007.py` writes the
tables of the sets. The chains are run with the command lines of `settings/commands_sensitivity_20261007.csv`
(column `command_line`), from the run root and with the limits of the `.env` file of the chain; every analysis is
summarised by `process(analysis, periods_csv, schedule_csv, ...)` of `code/postprocess_sensitivity_20261007_v4.py`,
which calls `code/derive_sensitivity_20261007_v4.py` and `code/check_sensitivity_20261007_v2.py`; then

    python combine_sensitivity_20261007.py SCHEDULE.csv   (from the working directory of the further analyses)
    python quantities_without_chains_20261007.py POSTERIOR.csv.gz ALIGNMENT.fasta PERIODS.csv PIPELINE_DIR OUTDIR

and `main(out_of_derive_sensitivity, periods_csv, pipeline, outdir)` of
`code/primary_derived_quantities_20261007.py`. `periods_csv` is `specs/periods_sensitivity_20261007.csv`.
`combine_sensitivity_20261007.py` reads the queue of the queue runner from `run/queue.csv` (its column `analysis`),
and `set_tables_20261007.py` takes the same file as `--queue`. The queue is not in the repository;
`settings/commands_sensitivity_20261007.csv` holds its analyses and chains row by row, in its columns
`order_in_the_queue`, `milestone`, `analysis`, `chain`, `genome_set`, `seed` and `command_line` (the queue names
them `order`, `milestone`, `analysis`, `chain`, `genome_set`, `seed` and `command`).

**6. Case counts and response events** (`06_case_counts_and_response_events/`): the repository holds no table of
case counts. They are regenerated from the public source (see 'The case counts' below), from the top folder of the
repository:

    git clone https://github.com/INRB-UMIE/BDBV2026-Data
    git -C BDBV2026-Data checkout 71e06bdaeda49c068ba499cab41359179cc4b9d0
    python 06_case_counts_and_response_events/regenerate_case_tables.py BDBV2026-Data

`regenerate_case_tables.py` was written for this repository. It builds the case series by the calls of the program
of the analysis (`case_series_update_20261007.py`): it checks that `build_case_series_20261002.py` is the script
whose SHA-256 that program holds, calls `read_zone_list(...)` and `build(...)` of that script on the folder given,
and writes the tables with `write_step1(...)` of the program of the analysis, to `results/case_counts_regenerated/`;
it copies the weekly table (`weekly_confirmed_cases_by_province_20261007.csv`) to `results/`; and it calls
`07_figures/add_case_counts_to_figure_1_values.py`, which writes the file of values of panels a to c of Figure 1
(`results/fig_ne_R_cases_20261007_values.csv`). It compares the SHA-256 of every table with that of the table of the
analysis and prints whether they are the same. `add_case_counts_to_figure_1_values.py` does the same for the file of
values; there the file of the analysis is meant in the form in which this repository holds such files (the cells
that named a stored version of a table empty, the tables named as the members of `results/` are named, three cells
of the column `note` shortened). The cumulative counts by province that the article prints are the last row of
`results/case_counts_regenerated/cumulative_confirmed_cases_by_province_20261007.csv`.

`case_series_update_20261007.py` is not run as a whole here. In the analysis it ran with `OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and `NUMEXPR_NUM_THREADS` set to 2, and beside the folder of the download
that holds the source and the fetch log (`--new-root`) it took inputs that are not in the repository: a folder with
the script, the list of health zones and tables against which it checks its output (`--old-root`, `--stored`,
`--earlier-suffix`), the folders of its output and of its working files (`--outdir`, `--workdir`), the file label
(`--suffix`), and the number of the situation report and the sentence of the article on the case counts
(`--report-of-the-draft`, `--draft-sentence`), from which it takes the number of the report and four counts and
holds them against the reports in the files. `regenerate_case_tables.py` does not call it with arguments; it imports
it and calls two of its functions, as said above.

**7. Figures** (`07_figures/`), from the top folder of the repository:

    MPLBACKEND=Agg PYTHONHASHSEED=0 OMP_NUM_THREADS=1 python 07_figures/draw_figure_1.py OUTDIR
    MPLBACKEND=Agg PYTHONHASHSEED=0 OMP_NUM_THREADS=1 python 07_figures/draw_figure_2.py OUTDIR

The first command draws Figure 1 with its four panels
(`OUTDIR/fig_ne_R_cases_short_version_zoom_update_20261007.png` and `.pdf`) from
`results/fig_ne_R_cases_20261007_values.csv`, which step 6 writes,
`results/fig_program_trajectories_pairs_A_and_B_estimated_20261007_values.csv` (the lines of BEAST X and of Delphy
with the smoothing estimated in panel a), and, for panel d, `results/primary_pairs_sensitivity_20261007.csv`,
`results/trajectories_delphy_sens_gridWeekly32_20261007.csv`,
`results/pairs_of_adjacent_intervals_delphy_sens_gridWeekly32_update_20261007_computed_20261008.csv` and
`results/figure1b_response_events_20261006.csv`. The program writes the values of panel d on every run
(`OUTDIR/fig_ne_R_cases_short_version_zoom_update_20261007_values_panel_d.csv`), and `draw_figure_1.py` says whether
they are those of `results/`. The second command draws Figure 2
(`OUTDIR/fig_R_robustness_focus_short_version_periods_20261007.png` and `.pdf`) from
`results/fig_R_robustness_focus_20261007_values.csv`. `draw_figure_1.py` and `draw_figure_2.py` were written for
this repository and only call the programs of the two figures with the arguments with which the pictures of the
article were drawn: `fig_ne_R_cases_short_version_zoom_update_20261007.py` draws Figure 1, through
`fig_ne_R_cases_short_version_event_legend_20261007.py`, `fig_ne_R_cases_short_version_20261007.py`,
`fig_ne_R_cases_lanes_update_20261007.py`, `code_main_figures_20261007/fig_ne_R_cases_20261007.py` (with
`main_fig_common_20261007.py`) and `rt_lib.py`; `fig_R_robustness_focus_short_version_periods_20261007.py` draws
Figure 2, through `code_main_figures_20261007/fig_R_robustness_focus_20261007.py`. The file of values of Figure 1
keeps the stored names of the periods, `W1`, `W2` and `W3` (in its columns `label`, `value` and `quantity_id`; the
column `stored_name` of `settings/periods_round_20261002.csv` gives them), and the picture prints `P1`, `P2` and
`P3`: `fig_ne_R_cases_short_version_event_legend_20261007.py` puts these names into the drawing (its constant
`NAMES`, line 61, and its function `periods_and_room`) and changes no cell of the file. In the same way the file of
values of Figure 2 keeps `W1` to `W3`, and `fig_R_robustness_focus_short_version_periods_20261007.py` gives its
picture the names `P1` to `P3`. `names_of_the_round_main_figures_20261007.json` gives the programs of Figure 1 the
names of the tables of `results/`. Drawn from this repository with the versions of the software given above, each
picture was byte for byte the picture of the article, as PNG and as PDF.

The table `results/pairs_of_adjacent_intervals_delphy_sens_gridWeekly32_update_20261007_computed_20261008.csv` was
computed on 8 October 2026 from the posterior of the analysis `delphy_sens_gridWeekly32`; the code that computed it
is
`05_sensitivity_analyses/lineage_of_pairs_of_adjacent_intervals_delphy_sens_gridWeekly32_update_20261007_computed_20261008.py`
(see 'Code that was not stored as a program').

The files of values are written by `code_main_figures_20261007/fig_ne_R_cases_20261007.py` and
`fig_R_robustness_focus_20261007.py` (panels a to c of Figure 1; Figure 2) and by
`code_supplementary_figures_20261007/fig_program_trajectories_20261007.py` (the lines of BEAST X and of Delphy with
the smoothing estimated in Figure 1) from the tables of `results/`:

    python fig_ne_R_cases_20261007.py INPUTS.json OUTDIR --interventions NAME.csv [--sources KEYS.csv]
    python fig_R_robustness_focus_20261007.py INPUTS.json [OUTDIR]
    python fig_program_trajectories_20261007.py INPUTS_SUPPLEMENT.json INPUTS_THIRD_TASK.json [OUTDIR]

`INPUTS.json` maps the name of every input file to its path. This step cannot be run from the repository alone:
beside the tables of `results/` and `settings/` the programs open the table of their inputs, the specification of
the numbers of the article, the table of the labels of the figures (`main_figures_tokens_20261004.csv`), the table
of the sources of the response events and the document of the conventions of the figures; these are working files of
the article and are not in the repository. In the analysis these programs were run with the environment variable
`MAIN_FIGURES_ROUND` set to the file of the names of the tables (`names_of_the_round_main_figures_20261007.json`),
and `fig_program_trajectories_20261007.py` with `SUPPLEMENT_FIGURES_ROUND` set to a second such file, which is not
in the repository.

## Code that was not stored as a program

Four tables of the repository were written by code that was typed interactively and not stored as a program. For
these the repository holds the code as it was typed, in files named `lineage_of_<table>.py`: the two schedules
(`settings/schedules/`), the table of the response events (`06_case_counts_and_response_events/`) and the table of
the pairs of adjacent intervals of the analysis `delphy_sens_gridWeekly32` (`05_sensitivity_analyses/`). Each was run
again afterwards with its input files under `inputs/` and wrote its table byte for byte. Most of these input files
are not in the repository, and `lineage_of_phase5_schedule_20261003.py` and
`lineage_of_figure1b_response_events_20261006.py` start from a previous state of their own table.

The code of the table of the pairs reads nine files from `inputs/`: an archive of the working tree
(`pipeline_round_20261002.tar.gz`), the alignment of the genome set `C00`, the tables of the posterior of the
primary analysis, `delphy_fixed_C00`, and of `delphy_sens_gridWeekly32`, the table of the intervals of
`delphy_sens_gridWeekly32`, the table of the pairs of the primary analysis (here
`results/primary_pairs_sensitivity_20261007.csv`), the table of the chains of `delphy_sens_gridWeekly32` (here
`settings/schedules/chains_delphy_sens_gridWeekly32_20261007.csv`), a first table of the pairs, from which it takes
the names of its columns and the entries `program` and `setting`, and a working file, from which it takes the
precision of the smoothing that the primary analysis holds fixed. Only the two tables named with their folders are in
the repository. The code makes a folder `work/zoom/tree` below the folder from which it is run, unpacks the archive
into it and loads `rt_lib.py` from it; that file stands here in
`04_summaries_and_reproduction_number/pipeline_round_20261002/work/pipeline_ext/`.

For `results/periods_round_20261002.csv`, `settings/periods_of_the_update_20261007.csv`,
`settings/ends_by_set_20261007.csv` (its column `computed_with` says how its dates were computed) and
`data/BDBV_post_update_seqset_accessions_20261008.csv` no such code is included.

## What is not included, and why

- **Genome sequences and alignments.** They are fetched from Pathoplexus by accession. Pathoplexus has terms of use
  for open and for restricted-use data (their addresses are in `data/snapshot_fetch_log_20261007.csv`), which a user
  of the genomes has to respect. Collection dates, places, terms of use and submitting groups of single genomes are
  not in the tables of the repository; they are fetched from Pathoplexus by accession. The one exception are the two
  XML files of BEAST X, which hold the collection date of each of their genomes as a setting of the analysis.
- **The program of the quality checks of the consensus sequences, and the tables of its outcome.** The program is
  not in the repository, and neither is the table that it writes with its outcome, which the programs that build the
  genome sets read and which their lines name `public_quality_screen_<date>.csv`. The table
  `genome_sets_by_called_fraction_20261007.csv`, which gives the outcome of the quality checks genome by genome, is
  not in the repository either. The genomes of every analysed set are given in `data/genomes_by_analysed_set.csv`,
  and every later step starts from them.
- **The output of the chains and the trees.** The chains need many processor hours: the table above gives the chains
  of every analysis and the steps of each chain, and the column `wall_time_s` of the tables of chains in
  `settings/schedules/` gives the wall time of every chain.
- **Delphy and BEAST X.** Their versions are given above; the launcher and the binaries are not part of the code.
- **The programs that went through the schedules and started the chains** (queue runners, launchers and the keeper
  of saved states), the queues and control files they read, and the programs that wrote these. The command line of
  every chain is in `settings/`; the wrapper that set the environment of a Delphy process is in `02_delphy/`, and
  the environment of every chain of BEAST X is in its `environment.txt`. A chain that was ended from outside would
  have been started again from state 0; none was.
- **The control files of the steps**, which gave the order of the programs of a step and the arguments of every
  call. 'The steps in their order, with the command of each' gives the calls with which they wrote the tables, for
  the programs that are here.
- **Checks, tests, second computations, comparisons, reports and logs**, the programs that assembled the tables of a
  step for the next, and the descriptions of working archives. A program of this kind is in the repository only where
  a program that is included imports it or reads what it writes (`sets_independent.py`,
  `check_sensitivity_20261007_v2.py`).
- **Other versions of programs.** Where a program exists in several versions, the repository holds the one that ran.
  Programs from which an included program was derived, and which its head may name, are not included.
- **Working files of the article**: its text, its templates and the tables of its numbers, among them the table of
  the labels of the figures (`main_figures_tokens_20261004.csv`), which two files of values of `results/` name as
  the source of a label.
- **Tables of case counts** (see 'The case counts').

## What runs from the repository alone, and what `results/` is for

From the repository alone, with the software above: the command of step 7 that draws Figure 2 from its file of
values in `results/`, and the command of step 4 that writes the two tables of the comparison of the two programs
again. Figure 1 needs, in addition, the case counts of the source below (step 6, then the first command of step 7).
Every other step needs the genomes (steps 1 to 5), the output of the chains (steps 4 and 5) or working files of the
article (the files of values of step 7). In step 1 two programs need, in addition, files of a control that is not in
the repository: `03_build_20261007.py` reads `out/control_facts.json` and stops without it, and `05_terms_of_use.py`
reads the folder `control_on_20261002/`.

`results/` can be deleted without touching any program of the analysis outside `07_figures/`: of the programs of the
other folders only `regenerate_case_tables.py` and `compare_programs_update.py`, both written for this repository,
write into it or read it. Without it the figures cannot be drawn, and the tables that the steps above write are no
longer at hand for comparison. The two schedules, `periods_round_20261002.csv`,
`commands_sensitivity_20261007.csv` and three of the tables of chains are in `settings/` as well as in `results/`.

## The case counts

The repository holds no table of case counts. The weekly counts of confirmed cases of panel c of Figure 1 come from
the situation reports of the INSP as transcribed in the public repository INRB-UMIE/BDBV2026-Data
(https://github.com/INRB-UMIE/BDBV2026-Data), state of 7 October 2026, commit 71e06bda, and step 6 regenerates them
from that commit. The licence of that repository is the MIT licence (its `LICENSE.md`); for the situation reports
its file `data/insp_sitrep/metadata.yaml` asks for reuse with attribution to the INSP and citation of the number and
date of the report, and to confirm the terms of distribution with the INSP before the counts are published
elsewhere.

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
- Adrienne Amuri-Aziza (INRB and Basic Sciences Department at the University of Kinshasa, Kinshasa DRC; Department of Clinical Sciences at the Institute of Tropical Medicine Antwerp Belgium and Department of Microbiology, Immunology and Transplantation KU Leuven, Belgium)
- Dav Ebengo (INRB, INOHA, Kinshasa, DRC)
- Emmanuel Lokilo-Lofiko (INRB, Kinshasa, DRC)
- Gradi Luakanda-Ndelemo (INRB, Kinshasa, DRC)
- Prince Akil (INRB, Kinshasa, DRC)
- Sifa Kavira (INRB, Kinshasa, DRC)
- Daniel Mukadi (INRB, University of Kinshasa, Kinshasa, DRC)
- Eddy Kinganda-Lusamaki (INRB, University of Kinshasa, Kinshasa, DRC; TransVIHMI, Université de Montpellier, INSERM, IRD, Montpellier, France)
- Pauline Muswamba (INRB, University of Kinshasa, Kinshasa, DRC)
- Placide Mbala-Kingebeni (INRB, University of Kinshasa, Kinshasa, DRC; Africa CDC)
- Steve Ahuka-Mundeke (INRB, University of Kinshasa, Kinshasa, DRC)
- Tania Bishola (INRB, University of Kinshasa, Kinshasa, DRC)
- Tony Wawina-Bokalanga (INRB, University of Kinshasa, Kinshasa, DRC; Department of Clinical Sciences, Institute of Tropical Medicine, Antwerp, Belgium)
- Dieudonne Mwamba (INSP, Kinshasa, DRC)
- Pierre Akilimali (INSP, University of Kinshasa, Kinshasa, DRC)
- Amadou Mouctar Diallo (WHO, DRC)
- Olga Ntumba (WHO, DRC)
- Nicksy Gumede (WHO Regional Office for Africa)

**Claude Science.** All analyses were performed with the assistance of Claude Science (Anthropic), which was used
to run the phylodynamic analyses, estimate reproduction numbers, and generate figures and tables; all outputs were
reviewed by the authors, and the code and parameter settings are available at
https://github.com/inrb-labgenpath/inrb_bdbv_claude_science.

**Acknowledgements and funding.** We gratefully acknowledge the laboratories, health authorities, and response teams
that generated the genome sequences and shared them rapidly through Pathoplexus, and the teams of the Institut
National de Santé Publique that compile the situation reports, and the INRB-UMIE team that transcribes them in
INRB-UMIE/BDBV2026-Data (https://github.com/INRB-UMIE/BDBV2026-Data). We thank the Ministry of Public Health, Hygiene
and Social Welfare of the DRC. We also thank Patrick Varilly and Andrew Rambaut for their comments and for reviewing
the analyses during the preparation of this post. The authors gratefully acknowledge the ongoing support of the
Africa Centers for Disease Control and Prevention (Africa CDC), the World Health Organization, and partner
non-governmental organizations. We also acknowledge the support provided by the University of Manitoba (Canada),
South African National Bioinformatics Institute (SANBI), Osaka Metropolitan University (OMU), Unité de Gestion du
Programme de Développement du Système de Santé (UG-PDSS), the Institute of Tropical Medicine (ITM) through Belgian
Directorate-general for Development Cooperation and Humanitarian Aid (DGD FA5 project), the Culmen International LLC,
the US CDC Atlanta, and the Agence Française de Développement through the AFROSCREEN project (grant agreement
CZZ3209), coordinated by ANRS Maladies Infectieuses émergentes in partnership with Institut Pasteur and Institut de
Recherche pour le Développement, the French Ministry of Europe and Foreign Affairs through FEF programme and support
provided by Institut de Recherche pour le Développement. AREBO project funded by ANRS-MIE. A.A.-A is supported by a
DGD sandwich PhD scholarship. P.M-K. acknowledges the support of the Wellcome Trust through the ARTIC Network (award
313694/Z/24/Z) and the Gates Foundation. This work was supported by Anthropic, which provided the credits and compute
used for the analyses. Claude Science carried out the phylodynamic runs, the estimation of reproduction numbers, and
the production of figures and tables under the authors’ direction, and the analyses reported here would not have been
completed in this time frame without it.

**Competing interests.** L.L. is an employee of Anthropic and holds equity in the company. Anthropic develops Claude
Science, the software used in this analysis. L.L. contributed to the analysis, interpretation, and writing of this
post; scientific outputs
were reviewed by all authors. The remaining authors declare no competing interests.

**Statement on continuing work and analyses prior to publication.** Please note that this data is based on work in
progress and should be considered preliminary. Our analyses are ongoing, and a publication communicating our findings
is in preparation. Sequences are publicly accessible under the Pathoplexus 'Restricted' licence and we would be
grateful if the terms of this were respected. If you intend to use our data prior to our publication, please contact
Dr Tony Wawina-Bokalanga (INRB, DRC) and/or Prof. Placide
Mbala-Kingebeni (INRB, DRC).
