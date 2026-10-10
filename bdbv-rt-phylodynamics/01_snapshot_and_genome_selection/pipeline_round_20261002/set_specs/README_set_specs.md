# Specifications of the genome sets C95, C90, C80 and C00 (round of October 2026)

Written on 2 October 2026 by `round_20261002/make_completeness_specs.py`, before the new snapshot existed
(analysis plan of 2 October 2026, section 2). Format: the set-specification format of `build_alignments.py --set-spec`.

| file | set | key `min_completeness` | called positions needed (of 18,900) |
|---|---|---|---|
| `alignment_set_spec_C95_20261002.json` | C95 | 0.95 | 17,955 |
| `alignment_set_spec_C90_20261002.json` | C90 | 0.9 | 17,010 |
| `alignment_set_spec_C80_20261002.json` | C80 | 0.8 | 15,120 |
| `alignment_set_spec_C00_20261002.json` | C00 | null | none |
| `alignment_set_spec_completeness_all_20261002.json` | the same four sets in one file, primary set C90 | | |

Every set has the filters of the set `primary` of the frozen specification of 24 Sep 2026
(`pipeline_20260924/alignment_set_spec_20260924.json`): `dated: true`, `exclude_classes: [E, E-, D, M, H]`,
`duplicate_rule: drop-discordant`, `data_use: exclude`. The sets differ in `min_completeness` only.

What is NOT in the specification and has to be given on the command line, as on 24 Sep:
`--completeness-precision exact` (the threshold acts on the number of called positions), `--duplicate-pairs FILE`
(table of the records that share a sample identifier), `--data-use-exclude LIST` (restricted-use genomes that are left out),
`--screen FILE` (table of `screen.py --rules v2`). `round_20261002/build_completeness_sets.py` gives them.

The top-level key `external_data_set` of the frozen specification is not copied: it describes the data set of another
analysis and enters no rule of membership. Without it the table of genome sets holds the yes/no column
`on_external_exclusion_table` in place of `state_in_external_analysis`.

One file per set: the pipeline writes the first reason for absence for the primary set of a specification only; a run with
`alignment_set_spec_<set>_20261002.json` gives the column `reason_not_in_<set>`.
