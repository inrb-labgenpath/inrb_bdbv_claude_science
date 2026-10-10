# Weekly confirmed cases by province and health zone, 2026 BDBV outbreak (public source)

State of the source: 1 October 2026. Tables built on 2 October 2026.

## Source

Situation reports of the Institut National de Sante Publique (INSP), Democratic Republic of the Congo (series SitRep MVE 2026), as transcribed in the public repository INRB-UMIE/BDBV2026-Data (https://github.com/INRB-UMIE/BDBV2026-Data), commit 8eb57154cf967e7dfdb7e690230f61f2a298d2b0 of 1 October 2026 (08:05 UTC), the most recent commit of the branch main when the files were retrieved on 2 October 2026 (11:14 to 11:15 UTC). Its message is "CI: release build-2026-10-01-b077a6c"; the commit "add sitRep 138" (352cf0e8) precedes it on 1 October 2026.

Files used:

- `data/insp_sitrep/processed/insp_sitrep__cumulative_confirmed_cases__daily.csv`: cumulative confirmed cases by health zone and report date (3,995 rows, 98 report dates);
- `data/insp_sitrep/processed/insp_sitrep__national_cumulative_confirmed_cases__daily.csv`: national cumulative confirmed cases (103 report dates);
- `data/shapefiles/DRC_Health_zones.dbf`: health zones with their provinces (519 health zones).

The SHA-256 of the three files is recorded in `case_series_checks_20261002.json`. The Git hash of every supplied file of the commit (39 files) equals the hash in the list of files of the commit.

## Most recent situation report

The most recent report date in both tables is **29 September 2026** (first report date: 14 May 2026). The folder of situation reports of the repository (`data/insp_sitrep/raw/`) lists 102 PDF files with 101 distinct report numbers between 001 and 138 (report 012 in two versions). The reports themselves were not among the files of this analysis and were not read; the link between report numbers and report dates was not verified.

## Licence and terms of use

- Repository licence: MIT (`LICENSE.md`: "Copyright (c) 2026 Kraemer Lab, University of Oxford").
- README of the repository: "The repository code is licensed under the terms in LICENSE. We do not claim ownership of or the right to license the third-party data or software tools used. Please pass forward any existing license/warranty/copyright information when redistributing."
- Metadata of the folder of situation reports (`data/insp_sitrep/metadata.yaml`): "INSP outbreak sitreps; reuse with attribution to INSP and citation of the specific report number and date. Confirm distribution terms with INSP before external republication."
- The README of the repository describes the epidemiological data as preliminary and work in progress.

## Processing (this work)

1. **Health zones and provinces.** Health zones were assigned to provinces with the health-zone table of the same repository. 67 of the 71 names of the source have exactly one province there; the other 4 were assigned by a fixed list: Gethy (not in the table under this spelling) to Ituri, as Gety; Lubunga (2 provinces in the table: Kasaï-Central, Tshopo) to Tshopo; Lubunga (Tshopo) (2 provinces in the table: Kasaï-Central, Tshopo) to Tshopo; Rumba (not in the table under this spelling) to Ituri, as Rimba. 6 spellings of a health zone were merged with another spelling (Gethy = Gety; Lubunga (Tshopo) = Lubunga; Makiso Kisangani = Makiso-Kisangani; Miti Murhesa = Miti-Murhesa; Nia Nia = Nia-Nia; Rumba = Rimba). After merging there are 65 health zones: 29 in Ituri, 17 in Nord-Kivu and 19 in 5 other provinces (Bas-Uele, Haut-Uele, Sud-Kivu, Sud-Ubangi, Tshopo). See `zone_province_map_20261002.csv`.
2. **Values repaired or set aside.** 1 date string with a stray character (`2026-06-25]`, Nyankunde) was read as 25 June 2026 and the row was used. 11 values that are not numbers ('ND' in 10 rows, '17-' in 1 row) were set aside. 42 rows whose name of the health zone is `NA` (cases that the source does not assign to a health zone; report dates 1 Jun to 20 Jul) were not attributed to a province; these cases enter through the national count only. 1 health zone (Jiba) has no numeric value. See `case_series_values_set_aside_20261002.csv`.
3. **One value per health zone and report date; carrying forward.** Where two spellings of a health zone report on the same date the larger value is used. The cumulative count of a health zone that is absent from a report is carried forward from its last report (0 before its first report). On 6 report dates of the national table (17 May, 31 May, 27 Jun, 8 Aug, 9 Aug, 10 Aug) the table of the health zones has no row, and all counts are carried forward. 1 report date of the table of the health zones (15 May) has no national count and is not a date of the series.
4. **Decreases are kept.** The cumulative count of a health zone falls from one report to the next 23 times in 15 health zones; the values were used as transcribed. The national cumulative count never falls.
5. **Weekly counts.** Cumulative counts by province were interpolated linearly between report dates (103 report dates, 14 May to 29 September 2026) and differenced at Sundays. This gives 19 complete weeks, weeks ending 24 May to 27 September 2026. The report of 29 Sep lies after the last Sunday and is used for the interpolation only. Weekly counts are rounded to one decimal.
6. **Cases attributed to health zones with delay.** Cases not attributed to a health zone at the time of a report (national count minus sum of the health zones) appear with delay in the counts of the health zones. Consecutive weeks were pooled until the pooled unattributed count was within 5 % of the pooled national count; within a pool the national weekly count is split by province in proportion to the pooled increments attributed to health zones (column `weeks_pooled_for_province_split`). Pools of more than one week: weeks ending 31 May to 21 Jun (4 weeks); weeks ending 28 Jun to 5 Jul (2 weeks); weeks ending 9 Aug to 16 Aug (2 weeks). National weekly totals are not affected. Because the split uses counts rounded to one decimal, a province column can differ by 0.1 from the attributed count of a week that is not pooled.
7. **Set of health zones.** The set of health zones is read from `zones_of_analysed_set_20261002.txt`: Bunia, Rwampara, Nizi, Mongbwalu, Lita, Mangala, Bambu (7 health zones, all in Ituri: yes). Its weekly count is given as attributed and with three treatments of the pooled weeks: A, the national weekly count multiplied by the share of the set among the attributed cases of the pool (the treatment of the provinces); B, every week of a pool replaced by the mean of the pool; M, the weeks ending 9 and 16 August replaced by their mean.
8. **Four-week periods.** Mean weekly counts are given for consecutive periods of four weeks counted from 25 May 2026; complete periods only (4 periods, to 13 September 2026). 2 complete weeks lie after the last complete period.
9. **Growth and reproduction number.** The growth rate r is the slope of an ordinary least-squares regression of the logarithm of the weekly count on time (weeks with a count of 0 or less would be left out; there was none). R = (1 + r s^2/m)^(m^2/s^2) for a gamma-distributed generation time with mean m = 15.3 days and standard deviation s = 9.3 days. The interval is the 95 % confidence interval of the slope, converted in the same way. It reflects the scatter of the weekly counts about the fitted line only; changes in testing, ascertainment and reporting are not modelled.
10. **Weeks after 13 September.** A range of weeks ending after 13 September 2026 is evaluated only if it holds at least 3 complete weeks. It holds 2 (20 Sep, 27 Sep); no growth rate and no reproduction number are given for it.
11. **Weeks that rest on few reports.** The week ending 27 Sep rests on 2 report dates (21 Sep, 26 Sep), fewer than any other week, and its last day lies between the reports of 26 Sep and 29 Sep; it is to be marked as provisional in figures. The week ending 20 Sep rests on 3 report dates.

## Cross-check within the repository

The national count was compared with the sum of the health zones at each of the 103 report dates (`case_series_crosscheck_20261002.csv`). The two are equal on 40 report dates and differ on 63, the last time on 17 August 2026; on all 26 report dates after it (19 Aug to 29 September 2026) they are equal. The largest differences are 476 cases on 10 Aug (national count above the health zones) and 33 cases on 28 May (health zones above the national count). The 19 weekly national counts add up to 8,106.3; the national cumulative count rises by 8,106.3 between the first and the last Sunday.

A second transcription of the situation reports was not among the files of this analysis; the counts were not compared with one, and they were not compared with the situation reports themselves.

Counts are by date of report, not by date of onset or infection.

## Files

- `weekly_confirmed_cases_by_province_20261002.csv`: weekly confirmed cases, national and by province (Ituri, Nord-Kivu, other provinces), weeks ending Sunday; columns raw_*: counts as attributed to health zones; weeks_pooled_for_province_split: size of the pool of weeks; report_dates_in_week.
- `cumulative_confirmed_cases_by_province_20261002.csv`: cumulative confirmed cases at every report date: sum of the health zones of Ituri, of Nord-Kivu and of the other provinces, sum of all health zones, national count, national minus health zones.
- `zone_province_map_20261002.csv`: every name of the source with its health zone and province.
- `mean_weekly_cases_by_zone_group_20261002.csv`: mean weekly cases by four-week period and group of health zones, counts as attributed.
- `mean_weekly_cases_seven_zones_20261002.csv`: mean weekly cases by four-week period for each health zone of the set, counts as attributed.
- `mean_weekly_cases_seven_zones_by_treatment_20261002.csv`: mean weekly cases of the set of health zones by four-week period: as attributed, treatment A, treatment B.
- `weekly_confirmed_cases_seven_zones_20261002.csv`: weekly cases of the set of health zones: as attributed and with the treatments A, B and M; pools of weeks.
- `case_growth_comparators_20261002.csv`: growth rate and reproduction number from weekly counts: DRC, Ituri, Nord-Kivu, the set of health zones, the other health zones of Ituri.
- `case_growth_seven_zones_treatments_20261002.csv`: the same for the set of health zones with each of the four treatments, full precision.
- `case_series_crosscheck_20261002.csv`: national count and sum of the health zones at every report date.
- `case_series_values_set_aside_20261002.csv`: every value of the source that was repaired, set aside or not attributed to a province.
- `case_series_checks_20261002.json`: counts and checks of the run, with the SHA-256 of the three input files.
- Code: `build_case_series_20261002.py` (reads the three files, writes all tables; no network access); `zones_of_analysed_set_20261002.txt` (the set of health zones).

Run: `python build_case_series_20261002.py --indir <folder with the files of the commit> --outdir <folder> --suffix 20261002 --zones-file zones_of_analysed_set_20261002.txt`. To use another set of health zones, replace the names in the text file; if the number of zones changes, pass `--zone-set-label` and `--zone-set-tag` as well.
