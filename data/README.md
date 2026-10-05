# Data provenance and replacement

`provisional_diameters.csv` contains the 52 measured `D_eq (nm)` values from the author-supplied `PrecipitationResults_Summary.xlsx`. Extraction uses Excel row 2 as the column header and rows 3 onward as measurements. Values are copied without filtering or rounding beyond the precision stored in the source.

| Aging time | Workbook particle sheet | Count | Report Table 4 count | Report mean (nm) |
|---|---|---:|---:|---:|
| 15 min | 15_min_0p25_h | 32 | 32 | 123.9 |
| 60 min | 60_min_1_h | 13 | 22 | 155.0 |
| 165 min | 165_min_2p75_h | 7 | 14 | 198.7 |

The workbook Summary lists means of 123.9, 188.3, and 221.1 nm. The raw measurements are authoritative for this runnable example; calculated statistics are saved alongside the plots. The example must not be represented as an exact reproduction of the final report's Figure 7 or Figure 8. The author plans to supply the final dataset.

`diameters_template.csv` is intentionally header-only. Add one particle per row, with `aging_time_min` and `diameter_nm`. Do not synthesize particle measurements from the report's summary statistics.

The script reads existing equivalent diameters; it does not segment the embedded TEM images. No outlier or morphology filtering is performed. In particular, all supplied measurements remain included even when they may represent irregular or merged particles.
