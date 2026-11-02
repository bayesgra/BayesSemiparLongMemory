# Supplied data

These CSV files are copied byte-for-byte from the user's uploaded `Github.zip`:

* `SP500.csv`: daily S&P 500 price observations, dated April 22, 2019–April 19, 2024. Dots indicating non-trading/missing observations are removed before computing consecutive observed-price log returns. The upload does not include a redistribution licence.
* `HadCRUT.5.0.2.0.analysis.summary_series.global.monthly.csv`: the supplied monthly global temperature anomaly record. The workflow uses first differences of the anomaly column, as in the old script. The upload does not include source licence metadata.
* `global_temp.csv`: a separate supplied annual record, preserved for provenance but not used by the analysis runner.

The Nile observations were not supplied as a CSV. The old script uses `longmemo::NileMin`, which is available if the optional R package `longmemo` is installed. The runner deliberately does not substitute R's built-in `Nile` data, which are a different record.

The complete manuscript supplement was not supplied. Nile and temperature preprocessing follows the legacy code and is identified as such, rather than claimed to reproduce unavailable supplementary analyses. Before public redistribution, retain or obtain the original data source and licensing attribution applicable to these datasets.
