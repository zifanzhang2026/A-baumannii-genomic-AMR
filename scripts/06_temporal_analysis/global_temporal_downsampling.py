#!/usr/bin/env python3
"""
Global equal-n temporal downsampling sensitivity analysis
=========================================================

Edit ONLY the USER SETTINGS section below, then run:

    python 02_global_temporal_downsampling_easy.py

Purpose
-------
Assess whether the global temporal pattern in acquired ARG burden is robust
to unequal numbers of genomes across collection periods.

This sensitivity analysis addresses unequal sample size across periods.
It does NOT correct geographic, surveillance, lineage, or host-composition bias.
"""

# ============================================================
# USER SETTINGS — EDIT THIS SECTION ONLY
# ============================================================

# 1. Final dataset used for the analysis (CSV or TSV).
#    Example:
# INPUT_FILE = "/public/home/yourname/project/final_19861_genomes.csv"
INPUT_FILE = "/path/to/your/final_dataset.csv"

# 2. Column containing collection year.
YEAR_COL = "year"

# 3. Column containing the acquired ARG count/burden for each genome.
#    Keep "NUM_FOUND" if this is the column used in your final dataset.
ARG_COUNT_COL = "NUM_FOUND"

# 4. Output directory.
#    Example:
# OUTPUT_DIR = "/public/home/yourname/project/downsampling_global"
OUTPUT_DIR = "./downsampling_global"

# 5. Number of repeated equal-n downsampling iterations.
N_RESAMPLES = 1000

# 6. Fixed random seed for reproducibility.
RANDOM_SEED = 20260916

# Fixed global periods used in the manuscript sensitivity analysis:
# <=2000 | 2001-2009 | 2010-2019 | >=2020
#
# Normally, DO NOT change these boundaries when reproducing the analysis.


# ============================================================
# ANALYSIS CODE — NO NEED TO EDIT BELOW THIS LINE
# ============================================================

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal


PERIOD_ORDER = ["<=2000", "2001-2009", "2010-2019", ">=2020"]


def read_table(path):
    """Read CSV or tab-delimited input based on filename extension."""
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"\nInput file was not found:\n{path}\n\n"
            "Please edit INPUT_FILE in the USER SETTINGS section."
        )

    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def assign_period(year):
    """Assign collection year to one of four fixed global periods."""
    if pd.isna(year):
        return pd.NA

    y = int(year)

    if y <= 2000:
        return "<=2000"
    elif y <= 2009:
        return "2001-2009"
    elif y <= 2019:
        return "2010-2019"
    else:
        return ">=2020"


def main():

    # --------------------------------------------------------
    # 1. Basic settings checks
    # --------------------------------------------------------

    if INPUT_FILE == "/path/to/your/final_dataset.csv":
        raise ValueError(
            "\nYou have not set INPUT_FILE yet.\n"
            "Open this script and replace:\n\n"
            'INPUT_FILE = "/path/to/your/final_dataset.csv"\n\n'
            "with the path to your final dataset."
        )

    if N_RESAMPLES < 1:
        raise ValueError("N_RESAMPLES must be >= 1.")

    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # 2. Read final dataset
    # --------------------------------------------------------

    df = read_table(INPUT_FILE)

    print("=" * 70)
    print("Global temporal downsampling sensitivity analysis")
    print("=" * 70)
    print(f"Input file:       {INPUT_FILE}")
    print(f"Input rows:       {len(df):,}")
    print(f"Year column:      {YEAR_COL}")
    print(f"ARG-count column: {ARG_COUNT_COL}")
    print()

    required_columns = [YEAR_COL, ARG_COUNT_COL]
    missing_columns = [c for c in required_columns if c not in df.columns]

    if missing_columns:
        raise ValueError(
            "Required column(s) not found in the input file: "
            + ", ".join(missing_columns)
            + "\n\nAvailable columns are:\n"
            + "\n".join(map(str, df.columns.tolist()))
        )

    # --------------------------------------------------------
    # 3. Prepare year and acquired-ARG burden
    # --------------------------------------------------------

    dat = df.copy()

    dat["_year"] = pd.to_numeric(dat[YEAR_COL], errors="coerce")
    dat["_arg_count"] = pd.to_numeric(dat[ARG_COUNT_COL], errors="coerce")

    n_input = len(dat)
    missing_year = int(dat["_year"].isna().sum())
    missing_arg = int(dat["_arg_count"].isna().sum())

    # Complete cases required only for year and ARG count.
    dat = dat.loc[
        dat["_year"].notna() &
        dat["_arg_count"].notna()
    ].copy()

    if len(dat) == 0:
        raise ValueError("No usable rows remain after checking year and ARG count.")

    if (dat["_arg_count"] < 0).any():
        raise ValueError(
            "Negative values were detected in the acquired ARG count column."
        )

    # Collection year should be integer-valued.
    non_integer_year = ~np.isclose(dat["_year"], np.round(dat["_year"]))

    if non_integer_year.any():
        examples = dat.loc[
            non_integer_year, YEAR_COL
        ].head(10).tolist()

        raise ValueError(
            "Non-integer collection years were detected. "
            f"Examples: {examples}"
        )

    dat["_year"] = dat["_year"].astype(int)

    # --------------------------------------------------------
    # 4. Assign fixed global periods
    # --------------------------------------------------------

    dat["period"] = dat["_year"].map(assign_period)

    period_counts = (
        dat["period"]
        .value_counts()
        .reindex(PERIOD_ORDER, fill_value=0)
    )

    if (period_counts == 0).any():
        absent_periods = period_counts[
            period_counts == 0
        ].index.tolist()

        raise ValueError(
            "At least one global period contains no usable genomes: "
            + ", ".join(absent_periods)
        )

    # Equal sample size = size of smallest period.
    n_min = int(period_counts.min())

    print("Available genomes by period:")
    for period in PERIOD_ORDER:
        print(f"  {period:10s}: {int(period_counts[period]):,}")

    print()
    print(f"Equal sample size per period per iteration: {n_min:,}")
    print(f"Number of iterations: {N_RESAMPLES:,}")
    print(f"Random seed: {RANDOM_SEED}")
    print()

    period_sample_sizes = pd.DataFrame({
        "period": PERIOD_ORDER,
        "available_n": [
            int(period_counts[p]) for p in PERIOD_ORDER
        ],
        "sampled_n_per_iteration": n_min,
    })

    period_sample_sizes.to_csv(
        output_dir / "period_sample_sizes.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # 5. Prepare period-specific arrays
    # --------------------------------------------------------

    values = {
        period: dat.loc[
            dat["period"] == period,
            "_arg_count"
        ].to_numpy(dtype=float)
        for period in PERIOD_ORDER
    }

    rng = np.random.default_rng(RANDOM_SEED)

    iteration_mean_rows = []
    kw_rows = []

    # --------------------------------------------------------
    # 6. Repeated equal-n downsampling
    # --------------------------------------------------------

    for iteration in range(1, N_RESAMPLES + 1):

        sampled = {}
        mean_row = {"iteration": iteration}

        for period in PERIOD_ORDER:

            indices = rng.choice(
                len(values[period]),
                size=n_min,
                replace=False,
            )

            sampled_values = values[period][indices]

            sampled[period] = sampled_values
            mean_row[period] = float(
                np.mean(sampled_values)
            )

        iteration_mean_rows.append(mean_row)

        # Secondary robustness statistic only.
        kw_result = kruskal(
            *(sampled[p] for p in PERIOD_ORDER)
        )

        kw_rows.append({
            "iteration": iteration,
            "H": float(kw_result.statistic),
            "p_value": float(kw_result.pvalue),
        })

    # --------------------------------------------------------
    # 7. Save iteration-level means
    # --------------------------------------------------------

    iteration_means = pd.DataFrame(
        iteration_mean_rows
    )

    iteration_means.to_csv(
        output_dir / "downsampling_iteration_means.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # 8. Summarize resampling distribution
    # --------------------------------------------------------

    summary_rows = []

    for period in PERIOD_ORDER:

        x = iteration_means[
            period
        ].to_numpy(dtype=float)

        summary_rows.append({
            "period": period,
            "available_n": int(period_counts[period]),
            "sampled_n_per_iteration": n_min,
            "n_resamples": N_RESAMPLES,
            "mean_of_resampled_means": float(np.mean(x)),
            "median_resampled_mean": float(np.median(x)),
            "percentile_2.5": float(np.percentile(x, 2.5)),
            "percentile_97.5": float(np.percentile(x, 97.5)),
        })

    downsampling_summary = pd.DataFrame(
        summary_rows
    )

    downsampling_summary.to_csv(
        output_dir / "downsampling_summary.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # 9. Save secondary Kruskal-Wallis results
    # --------------------------------------------------------

    kw_results = pd.DataFrame(kw_rows)

    kw_results.to_csv(
        output_dir / "downsampling_kruskal_wallis.tsv",
        sep="\t",
        index=False,
    )

    proportion_kw_p_lt_005 = float(
        (kw_results["p_value"] < 0.05).mean()
    )

    # --------------------------------------------------------
    # 10. Human-readable analysis summary
    # --------------------------------------------------------

    summary_file = output_dir / "analysis_summary.txt"

    with open(
        summary_file,
        "w",
        encoding="utf-8",
    ) as f:

        f.write(
            "Global equal-n temporal downsampling sensitivity analysis\n"
        )
        f.write(
            "========================================================\n"
        )

        f.write(f"Input file:\t{INPUT_FILE}\n")
        f.write(f"Input rows:\t{n_input}\n")
        f.write(
            f"Rows with missing/unparseable year:\t"
            f"{missing_year}\n"
        )
        f.write(
            f"Rows with missing/unparseable ARG count:\t"
            f"{missing_arg}\n"
        )
        f.write(f"Rows used:\t{len(dat)}\n")

        f.write(
            f"Number of resamples:\t{N_RESAMPLES}\n"
        )
        f.write(
            f"Random seed:\t{RANDOM_SEED}\n"
        )

        f.write(
            "Sampling:\tequal n across periods; "
            "without replacement within each iteration\n"
        )

        f.write(
            f"n per period per iteration:\t{n_min}\n"
        )

        f.write(
            "Periods:\t"
            "<=2000 | 2001-2009 | 2010-2019 | >=2020\n"
        )

        f.write("\nAvailable sample sizes:\n")

        for period in PERIOD_ORDER:
            f.write(
                f"{period}\t"
                f"{int(period_counts[period])}\n"
            )

        f.write(
            "\nSecondary robustness statistic:\n"
        )

        f.write(
            "Kruskal-Wallis p<0.05 proportion "
            "across resamples:\t"
            f"{proportion_kw_p_lt_005:.6f}\n"
        )

        f.write(
            "\nInterpretation note:\t"
            "This downsampling analysis addresses unequal "
            "sample size across temporal periods; it does not "
            "correct surveillance, geographic, lineage, or "
            "host-composition bias.\n"
        )

        f.write(
            "Inference note:\t"
            "The distribution of resampled mean ARG burdens "
            "is the primary downsampling output. "
            "Mean p-values across resamples are not used.\n"
        )

    # --------------------------------------------------------
    # 11. Print final results
    # --------------------------------------------------------

    print("=" * 70)
    print("Downsampling completed successfully.")
    print("=" * 70)

    print("\nSummary of resampled mean acquired ARG burden:\n")

    print(
        downsampling_summary.to_string(
            index=False
        )
    )

    print(
        "\nProportion of iterations with "
        "Kruskal-Wallis p < 0.05: "
        f"{proportion_kw_p_lt_005:.3f}"
    )

    print(
        f"\nOutput directory:\n{output_dir.resolve()}"
    )

    print(
        "\nIMPORTANT:"
        "\nThis analysis is a sensitivity analysis for unequal "
        "sample size across periods."
        "\nDo not describe it as correction for surveillance "
        "or compositional bias."
    )


if __name__ == "__main__":
    main()
