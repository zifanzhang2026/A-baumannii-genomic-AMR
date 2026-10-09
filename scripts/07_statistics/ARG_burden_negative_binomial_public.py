# ============================================================
# Acquired ARG burden analysis
#
# Outcome:
#   NUM_FOUND = number of acquired resistance genes per genome
#
# Model:
#   NUM_FOUND ~ ST_group + host_group + continent + year
#
# Main outputs:
#   IRR
#   95% CI
#   P value
#   global P values
#   Poisson overdispersion diagnostics
#   Negative Binomial model
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import statsmodels.api as sm
import statsmodels.formula.api as smf

from scipy.stats import chi2
from patsy.contrasts import Treatment


# ============================================================
# 1. Input files
# ============================================================

# Original acquired ARG file
# The second column NUM_FOUND is the total number of
# acquired resistance genes detected in each genome.
arg_file = "acquired_resistance_genes.csv"

# Previously generated master dataset
merged_file = "merged_ARG_metadata.csv"


# Output directory
output_dir = Path("ARG_count_model_results")
output_dir.mkdir(exist_ok=True)


# ============================================================
# 2. Analysis settings
# ============================================================

# Reference groups
ST_REFERENCE = "Others"
HOST_REFERENCE = "Human"
CONTINENT_REFERENCE = "North America"

# Year scaling:
# IRR will be interpreted per 5-year increase.
YEAR_CENTER = 2015
YEAR_SCALE = 5


# Six valid continents
valid_continents = [
    "Africa",
    "Asia",
    "Europe",
    "North America",
    "South America",
    "Oceania"
]


# A pragmatic threshold used only to describe the Poisson overdispersion diagnostic.
# It does NOT determine the manuscript primary model.
OVERDISPERSION_THRESHOLD = 1.5


# ============================================================
# 3. Function for reading CSV files
# ============================================================

def read_csv_auto(file_path):

    encodings = [
        "utf-8-sig",
        "utf-8",
        "gb18030",
        "gbk"
    ]

    for encoding in encodings:

        try:

            df = pd.read_csv(
                file_path,
                encoding=encoding,
                low_memory=False
            )

            print(
                f"Successfully read: {file_path}\n"
                f"Encoding: {encoding}\n"
                f"Rows: {df.shape[0]:,}\n"
                f"Columns: {df.shape[1]:,}\n"
            )

            return df

        except UnicodeDecodeError:
            continue

    raise ValueError(
        f"Unable to determine encoding for: {file_path}"
    )


# ============================================================
# 4. Read files
# ============================================================

arg = read_csv_auto(arg_file)
merged = read_csv_auto(merged_file)


# ============================================================
# 5. Clean column names
# ============================================================

arg.columns = (
    arg.columns
    .astype(str)
    .str.strip()
)

merged.columns = (
    merged.columns
    .astype(str)
    .str.strip()
)


# ============================================================
# 6. Check required columns
# ============================================================

required_arg_columns = [
    "Genome_ID",
    "NUM_FOUND"
]

required_merged_columns = [
    "Genome_ID",
    "ST_group",
    "host_group",
    "continent",
    "year"
]


missing_arg_columns = [
    col
    for col in required_arg_columns
    if col not in arg.columns
]

missing_merged_columns = [
    col
    for col in required_merged_columns
    if col not in merged.columns
]


if missing_arg_columns:

    raise ValueError(
        "Missing columns in acquired_resistance_genes.csv:\n"
        + "\n".join(missing_arg_columns)
    )


if missing_merged_columns:

    raise ValueError(
        "Missing columns in merged_ARG_metadata.csv:\n"
        + "\n".join(missing_merged_columns)
    )


# ============================================================
# 7. Clean Genome_ID
# ============================================================

arg["Genome_ID"] = (
    arg["Genome_ID"]
    .astype("string")
    .str.strip()
)

merged["Genome_ID"] = (
    merged["Genome_ID"]
    .astype("string")
    .str.strip()
)


# ============================================================
# 8. Check duplicate Genome_ID
# ============================================================

def check_duplicate_id(df, name):

    duplicates = df.loc[
        df["Genome_ID"].duplicated(
            keep=False
        ),
        "Genome_ID"
    ]

    if len(duplicates) > 0:

        print(
            f"\nDuplicate Genome_ID found in {name}:"
        )

        print(
            duplicates
            .value_counts()
            .head(20)
        )

        raise ValueError(
            f"Genome_ID is not unique in {name}."
        )

    else:

        print(
            f"{name}: "
            f"{df['Genome_ID'].nunique():,} unique Genome_IDs"
        )


check_duplicate_id(
    arg,
    "acquired_resistance_genes.csv"
)

check_duplicate_id(
    merged,
    "merged_ARG_metadata.csv"
)


# ============================================================
# 9. Extract original NUM_FOUND
# ============================================================

arg_count = arg[
    [
        "Genome_ID",
        "NUM_FOUND"
    ]
].copy()


# Convert NUM_FOUND to numeric
arg_count["NUM_FOUND"] = pd.to_numeric(
    arg_count["NUM_FOUND"],
    errors="coerce"
)


# ============================================================
# 10. Validate NUM_FOUND
# ============================================================

print("\n========================================")
print("Checking NUM_FOUND")
print("========================================")


missing_num_found = (
    arg_count["NUM_FOUND"]
    .isna()
    .sum()
)

print(
    f"Missing NUM_FOUND: "
    f"{missing_num_found:,}"
)


# Negative counts should not exist
negative_counts = (
    arg_count["NUM_FOUND"] < 0
).sum()

print(
    f"Negative NUM_FOUND values: "
    f"{negative_counts:,}"
)


if negative_counts > 0:

    raise ValueError(
        "NUM_FOUND contains negative values."
    )


# Check whether counts are integers
non_integer = arg_count.loc[
    arg_count["NUM_FOUND"].notna()
    &
    (
        np.abs(
            arg_count["NUM_FOUND"]
            -
            np.round(
                arg_count["NUM_FOUND"]
            )
        )
        >
        1e-8
    )
]


if len(non_integer) > 0:

    raise ValueError(
        "NUM_FOUND contains non-integer values."
    )


# Convert to nullable integer
arg_count["NUM_FOUND"] = (
    arg_count["NUM_FOUND"]
    .round()
    .astype("Int64")
)


# ============================================================
# 11. Remove NUM_FOUND if it already exists in merged dataset
#
# We deliberately use NUM_FOUND from the ORIGINAL
# acquired_resistance_genes.csv.
# ============================================================

if "NUM_FOUND" in merged.columns:

    print(
        "\nNUM_FOUND already exists in merged dataset."
        "\nIt will be removed and replaced with the "
        "original NUM_FOUND."
    )

    merged = merged.drop(
        columns=["NUM_FOUND"]
    )


# ============================================================
# 12. Merge NUM_FOUND into master dataset
# ============================================================

data = pd.merge(
    merged,
    arg_count,
    on="Genome_ID",
    how="left",
    validate="one_to_one",
    indicator="_NUM_FOUND_merge"
)


print("\n========================================")
print("NUM_FOUND merge status")
print("========================================")

print(
    data["_NUM_FOUND_merge"]
    .value_counts(
        dropna=False
    )
)


n_missing_after_merge = (
    data["NUM_FOUND"]
    .isna()
    .sum()
)


print(
    f"\nGenomes without NUM_FOUND after merge: "
    f"{n_missing_after_merge:,}"
)


# IMPORTANT:
# Do NOT convert missing NUM_FOUND to 0.
# Missing means no matching information, not absence of ARGs.


# ============================================================
# 13. Save master dataset with NUM_FOUND
# ============================================================

master_output = (
    output_dir
    /
    "merged_ARG_metadata_with_NUM_FOUND.csv"
)


data.to_csv(
    master_output,
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 14. Overall descriptive statistics
#
# Use all genomes with valid NUM_FOUND,
# not only regression-eligible genomes.
# ============================================================

all_counts = (
    data["NUM_FOUND"]
    .dropna()
    .astype(float)
)


n_total_count = len(all_counts)

mean_count = all_counts.mean()
variance_count = all_counts.var(ddof=1)
sd_count = all_counts.std(ddof=1)

median_count = all_counts.median()

q1 = all_counts.quantile(0.25)
q3 = all_counts.quantile(0.75)

minimum = all_counts.min()
maximum = all_counts.max()

zero_n = int(
    (all_counts == 0).sum()
)

zero_percent = (
    zero_n
    /
    n_total_count
    *
    100
)


variance_mean_ratio = (
    variance_count
    /
    mean_count
    if mean_count > 0
    else np.nan
)


descriptive_summary = pd.DataFrame(
    {
        "Statistic": [
            "N",
            "Mean",
            "Variance",
            "SD",
            "Median",
            "Q1",
            "Q3",
            "Minimum",
            "Maximum",
            "Zero_count",
            "Zero_percent",
            "Variance_to_mean_ratio"
        ],
        "Value": [
            n_total_count,
            mean_count,
            variance_count,
            sd_count,
            median_count,
            q1,
            q3,
            minimum,
            maximum,
            zero_n,
            zero_percent,
            variance_mean_ratio
        ]
    }
)


descriptive_summary.to_csv(
    output_dir
    /
    "01_NUM_FOUND_descriptive_summary.csv",
    index=False,
    encoding="utf-8-sig"
)


print("\n========================================")
print("Overall NUM_FOUND distribution")
print("========================================")

print(
    descriptive_summary.to_string(
        index=False
    )
)


# ============================================================
# 15. Histogram of NUM_FOUND
# ============================================================

plt.figure(
    figsize=(8, 5)
)


# Integer-centered bins
hist_min = int(all_counts.min())
hist_max = int(all_counts.max())

bins = np.arange(
    hist_min - 0.5,
    hist_max + 1.5,
    1
)


plt.hist(
    all_counts,
    bins=bins
)

plt.xlabel(
    "Number of acquired resistance genes per genome"
)

plt.ylabel(
    "Number of genomes"
)

plt.title(
    "Distribution of acquired ARG counts"
)

plt.tight_layout()


hist_file = (
    output_dir
    /
    "02_NUM_FOUND_distribution.png"
)


plt.savefig(
    hist_file,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 16. Prepare year
# ============================================================

data["year"] = pd.to_numeric(
    data["year"],
    errors="coerce"
)


# ============================================================
# 17. Create regression dataset
#
# Same complete-case principle as the previous
# logistic regression analysis.
#
# ST:
#   Others
#   ST2
#
# Host:
#   Human
#   Animal_Environment
#
# Continent:
#   6 continents
#
# Year:
#   valid
#
# NUM_FOUND:
#   valid
# ============================================================

analysis_df = data[
    data["NUM_FOUND"].notna()
    &
    data["ST_group"].isin(
        [
            "Others",
            "ST2"
        ]
    )
    &
    data["host_group"].isin(
        [
            "Human",
            "Animal_Environment"
        ]
    )
    &
    data["continent"].isin(
        valid_continents
    )
    &
    data["year"].notna()
].copy()


analysis_df["NUM_FOUND"] = (
    analysis_df["NUM_FOUND"]
    .astype(int)
)


# ============================================================
# 18. Create year5
#
# Example:
#
# 2015 -> 0
# 2020 -> 1
# 2010 -> -1
#
# IRR = effect per 5-year increase
# ============================================================

analysis_df["year5"] = (
    analysis_df["year"]
    -
    YEAR_CENTER
) / YEAR_SCALE


print("\n========================================")
print("Regression dataset")
print("========================================")

print(
    f"Total genomes in master dataset: "
    f"{len(data):,}"
)

print(
    f"Genomes included in count regression: "
    f"{len(analysis_df):,}"
)

print(
    f"Excluded from count regression: "
    f"{len(data) - len(analysis_df):,}"
)


# ============================================================
# 19. Check predictor distributions
# ============================================================

print("\nST_group:")
print(
    analysis_df["ST_group"]
    .value_counts()
)

print("\nhost_group:")
print(
    analysis_df["host_group"]
    .value_counts()
)

print("\ncontinent:")
print(
    analysis_df["continent"]
    .value_counts()
)


# ============================================================
# 20. Descriptive ARG-count statistics by groups
# ============================================================

def group_count_summary(
    df,
    variable
):

    temp = (
        df
        .groupby(
            variable,
            observed=True
        )["NUM_FOUND"]
        .agg(
            N="count",
            Mean="mean",
            SD="std",
            Median="median",
            Minimum="min",
            Maximum="max"
        )
        .reset_index()
    )

    q1 = (
        df
        .groupby(
            variable,
            observed=True
        )["NUM_FOUND"]
        .quantile(0.25)
        .reset_index(
            name="Q1"
        )
    )

    q3 = (
        df
        .groupby(
            variable,
            observed=True
        )["NUM_FOUND"]
        .quantile(0.75)
        .reset_index(
            name="Q3"
        )
    )

    temp = temp.merge(
        q1,
        on=variable
    )

    temp = temp.merge(
        q3,
        on=variable
    )

    temp.insert(
        0,
        "Variable",
        variable
    )

    temp = temp.rename(
        columns={
            variable:
            "Level"
        }
    )

    return temp


group_summaries = []


for variable in [
    "ST_group",
    "host_group",
    "continent"
]:

    group_summaries.append(
        group_count_summary(
            analysis_df,
            variable
        )
    )


group_summary_df = pd.concat(
    group_summaries,
    ignore_index=True
)


group_summary_df.to_csv(
    output_dir
    /
    "03_NUM_FOUND_group_summary.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 21. Model formula
#
# Reference:
#
# ST:
#   Others
#
# Host:
#   Human
#
# Continent:
#   North America
#
# Year:
#   per 5-year increase
# ============================================================

formula = f"""
NUM_FOUND
~
C(
    ST_group,
    Treatment(reference="{ST_REFERENCE}")
)
+
C(
    host_group,
    Treatment(reference="{HOST_REFERENCE}")
)
+
C(
    continent,
    Treatment(reference="{CONTINENT_REFERENCE}")
)
+
year5
"""


# Remove line breaks for cleaner display
formula = " ".join(
    formula.split()
)


print("\n========================================")
print("Regression formula")
print("========================================")

print(formula)


# ============================================================
# 22. Fit Poisson regression
# ============================================================

poisson_model = smf.glm(
    formula=formula,
    data=analysis_df,
    family=sm.families.Poisson()
)


poisson_result = poisson_model.fit()


print("\n========================================")
print("Poisson regression completed")
print("========================================")


# ============================================================
# 23. Poisson overdispersion diagnostic
#
# Pearson dispersion =
#
# sum(Pearson residual^2) / residual df
#
# Approximately 1:
#   consistent with Poisson assumption
#
# Clearly > 1:
#   overdispersion
# ============================================================

pearson_chi2 = np.sum(
    poisson_result.resid_pearson ** 2
)

poisson_df_resid = (
    poisson_result.df_resid
)

poisson_dispersion = (
    pearson_chi2
    /
    poisson_df_resid
)


print(
    f"Poisson Pearson chi-square: "
    f"{pearson_chi2:.3f}"
)

print(
    f"Poisson residual df: "
    f"{poisson_df_resid:.0f}"
)

print(
    f"Poisson dispersion statistic: "
    f"{poisson_dispersion:.3f}"
)


if poisson_dispersion > OVERDISPERSION_THRESHOLD:

    print(
        "\nClear overdispersion detected."
    )

    print(
        "Negative Binomial regression is preferred."
    )

else:

    print(
        "\nNo strong overdispersion detected "
        "using the predefined threshold."
    )


# ============================================================
# 24. Fit Negative Binomial regression
#
# NB2:
#
# Var(Y) = mu + alpha * mu^2
#
# alpha is estimated from the data.
# ============================================================

nb_model = smf.negativebinomial(
    formula=formula,
    data=analysis_df,
    loglike_method="nb2"
)


nb_result = nb_model.fit(
    disp=False,
    maxiter=500
)


print("\n========================================")
print("Negative Binomial regression completed")
print("========================================")

print(
    f"Converged: "
    f"{nb_result.mle_retvals.get('converged', True)}"
)


# ============================================================
# 25. Extract estimated alpha
# ============================================================

if "alpha" in nb_result.params.index:

    alpha_estimate = (
        nb_result.params["alpha"]
    )

    alpha_pvalue = (
        nb_result.pvalues["alpha"]
    )

else:

    alpha_estimate = np.nan
    alpha_pvalue = np.nan


print(
    f"Estimated NB alpha: "
    f"{alpha_estimate:.6f}"
)

print(
    f"Alpha P value: "
    f"{alpha_pvalue:.6g}"
)


# ============================================================
# 26. Model comparison
# ============================================================

model_comparison = pd.DataFrame(
    {
        "Model": [
            "Poisson",
            "Negative Binomial"
        ],
        "Log_likelihood": [
            poisson_result.llf,
            nb_result.llf
        ],
        "AIC": [
            poisson_result.aic,
            nb_result.aic
        ],
        "BIC": [
            poisson_result.bic_llf,
            nb_result.bic
        ],
        "Poisson_dispersion": [
            poisson_dispersion,
            np.nan
        ],
        "NB_alpha": [
            np.nan,
            alpha_estimate
        ]
    }
)


model_comparison.to_csv(
    output_dir
    /
    "04_model_comparison.csv",
    index=False,
    encoding="utf-8-sig"
)


print("\n========================================")
print("Model comparison")
print("========================================")

print(
    model_comparison.to_string(
        index=False
    )
)


# ============================================================
# 27. Primary model used in the manuscript
#
# Poisson regression is retained only as an overdispersion
# diagnostic and for descriptive model comparison.
#
# The manuscript's primary count model is fixed as
# Negative Binomial (NB2), consistent with the reported methods.
# ============================================================

selected_model_name = "Negative Binomial regression"
selected_result = nb_result

print("\n========================================")
print("Primary model")
print("========================================")
print(selected_model_name)
print(
    "Poisson results are retained only for overdispersion "
    "diagnostics and model comparison."
)


# ============================================================
# 28. Function to convert coefficients to IRR
# ============================================================

def clean_parameter_name(
    parameter
):

    if parameter == "year5":

        return (
            "Collection year: "
            "per 5-year increase"
        )


    if (
        "ST_group"
        in parameter
    ):

        if "[T.ST2]" in parameter:

            return (
                "ST2 vs Others"
            )


    if (
        "host_group"
        in parameter
    ):

        if (
            "[T.Animal_Environment]"
            in parameter
        ):

            return (
                "Animal/Environment vs Human"
            )


    if (
        "continent"
        in parameter
    ):

        if "[T.Africa]" in parameter:

            return (
                "Africa vs North America"
            )

        if "[T.Asia]" in parameter:

            return (
                "Asia vs North America"
            )

        if "[T.Europe]" in parameter:

            return (
                "Europe vs North America"
            )

        if "[T.Oceania]" in parameter:

            return (
                "Oceania vs North America"
            )

        if (
            "[T.South America]"
            in parameter
        ):

            return (
                "South America vs North America"
            )


    return parameter


def extract_irr_table(
    result,
    model_name
):

    params = result.params.copy()
    pvalues = result.pvalues.copy()

    conf = result.conf_int()
    conf.columns = [
        "CI_low_beta",
        "CI_high_beta"
    ]


    output = pd.DataFrame(
        {
            "Parameter":
                params.index,

            "Beta":
                params.values,

            "P_value":
                pvalues.values
        }
    )


    output = output.merge(
        conf,
        left_on="Parameter",
        right_index=True,
        how="left"
    )


    # Remove intercept
    output = output[
        output["Parameter"]
        !=
        "Intercept"
    ].copy()


    # Negative Binomial includes alpha.
    # alpha is NOT a covariate IRR.
    output = output[
        output["Parameter"]
        !=
        "alpha"
    ].copy()


    output["Comparison"] = (
        output["Parameter"]
        .apply(
            clean_parameter_name
        )
    )


    output["IRR"] = np.exp(
        output["Beta"]
    )

    output["CI_low"] = np.exp(
        output["CI_low_beta"]
    )

    output["CI_high"] = np.exp(
        output["CI_high_beta"]
    )


    output.insert(
        0,
        "Model",
        model_name
    )


    output = output[
        [
            "Model",
            "Comparison",
            "IRR",
            "CI_low",
            "CI_high",
            "P_value",
            "Beta",
            "Parameter"
        ]
    ]


    return output


# ============================================================
# 29. Export Poisson IRR
# ============================================================

poisson_irr = extract_irr_table(
    poisson_result,
    "Poisson regression"
)


poisson_irr.to_csv(
    output_dir
    /
    "05_poisson_IRR.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 30. Export Negative Binomial IRR
# ============================================================

nb_irr = extract_irr_table(
    nb_result,
    "Negative Binomial regression"
)


nb_irr.to_csv(
    output_dir
    /
    "06_negative_binomial_IRR.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 31. Export selected-model IRR
# ============================================================

selected_irr = extract_irr_table(
    selected_result,
    selected_model_name
)


selected_irr.to_csv(
    output_dir
    /
    "07_PRIMARY_MODEL_IRR.csv",
    index=False,
    encoding="utf-8-sig"
)


print("\n========================================")
print("Primary model IRR")
print("========================================")

print(
    selected_irr[
        [
            "Comparison",
            "IRR",
            "CI_low",
            "CI_high",
            "P_value"
        ]
    ].to_string(
        index=False
    )
)


# ============================================================
# 32. Global likelihood-ratio tests
#
# This is especially useful for continent,
# because continent has multiple levels.
#
# Full model:
#
# NUM_FOUND ~ ST + Host + Continent + Year
#
# For each variable, refit a reduced Negative Binomial (NB2) model
# without that entire variable.
# ============================================================

terms = {
    "ST_group":
        f'C(ST_group, Treatment(reference="{ST_REFERENCE}"))',

    "host_group":
        f'C(host_group, Treatment(reference="{HOST_REFERENCE}"))',

    "continent":
        f'C(continent, Treatment(reference="{CONTINENT_REFERENCE}"))',

    "year":
        "year5"
}


def fit_selected_model(
    formula_string
):

    # The manuscript primary model is fixed as Negative Binomial (NB2).
    # Reduced models for global likelihood-ratio tests use the same
    # likelihood family as the full model.
    model = smf.negativebinomial(
        formula=formula_string,
        data=analysis_df,
        loglike_method="nb2"
    )

    result = model.fit(
        disp=False,
        maxiter=500
    )

    return result


full_result = selected_result

global_tests = []


for variable, term in terms.items():

    reduced_terms = [
        value
        for key, value in terms.items()
        if key != variable
    ]


    reduced_formula = (
        "NUM_FOUND ~ "
        +
        " + ".join(
            reduced_terms
        )
    )


    reduced_result = fit_selected_model(
        reduced_formula
    )


    lr_statistic = (
        2
        *
        (
            full_result.llf
            -
            reduced_result.llf
        )
    )


    # Difference in number of fitted parameters.
    # For Negative Binomial, alpha appears in both
    # models and therefore cancels out.
    df_difference = (
        len(full_result.params)
        -
        len(reduced_result.params)
    )


    global_p = chi2.sf(
        lr_statistic,
        df_difference
    )


    global_tests.append(
        {
            "Variable":
                variable,

            "LR_statistic":
                lr_statistic,

            "df":
                df_difference,

            "Global_P_value":
                global_p,

            "Model":
                selected_model_name
        }
    )


global_tests_df = pd.DataFrame(
    global_tests
)


global_tests_df.to_csv(
    output_dir
    /
    "08_global_likelihood_ratio_tests.csv",
    index=False,
    encoding="utf-8-sig"
)


print("\n========================================")
print("Global likelihood-ratio tests")
print("========================================")

print(
    global_tests_df.to_string(
        index=False
    )
)


# ============================================================
# 33. Expected ARG count for each genome
#
# These are fitted expected counts from the model.
# They are NOT an external-validation performance measure.
# ============================================================

analysis_df[
    "Predicted_ARG_count"
] = selected_result.predict(
    analysis_df
)


prediction_output = analysis_df[
    [
        "Genome_ID",
        "NUM_FOUND",
        "Predicted_ARG_count",
        "ST_group",
        "host_group",
        "continent",
        "year"
    ]
].copy()


prediction_output.to_csv(
    output_dir
    /
    "09_genome_level_expected_ARG_counts.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 34. Basic model-fit comparison:
# observed vs fitted values
# ============================================================

observed = (
    analysis_df["NUM_FOUND"]
    .astype(float)
)

predicted = (
    analysis_df[
        "Predicted_ARG_count"
    ]
    .astype(float)
)


mae = np.mean(
    np.abs(
        observed
        -
        predicted
    )
)


rmse = np.sqrt(
    np.mean(
        (
            observed
            -
            predicted
        ) ** 2
    )
)


mean_observed = observed.mean()
mean_predicted = predicted.mean()


fit_summary = pd.DataFrame(
    {
        "Metric": [
            "Observed_mean_ARG_count",
            "Predicted_mean_ARG_count",
            "MAE_in_sample",
            "RMSE_in_sample"
        ],
        "Value": [
            mean_observed,
            mean_predicted,
            mae,
            rmse
        ]
    }
)


fit_summary.to_csv(
    output_dir
    /
    "10_model_fit_summary.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 35. Save regression dataset
# ============================================================

analysis_df.to_csv(
    output_dir
    /
    "11_ARG_count_regression_dataset.csv",
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 36. Save text summaries
# ============================================================

with open(
    output_dir
    /
    "12_poisson_model_summary.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        poisson_result.summary()
        .as_text()
    )


with open(
    output_dir
    /
    "13_negative_binomial_model_summary.txt",
    "w",
    encoding="utf-8"
) as f:

    f.write(
        nb_result.summary()
        .as_text()
    )


# ============================================================
# 37. Final screen output
# ============================================================

print("\n\n========================================")
print("Analysis completed")
print("========================================")

print(
    f"Regression N: "
    f"{len(analysis_df):,}"
)

print(
    f"Mean NUM_FOUND: "
    f"{analysis_df['NUM_FOUND'].mean():.3f}"
)

print(
    f"Variance NUM_FOUND: "
    f"{analysis_df['NUM_FOUND'].var(ddof=1):.3f}"
)

print(
    f"Poisson dispersion: "
    f"{poisson_dispersion:.3f}"
)

print(
    f"Negative Binomial alpha: "
    f"{alpha_estimate:.6f}"
)

print(
    f"Poisson AIC: "
    f"{poisson_result.aic:.3f}"
)

print(
    f"Negative Binomial AIC: "
    f"{nb_result.aic:.3f}"
)

print(
    f"Selected primary model: "
    f"{selected_model_name}"
)


print("\nOutput directory:")

print(
    output_dir.resolve()
)


print("\nGenerated files:")

for file in sorted(
    output_dir.iterdir()
):

    print(file.name)