import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

# Firth logistic regression
from firthmodels import FirthLogisticRegression
from firthmodels.adapters.statsmodels import FirthLogit


# ============================================================
# 1. 输入文件和输出目录
# ============================================================

input_file = "merged_ARG_metadata.csv"

output_dir = Path("FINAL_regression_results")
output_dir.mkdir(exist_ok=True)


# ============================================================
# 2. 12 个目标 ARG
# ============================================================

# ------------------------------------------------------------
# 8 个普通 multivariable logistic regression
# ------------------------------------------------------------

standard_genes = [
    "blaOXA-23-like",
    "tet(B)",
    "armA",
    "sul2",
    "sul1",
    "TEM-12",
    "APH(3')-VIa",
    "AAC(3)-Ia",
]


# ------------------------------------------------------------
# 4 个 Firth logistic regression
#
# 原因：
# 某些 ST / Host / Continent 组阳性病例过少或为 0
# ------------------------------------------------------------

firth_genes = [
    "blaOXA-24-like",
    "blaOXA-58-like",
    "NDM-1",
    "AAC(6')-Ib7",
]


all_genes = standard_genes + firth_genes


# ============================================================
# 3. Reference groups
# ============================================================

# ST:
# Others = reference
# 所以得到 ST2 vs Others

# Host:
# Human = reference
# 所以得到 Animal/Environment vs Human

# Continent:
# North America = reference

# Year:
# 连续变量，每增加 5 年


# ============================================================
# 4. 读取数据
# ============================================================

df = pd.read_csv(
    input_file,
    encoding="utf-8-sig",
    low_memory=False
)


print("=" * 70)
print("Original dataset")
print("=" * 70)

print(f"Total genomes : {len(df):,}")
print(f"Total columns : {df.shape[1]:,}")


# ============================================================
# 5. 检查必要列
# ============================================================

required_metadata = [
    "Genome_ID",
    "ST_group",
    "host_group",
    "continent",
    "year",
]


missing_columns = [
    col
    for col in required_metadata + all_genes
    if col not in df.columns
]


if missing_columns:

    raise ValueError(
        "以下列在文件中不存在：\n"
        + "\n".join(missing_columns)
    )


# ============================================================
# 6. Year 转换为 numeric
# ============================================================

df["year"] = pd.to_numeric(
    df["year"],
    errors="coerce"
)


# ============================================================
# 7. 建立正式 regression dataset
#
# 只把以下样本用于模型：
#
# ST:
#   ST2
#   Others
#
# Host:
#   Human
#   Animal_Environment
#
# Continent:
#   六个大洲
#
# Year:
#   有有效年份
#
# Unknown 并没有从 master dataset 删除，
# 只是没有进入本次多因素回归。
# ============================================================

valid_continents = [
    "North America",
    "Asia",
    "Europe",
    "Africa",
    "Oceania",
    "South America",
]


analysis_df = df[
    df["ST_group"].isin(
        ["ST2", "Others"]
    )
    &
    df["host_group"].isin(
        ["Human", "Animal_Environment"]
    )
    &
    df["continent"].isin(
        valid_continents
    )
    &
    df["year"].notna()
].copy()


print("\n" + "=" * 70)
print("Regression dataset")
print("=" * 70)

print(
    f"Included genomes : "
    f"{len(analysis_df):,}"
)

print(
    f"Not included     : "
    f"{len(df) - len(analysis_df):,}"
)


# ============================================================
# 8. 检查 ARG 是否确实为 0/1
# ============================================================

for gene in all_genes:

    analysis_df[gene] = pd.to_numeric(
        analysis_df[gene],
        errors="coerce"
    )

    invalid = analysis_df.loc[
        analysis_df[gene].notna()
        &
        ~analysis_df[gene].isin([0, 1]),
        gene
    ].unique()

    if len(invalid) > 0:

        raise ValueError(
            f"{gene} 存在非 0/1 值：{invalid}"
        )


# ============================================================
# 9. Year 转换成每 5 年
#
# 中心化只是为了数值稳定和截距解释。
# 不改变 year 的 OR。
# ============================================================

YEAR_CENTER = float(
    analysis_df["year"].median()
)


analysis_df["Year_per_5_years"] = (
    analysis_df["year"] - YEAR_CENTER
) / 5.0


print(
    f"\nYear centering value = "
    f"{YEAR_CENTER:g}"
)

print(
    "Year effect = adjusted OR "
    "per 5-year increase."
)


# ============================================================
# 10. 构建设计矩阵
#
# 不使用自动 dummy coding，
# 明确写出每一个比较方向，
# 避免 reference 搞错。
# ============================================================

def build_design_matrix(data):

    X = pd.DataFrame(
        index=data.index
    )

    # --------------------------------------------------------
    # ST
    # Others = reference
    # --------------------------------------------------------

    X["ST2_vs_Others"] = (
        data["ST_group"] == "ST2"
    ).astype(float)


    # --------------------------------------------------------
    # Host
    # Human = reference
    # --------------------------------------------------------

    X["AnimalEnv_vs_Human"] = (
        data["host_group"]
        ==
        "Animal_Environment"
    ).astype(float)


    # --------------------------------------------------------
    # Continent
    # North America = reference
    # --------------------------------------------------------

    X["Asia_vs_NorthAmerica"] = (
        data["continent"] == "Asia"
    ).astype(float)


    X["Europe_vs_NorthAmerica"] = (
        data["continent"] == "Europe"
    ).astype(float)


    X["Africa_vs_NorthAmerica"] = (
        data["continent"] == "Africa"
    ).astype(float)


    X["Oceania_vs_NorthAmerica"] = (
        data["continent"] == "Oceania"
    ).astype(float)


    X["SouthAmerica_vs_NorthAmerica"] = (
        data["continent"]
        ==
        "South America"
    ).astype(float)


    # --------------------------------------------------------
    # Year
    # --------------------------------------------------------

    X["Year_per_5_years"] = (
        data["Year_per_5_years"]
        .astype(float)
    )


    return X.astype(float)


# ============================================================
# 11. Predictor groups
#
# 用于 overall/global test。
#
# Continent 一共有 5 个 dummy，
# 所以 global test 是 5 df。
# ============================================================

predictor_groups = {

    "ST": [
        "ST2_vs_Others"
    ],

    "Host": [
        "AnimalEnv_vs_Human"
    ],

    "Continent": [
        "Asia_vs_NorthAmerica",
        "Europe_vs_NorthAmerica",
        "Africa_vs_NorthAmerica",
        "Oceania_vs_NorthAmerica",
        "SouthAmerica_vs_NorthAmerica",
    ],

    "Year": [
        "Year_per_5_years"
    ],
}


# ============================================================
# 12. 输出标签
# ============================================================

effect_labels = {

    "ST2_vs_Others": {
        "Predictor": "ST",
        "Comparison": "ST2 vs Others",
        "Reference": "Others",
    },

    "AnimalEnv_vs_Human": {
        "Predictor": "Host",
        "Comparison":
            "Animal/Environment vs Human",
        "Reference": "Human",
    },

    "Asia_vs_NorthAmerica": {
        "Predictor": "Continent",
        "Comparison":
            "Asia vs North America",
        "Reference": "North America",
    },

    "Europe_vs_NorthAmerica": {
        "Predictor": "Continent",
        "Comparison":
            "Europe vs North America",
        "Reference": "North America",
    },

    "Africa_vs_NorthAmerica": {
        "Predictor": "Continent",
        "Comparison":
            "Africa vs North America",
        "Reference": "North America",
    },

    "Oceania_vs_NorthAmerica": {
        "Predictor": "Continent",
        "Comparison":
            "Oceania vs North America",
        "Reference": "North America",
    },

    "SouthAmerica_vs_NorthAmerica": {
        "Predictor": "Continent",
        "Comparison":
            "South America vs North America",
        "Reference": "North America",
    },

    "Year_per_5_years": {
        "Predictor": "Year",
        "Comparison":
            "Per 5-year increase",
        "Reference": "Continuous",
    },
}


# ============================================================
# 13. 格式化 P value
# ============================================================

def format_p_value(p):

    if pd.isna(p):
        return ""

    if p < 0.001:
        return "<0.001"

    return f"{p:.3f}"


# ============================================================
# 14. 普通 Logistic 的 global likelihood-ratio test
# ============================================================

def standard_global_lr_test(
    y,
    X,
    full_result,
    columns_to_drop
):

    reduced_X = X.drop(
        columns=columns_to_drop
    )

    reduced_X = sm.add_constant(
        reduced_X,
        has_constant="add"
    )


    reduced_model = sm.Logit(
        y,
        reduced_X
    )


    reduced_result = reduced_model.fit(
        disp=False,
        maxiter=500
    )


    LR = 2 * (
        full_result.llf
        -
        reduced_result.llf
    )


    LR = max(
        float(LR),
        0.0
    )


    df_test = len(
        columns_to_drop
    )


    p_value = chi2.sf(
        LR,
        df_test
    )


    return (
        LR,
        df_test,
        p_value
    )


# ============================================================
# 15. Firth estimator
#
# 用于 global penalized likelihood-ratio tests。
# ============================================================

def fit_firth_estimator(
    X,
    y
):

    model = FirthLogisticRegression(

        fit_intercept=True,

        backend="auto",

        max_iter=100,

        # 比默认容差更严格一点
        gtol=1e-7,

        xtol=1e-7,

        max_halfstep=50,

    )


    model.fit(
        X,
        y
    )


    if not model.converged_:

        raise RuntimeError(
            "Firth model did not converge."
        )


    return model


# ============================================================
# 16. Firth global penalized likelihood-ratio test
#
# 完整模型与 reduced model 的
# penalized log-likelihood 比较。
# ============================================================

def firth_global_plr_test(
    y,
    X,
    full_model,
    columns_to_drop
):

    reduced_X = X.drop(
        columns=columns_to_drop
    )


    reduced_model = fit_firth_estimator(
        reduced_X,
        y
    )


    LR = 2 * (
        full_model.loglik_
        -
        reduced_model.loglik_
    )


    LR = max(
        float(LR),
        0.0
    )


    df_test = len(
        columns_to_drop
    )


    p_value = chi2.sf(
        LR,
        df_test
    )


    return (
        LR,
        df_test,
        p_value
    )


# ============================================================
# 17. 结果容器
# ============================================================

effect_results = []

global_results = []

diagnostic_results = []


# ============================================================
# 18. 循环跑 12 个 ARG
# ============================================================

for gene in all_genes:

    print("\n" + "=" * 70)
    print(f"ARG: {gene}")
    print("=" * 70)


    # --------------------------------------------------------
    # 当前基因的数据
    # --------------------------------------------------------

    gene_df = analysis_df[
        [
            gene,
            "ST_group",
            "host_group",
            "continent",
            "Year_per_5_years",
        ]
    ].dropna().copy()


    y = gene_df[
        gene
    ].astype(float)


    X = build_design_matrix(
        gene_df
    )


    n_total = len(
        gene_df
    )

    n_positive = int(
        (y == 1).sum()
    )

    n_negative = int(
        (y == 0).sum()
    )


    prevalence = (
        n_positive
        /
        n_total
        *
        100
    )


    print(
        f"N = {n_total:,}"
    )

    print(
        f"Positive = {n_positive:,}"
    )

    print(
        f"Negative = {n_negative:,}"
    )

    print(
        f"Prevalence = "
        f"{prevalence:.2f}%"
    )


    # ========================================================
    # A. 普通 Logistic Regression
    # ========================================================

    if gene in standard_genes:

        print(
            "Model: Standard multivariable logistic regression"
        )


        X_const = sm.add_constant(
            X,
            has_constant="add"
        )


        with warnings.catch_warnings():

            warnings.simplefilter(
                "ignore"
            )

            full_model = sm.Logit(
                y,
                X_const
            )

            full_result = full_model.fit(
                disp=False,
                maxiter=500
            )


        converged = (
            full_result
            .mle_retvals
            .get(
                "converged",
                True
            )
        )


        if not converged:

            raise RuntimeError(
                f"{gene}: standard logistic "
                f"did not converge."
            )


        # ----------------------------------------------------
        # coefficient-level:
        #
        # Beta
        # aOR
        # Wald 95% CI
        # Wald P
        # ----------------------------------------------------

        ci = full_result.conf_int(
            alpha=0.05
        )


        for term in X.columns:

            beta = float(
                full_result.params[
                    term
                ]
            )


            lower_beta = float(
                ci.loc[
                    term,
                    0
                ]
            )


            upper_beta = float(
                ci.loc[
                    term,
                    1
                ]
            )


            p_value = float(
                full_result.pvalues[
                    term
                ]
            )


            info = effect_labels[
                term
            ]


            OR = np.exp(
                beta
            )

            OR_lower = np.exp(
                lower_beta
            )

            OR_upper = np.exp(
                upper_beta
            )


            effect_results.append({

                "ARG":
                    gene,

                "Model":
                    "Standard logistic regression",

                "Predictor":
                    info["Predictor"],

                "Comparison":
                    info["Comparison"],

                "Reference":
                    info["Reference"],

                "Beta":
                    beta,

                "Adjusted_OR":
                    OR,

                "CI95_lower":
                    OR_lower,

                "CI95_upper":
                    OR_upper,

                "P_value":
                    p_value,

                "P_method":
                    "Wald",

                "CI_method":
                    "Wald 95% CI",
            })


        # ----------------------------------------------------
        # global likelihood-ratio tests
        # ----------------------------------------------------

        for predictor, columns in (
            predictor_groups.items()
        ):

            LR, df_test, p_global = (
                standard_global_lr_test(
                    y,
                    X,
                    full_result,
                    columns
                )
            )


            global_results.append({

                "ARG":
                    gene,

                "Model":
                    "Standard logistic regression",

                "Predictor":
                    predictor,

                "LR_statistic":
                    LR,

                "df":
                    df_test,

                "Global_P":
                    p_global,

                "Global_test_method":
                    "Likelihood-ratio test",
            })


        diagnostic_results.append({

            "ARG":
                gene,

            "Model":
                "Standard logistic regression",

            "N":
                n_total,

            "Positive":
                n_positive,

            "Negative":
                n_negative,

            "Prevalence_percent":
                prevalence,

            "Converged":
                converged,

        })


    # ========================================================
    # B. Firth Logistic Regression
    # ========================================================

    else:

        print(
            "Model: Firth penalized logistic regression"
        )


        # ----------------------------------------------------
        # Statsmodels-style adapter
        #
        # pl=True:
        #
        # P:
        # Penalized likelihood-ratio test
        #
        # CI:
        # Profile penalized-likelihood CI
        # ----------------------------------------------------

        X_const = sm.add_constant(
            X,
            has_constant="add"
        )


        firth_result = FirthLogit(
            y,
            X_const
        ).fit(
            pl=True
        )


        # ----------------------------------------------------
        # 把各种返回类型统一成 pandas
        # ----------------------------------------------------

        if isinstance(
            firth_result.params,
            pd.Series
        ):

            params = (
                firth_result.params.copy()
            )

        else:

            params = pd.Series(
                np.asarray(
                    firth_result.params
                ),
                index=X_const.columns
            )


        if isinstance(
            firth_result.pvalues,
            pd.Series
        ):

            pvalues = (
                firth_result.pvalues.copy()
            )

        else:

            pvalues = pd.Series(
                np.asarray(
                    firth_result.pvalues
                ),
                index=X_const.columns
            )


        ci_raw = (
            firth_result.conf_int()
        )


        if isinstance(
            ci_raw,
            pd.DataFrame
        ):

            ci_firth = (
                ci_raw.copy()
            )

            # 保证顺序和设计矩阵一致
            ci_firth.index = (
                X_const.columns
            )

            ci_firth.columns = [
                "lower",
                "upper"
            ]

        else:

            ci_firth = pd.DataFrame(

                np.asarray(
                    ci_raw
                ),

                index=X_const.columns,

                columns=[
                    "lower",
                    "upper"
                ]
            )


        # ----------------------------------------------------
        # 提取 coefficient-level Firth results
        # ----------------------------------------------------

        for term in X.columns:

            beta = float(
                params[
                    term
                ]
            )


            lower_beta = float(
                ci_firth.loc[
                    term,
                    "lower"
                ]
            )


            upper_beta = float(
                ci_firth.loc[
                    term,
                    "upper"
                ]
            )


            p_value = float(
                pvalues[
                    term
                ]
            )


            info = effect_labels[
                term
            ]


            OR = np.exp(
                beta
            )

            OR_lower = np.exp(
                lower_beta
            )

            OR_upper = np.exp(
                upper_beta
            )


            effect_results.append({

                "ARG":
                    gene,

                "Model":
                    "Firth logistic regression",

                "Predictor":
                    info["Predictor"],

                "Comparison":
                    info["Comparison"],

                "Reference":
                    info["Reference"],

                "Beta":
                    beta,

                "Adjusted_OR":
                    OR,

                "CI95_lower":
                    OR_lower,

                "CI95_upper":
                    OR_upper,

                "P_value":
                    p_value,

                "P_method":
                    "Penalized likelihood-ratio",

                "CI_method":
                    "Profile penalized-likelihood 95% CI",
            })


        # ----------------------------------------------------
        # 再使用 sklearn-style estimator
        #
        # 获取 penalized log-likelihood，
        # 做 ST / Host / Continent / Year global PLR。
        # ----------------------------------------------------

        full_firth_model = (
            fit_firth_estimator(
                X,
                y
            )
        )


        for predictor, columns in (
            predictor_groups.items()
        ):

            LR, df_test, p_global = (
                firth_global_plr_test(
                    y,
                    X,
                    full_firth_model,
                    columns
                )
            )


            global_results.append({

                "ARG":
                    gene,

                "Model":
                    "Firth logistic regression",

                "Predictor":
                    predictor,

                "LR_statistic":
                    LR,

                "df":
                    df_test,

                "Global_P":
                    p_global,

                "Global_test_method":
                    "Penalized likelihood-ratio test",
            })


        diagnostic_results.append({

            "ARG":
                gene,

            "Model":
                "Firth logistic regression",

            "N":
                n_total,

            "Positive":
                n_positive,

            "Negative":
                n_negative,

            "Prevalence_percent":
                prevalence,

            "Converged":
                bool(
                    full_firth_model.converged_
                ),

        })


# ============================================================
# 19. 生成结果表
# ============================================================

effects_df = pd.DataFrame(
    effect_results
)


global_df = pd.DataFrame(
    global_results
)


diagnostics_df = pd.DataFrame(
    diagnostic_results
)


# ============================================================
# 20. 对 48 个主要 association tests 做 BH-FDR
#
# 12 ARG × 4 predictors
#
# ST
# Host
# Continent
# Year
# ============================================================

if len(global_df) != 48:

    print(
        "\nWARNING:"
        f" Expected 48 global tests, "
        f"but found {len(global_df)}."
    )


reject, q_values, _, _ = (
    multipletests(
        global_df[
            "Global_P"
        ].values,
        alpha=0.05,
        method="fdr_bh"
    )
)


global_df[
    "FDR_q_value"
] = q_values


global_df[
    "FDR_significant_0.05"
] = reject


# ============================================================
# 21. P/q 显示格式
#
# 原始 numeric P 仍然保留，
# 这两列只是方便看论文结果。
# ============================================================

effects_df[
    "P_display"
] = effects_df[
    "P_value"
].apply(
    format_p_value
)


global_df[
    "Global_P_display"
] = global_df[
    "Global_P"
].apply(
    format_p_value
)


global_df[
    "FDR_q_display"
] = global_df[
    "FDR_q_value"
].apply(
    format_p_value
)


# ============================================================
# 22. 检查 P 与 CI 是否一致
#
# Firth:
# PLR P + profile PL CI
#
# 理论上在 0.05 水平应基本一致。
# ============================================================

effects_df[
    "CI_excludes_1"
] = (
    (
        effects_df[
            "CI95_lower"
        ] > 1
    )
    |
    (
        effects_df[
            "CI95_upper"
        ] < 1
    )
)


effects_df[
    "P_significant_0.05"
] = (
    effects_df[
        "P_value"
    ] < 0.05
)


effects_df[
    "P_CI_consistent"
] = (
    effects_df[
        "CI_excludes_1"
    ]
    ==
    effects_df[
        "P_significant_0.05"
    ]
)


# ============================================================
# 23. 把 global P 和 FDR 加入 effect-level 表
#
# 这样最终一张表同时可以看到：
#
# aOR
# 95% CI
# coefficient P
# predictor global P
# global FDR q
# ============================================================

global_for_merge = global_df[
    [
        "ARG",
        "Predictor",
        "Global_P",
        "FDR_q_value",
        "FDR_significant_0.05",
    ]
].copy()


final_effects = effects_df.merge(

    global_for_merge,

    on=[
        "ARG",
        "Predictor"
    ],

    how="left"
)


# ============================================================
# 24. 排序
# ============================================================

gene_order = {
    gene: i
    for i, gene
    in enumerate(
        all_genes
    )
}


predictor_order = {
    "ST": 1,
    "Host": 2,
    "Continent": 3,
    "Year": 4,
}


comparison_order = {

    "ST2 vs Others": 1,

    "Animal/Environment vs Human": 1,

    "Asia vs North America": 1,
    "Europe vs North America": 2,
    "Africa vs North America": 3,
    "Oceania vs North America": 4,
    "South America vs North America": 5,

    "Per 5-year increase": 1,
}


final_effects[
    "_gene_order"
] = (
    final_effects["ARG"]
    .map(gene_order)
)


final_effects[
    "_predictor_order"
] = (
    final_effects[
        "Predictor"
    ].map(
        predictor_order
    )
)


final_effects[
    "_comparison_order"
] = (
    final_effects[
        "Comparison"
    ].map(
        comparison_order
    )
)


final_effects = (
    final_effects
    .sort_values(
        [
            "_gene_order",
            "_predictor_order",
            "_comparison_order",
        ]
    )
    .drop(
        columns=[
            "_gene_order",
            "_predictor_order",
            "_comparison_order",
        ]
    )
)


global_df[
    "_gene_order"
] = (
    global_df["ARG"]
    .map(
        gene_order
    )
)


global_df[
    "_predictor_order"
] = (
    global_df[
        "Predictor"
    ].map(
        predictor_order
    )
)


global_df = (
    global_df
    .sort_values(
        [
            "_gene_order",
            "_predictor_order",
        ]
    )
    .drop(
        columns=[
            "_gene_order",
            "_predictor_order",
        ]
    )
)


diagnostics_df[
    "_gene_order"
] = (
    diagnostics_df[
        "ARG"
    ].map(
        gene_order
    )
)


diagnostics_df = (
    diagnostics_df
    .sort_values(
        "_gene_order"
    )
    .drop(
        columns="_gene_order"
    )
)


# ============================================================
# 25. 四舍五入
#
# P value 不做 round，
# 防止非常小的 P 被变成 0。
# ============================================================

round_columns = [
    "Beta",
    "Adjusted_OR",
    "CI95_lower",
    "CI95_upper",
]


for col in round_columns:

    final_effects[col] = (
        final_effects[col]
        .astype(float)
        .round(4)
    )


global_df[
    "LR_statistic"
] = (
    global_df[
        "LR_statistic"
    ]
    .astype(float)
    .round(4)
)


diagnostics_df[
    "Prevalence_percent"
] = (
    diagnostics_df[
        "Prevalence_percent"
    ]
    .round(4)
)


# ============================================================
# 26. 保存最终结果
# ============================================================

effects_file = (
    output_dir
    /
    "01_FINAL_12_ARG_adjusted_OR.csv"
)


global_file = (
    output_dir
    /
    "02_FINAL_12_ARG_global_tests_FDR.csv"
)


diagnostics_file = (
    output_dir
    /
    "03_FINAL_model_diagnostics.csv"
)


final_effects.to_csv(
    effects_file,
    index=False,
    encoding="utf-8-sig"
)


global_df.to_csv(
    global_file,
    index=False,
    encoding="utf-8-sig"
)


diagnostics_df.to_csv(
    diagnostics_file,
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 27. 单独输出 4 个 Firth 基因
#
# 方便检查旧结果和新结果的变化
# ============================================================

firth_only = final_effects[
    final_effects[
        "ARG"
    ].isin(
        firth_genes
    )
].copy()


firth_file = (
    output_dir
    /
    "04_FIRTH_ONLY_profile_CI_PLR_P.csv"
)


firth_only.to_csv(
    firth_file,
    index=False,
    encoding="utf-8-sig"
)


# ============================================================
# 28. 屏幕输出
# ============================================================

print("\n\n" + "=" * 70)
print("MODEL DIAGNOSTICS")
print("=" * 70)

print(
    diagnostics_df.to_string(
        index=False
    )
)


print("\n\n" + "=" * 70)
print("FOUR FIRTH MODELS")
print("=" * 70)

print(
    firth_only[
        [
            "ARG",
            "Predictor",
            "Comparison",
            "Adjusted_OR",
            "CI95_lower",
            "CI95_upper",
            "P_display",
            "P_CI_consistent",
        ]
    ].to_string(
        index=False
    )
)


print("\n\n" + "=" * 70)
print("GLOBAL TESTS + BH-FDR")
print("=" * 70)

print(
    global_df[
        [
            "ARG",
            "Predictor",
            "Global_P_display",
            "FDR_q_display",
            "FDR_significant_0.05",
        ]
    ].to_string(
        index=False
    )
)


print("\n\n" + "=" * 70)
print("ANALYSIS COMPLETED")
print("=" * 70)

print(
    f"Results directory:\n"
    f"{output_dir.resolve()}"
)

print(
    "\nGenerated files:"
)

print(
    f"1. {effects_file.name}"
)

print(
    f"2. {global_file.name}"
)

print(
    f"3. {diagnostics_file.name}"
)

print(
    f"4. {firth_file.name}"
)