# Reproducibility materials for the *Acinetobacter baumannii* genomic AMR study

## 1. Scope

This package contains the processed tables and the representative custom scripts used for assembly quality control, species verification, genomic de-replication, ARG/MGE annotation and post-processing, MLST, temporal sensitivity analysis, multivariable regression, and figure generation.

Raw public genome assemblies (FASTA) are **not** redistributed. All file names, row counts, column counts, thresholds, and argument signatures recorded below were verified directly against the contents of this package.

### 1.1 Verification basis

Every quantitative statement in this README was checked against the deposited files:

| Statement type                         | Verification method                                              |
| -------------------------------------- | ---------------------------------------------------------------- |
| Row/column counts                      | Full read of each `data/*.csv` with a CSV parser (quoting-aware) |
| Accession-set identity                 | Exact set comparison of `Genome_ID` / `Assembly` columns         |
| Script arguments, defaults, thresholds | Direct inspection of the deposited script source                 |
| Figure→file mapping                    | Direct inspection of each deposited notebook's code cells        |

### 1.2 Code availability convention

The package is organised around **distinct analytical workflows**, not one file per manuscript panel. Where several figures use the same calculations and differ only in stratification or plotted outcome, one representative implementation is deposited. The same applies to repeated data-processing operations: one representative implementation is sufficient when the underlying logic is unchanged across genes, groups, or distance cutoffs.

Key harmonisation and post-processing rules used in this study are:

```text
host_1 grouping:   Human / Animal / Other
ARG matrix:        presence = 1, absence = 0
ARG-MGE pairing:   same Genome_ID + same contig + MGE within ARG +/- 5 kb
```

The ARG-MGE pairing rule is already **baked into** the deposited `ARG_MGE_*_filtered.csv` and `ARGs_MGEs_merged.csv` tables; the code that applied it is not part of this package (see Section 10).

## 2. Repository structure

Actual deposited layout, verified from the archive:

```text
A.baumannii/
├── README.md
├── data/                                  26 CSV tables, including
│   ├── metadata.csv                       final accession-level dataset (19,861 x 11)
│   └── metadata before mash.csv           pre-de-replication metadata (39,932 x 6)
└── scripts/
    ├── 01_QC/
    │   ├── qc_assemblies.sh               generic QC driver (env-var configurable)
    │   └── xaa_qc_assemblies.sh           SLURM wrapper, production invocation
    ├── 02_dereplication/
    │   └── country_hostgroup_year_dereplication.sh
    ├── 03_ARGs_annotation/
    │   └── card.sh                        ABRicate / CARD ARG annotation
    ├── 04_MGEs_annotation/
    │   ├── MGEs.sh                        BLASTn MGE annotation
    │   └── MGEs_FINAL_99perc_trim.fasta   MGE reference database (2,714 sequences)
    ├── 05_MLST/
    │   └── MLST.sh                        mlst typing, Pasteur scheme
    ├── 06_temporal_analysis/
    │   └── global_temporal_downsampling.py
    ├── 07_statistics/
    │   ├── ARG_burden_negative_binomial_public.py
    │   └── ARG_logistic_firth_regression_clean.py
    ├── 08_Figures/                        15 notebooks + 1 Python script
    └── standard_world_map-main/           world map shapefile + LICENSE + README
```

Two files are easy to confuse: `data/metadata.csv` is the **final** 19,861-genome accession table described in Section 3, while `data/metadata before mash.csv` is the **pre-de-replication** 39,932-genome table. An earlier draft of this package used the name `final_accession_metadata.csv` for the former and a 4-column, 39,932-row `metadata.csv` for a copy of the latter; only the two tables named above are deposited now.

The directory numbering now follows the analysis order: QC, de-replication, ARG annotation, MGE annotation, MLST, temporal analysis, statistics, figures. The path quoted in each section below is the deposited path.

## 3. Accession-level dataset (`data/metadata.csv`)

The final accession-level dataset is deposited as `data/metadata.csv`. Verified dimensions:

```text
Rows:                        19,861
Columns:                     11
Unique Assembly accessions:  19,861
Unique #BioSample values:    19,861
File size:                   1,903,685 bytes
Encoding:                    UTF-8 without BOM
Delimiter:                   comma
```

Columns, with the value counts verified from the deposited file:

| Column        | Distinct | Empty | Description                                                             |
| ------------- | --------:| -----:| ----------------------------------------------------------------------- |
| `#BioSample`  | 19,861   | 0     | NCBI BioSample accession; unique, one per genome                        |
| `Assembly`    | 19,861   | 0     | NCBI Assembly accession; the study-wide genome identifier               |
| `host`        | 25       | 0     | collection host, kept at its original granularity (see Section 3.3)     |
| `country`     | 106      | 0     | country of collection                                                   |
| `continent`   | 6        | 0     | continent of collection                                                 |
| `SNP cluster` | 1,732    | 3,536 | NCBI Pathogen Detection SNP cluster ID, e.g. `PDS000111232.1`           |
| `Length`      | 19,452   | 0     | assembly length in bp; range 3,601,717-4,299,988                        |
| `BioProject`  | 1,209    | 0     | NCBI BioProject accession                                               |
| `Contigs`     | 300      | 0     | number of contigs; range 1-300                                          |
| `year`        | 54       | 0     | collection year; range 1900-2026                                        |
| `ST`          | 598      | 0     | MLST sequence type, Pasteur scheme; `-` marks an unknown type (n = 980) |

Distributions of the two categorical columns used for stratification:

```text
continent:  North America 8,523 | Asia 6,936 | Europe 2,969 | South America 564
            Africa 485 | Oceania 384

host (top): human 17,719 | others 1,561 | bird 251 | cattle 109 | cat 26
            horse 26 | duck 25 | Environment 22 | wolf 22 | pig 19 | ...
```

Note the units of the two QC-related columns: `Length` is the assembly length in bp (so the 3.6-4.3 Mb QC window in Section 5 appears as 3,600,000-4,300,000 here), and `Contigs` is a count.

### 3.1 Verified relations to the other tables

The deposited file was validated against every other table by exact accession join:

```text
Assembly set identical to the Genome_ID set of acquired_resistance_genes.csv:      True (19,861)
Assembly set identical to the Genome_ID set of merged_ARG_metadata.csv:            True (19,861)
Assembly set identical to the Genome_ID set of acquired_ARGs_sum_with_metadata.csv: True (19,861)
Assembly set is a subset of metadata before mash.csv (39,932 accessions):          True
Assembly set is a subset of antibiotic_classes_matched_metadata.csv (19,861):      True
Genomes present here but missing from ARGs_MGEs_merged.csv:                        18
Genomes present here but missing from MGE_count_AB.csv:                            18
```

Field-by-field agreement with the source tables, over all 19,861 accessions:

```text
country    vs merged_ARG_metadata.csv / metadata before mash.csv:   0 mismatches
continent  vs merged_ARG_metadata.csv:                             0 mismatches
year       vs merged_ARG_metadata.csv / metadata before mash.csv:   0 mismatches
ST         vs antibiotic_classes_matched_metadata.csv:             0 mismatches
ST         vs acquired_ARGs_sum_with_metadata.csv:                 0 mismatches
#BioSample vs metadata before mash.csv:                            0 mismatches
```

The two columns with non-zero differences are covered in Section 3.2 (`ST`) and Section 3.3 (`host`).

### 3.2 Important note on the `ST` column

The column named `ST` in `merged_ARG_metadata.csv` is **not** the MLST sequence type. It is a pre-computed three-level grouping variable with exactly three values:

```text
ST = "2"        n = 11,706   -> ST2
ST = "Others"   n =  7,175
ST = "Unknown"  n =    980
```

It is accompanied by an explicit `ST_group` column with the values `ST2` / `Others` / `Unknown`.

The `ST` column of `metadata.csv`, `antibiotic_classes_matched_metadata.csv`, and `acquired_ARGs_sum_with_metadata.csv` is the **real** MLST sequence type: 598 named STs plus the `-` code for an unknown type, i.e. 599 distinct values, with `-` occurring 980 times. That count of unknowns matches the `Unknown` count in the `merged_ARG_metadata.csv` grouping variable exactly, and the two representations agree for all 19,861 accessions. When reporting ST distributions, use the `ST` column of `metadata.csv`; when reporting the ST2-versus-others regression term, use `ST_group`.

### 3.3 Important note on the `host` column

`metadata.csv` carries a single column named `host`, and it is the **most granular** host representation in the package. It is not the harmonised three-category variable used for de-replication stratification (Section 6). Its values and counts are:

```text
human            17,719
others            1,561
bird                251
cattle              109
cat                  26
horse                26
duck                 25
Environment          22
wolf                 22
pig                  19
earthworm            18
chicken              14
dog                  12
Lizard               11
sheep                 6
snake                 4
Soricidae             3
mice                  3
goose                 2
rabbit                2
fish                  2
turtle                1
ferret                1
parrot                1
Equus caballus        1
```

Relations to the other tables, verified accession by accession:

- Against `metadata before mash.csv` and `antibiotic_classes_matched_metadata.csv` (`host_1` in both), there are exactly 22 differences, all of one kind: a record that is `others` in those tables is `Environment` here.
- Against `antibiotic_classes_matched_metadata.csv` (`Host`), there are 59 differences, again all of the same kind: records collapsed to `others` there are resolved here to `Environment` (22), `Lizard` (11), `sheep` (6), `snake` (4), `Soricidae` (3), `mice` (3), `goose` (2), `rabbit` (2), `fish` (2), `turtle`, `ferret`, `parrot`, `Equus caballus` (1 each).
- Against `merged_ARG_metadata.csv` (`host_1`), there are 559 differences, because that table collapses all non-human records to `Animals`; here the animal terms are resolved individually.
- The `host_group` values actually used in the regression are `Human` / `Animal_Environment` / `Unknown`, and their counts are exactly the harmonisation 17,719 / 581 / 1,561 reported in Section 6. The 581 therefore decomposes as 251 `bird` + 109 `cattle` + 26 `cat` + 26 `horse` + 25 `duck` + 22 `wolf` + 19 `pig` + 18 `earthworm` + 14 `chicken` + 12 `dog` + 11 `Lizard` + 6 `sheep` + 4 `snake` + 3 `Soricidae` + 3 `mice` + 2 `goose` + 2 `rabbit` + 2 `fish` + `turtle` + `ferret` + `parrot` + `Equus caballus` + the 22 `Environment` records.

The harmonisation rule that reproduces the counts used downstream is therefore:

```text
human                              -> Human
others, missing, unspecified       -> Other-or-unspecified
Environment and every animal term  -> Animal_Environment
```

Because the manuscript text uses "Human / Animal / Other", the mapping from these 25 raw values to the three reported categories must be stated explicitly rather than left implicit.

### 3.4 Relation to the richer source metadata

`Fig.1_continent_and country_distribution.ipynb` stores a cached output from a **20-column** version of this metadata table:

```text
#BioSample, Assembly, host, host_1, Location, country, continent, Isolation source,
SNP cluster, Length, Virulence genotypes, AMR genotypes, AST phenotypes,
Isolation type, PFGE primary enzyme pattern, BioProject, Contigs,
Collection date, year, ST
```

The deposited `metadata.csv` is an 11-column selection from that table. Nine columns are not deposited: `host_1`, `Location`, `Isolation source`, `Virulence genotypes`, `AMR genotypes`, `AST phenotypes`, `Isolation type`, `PFGE primary enzyme pattern`, and `Collection date`. The `host_1` and `Location` fields are recoverable from `metadata before mash.csv` (which holds `host_1` and `Location` for all 39,932 pre-de-replication accessions, with 0 `#BioSample` mismatches against the final set). The isolation-source, genotype, phenotype, PFGE, and collection-date fields exist nowhere in this package, so no analysis that depends on them is reproducible from these materials.

## 4. Deposited data tables

All 26 files in `data/` were profiled in full. Row counts exclude the header row.

| File                                      | Rows   | Columns | Main role                                   |
| ----------------------------------------- | ------:| -------:| ------------------------------------------- |
| `metadata.csv`                            | 19,861 | 11      | final accession-level dataset (Section 3)   |
| `metadata before mash.csv`                | 39,932 | 6       | pre-de-replication metadata (with Location) |
| `acquired_resistance_genes.csv`           | 19,861 | 244     | ARG count and binary ARG matrix             |
| `merged_ARG_metadata.csv`                 | 19,861 | 227     | regression metadata + ARG matrix            |
| `acquired_ARGs_sum_with_metadata.csv`     | 19,861 | 23      | ARG burden/time and figure support          |
| `ARGs_matched_representatives.csv`        | 19,861 | 423     | figure-level ARG presence/absence           |
| `antibiotic_classes_matched_metadata.csv` | 19,861 | 25      | antibiotic-class figure metadata + true ST  |
| `carbapenem_matched_representatives.csv`  | 19,861 | 24      | Fig. 2 carbapenem profiles                  |
| `ARGs_MGEs_merged.csv`                    | 19,843 | 56      | ARG-MGE Spearman analysis                   |
| `MGE_count_AB.csv`                        | 19,843 | 90      | MGE count/network support                   |
| `ARG_MGE_1_filtered.csv`                  | 38,977 | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_2_filtered.csv`                  | 47,821 | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_3_filtered.csv`                  | 46,906 | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_4_filtered.csv`                  | 41,844 | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_5_filtered.csv`                  | 49,859 | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_6_filtered.csv`                  | 37,763 | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_7_filtered.csv`                  | 41,429 | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_8_filtered.csv`                  | 4,111  | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_9_filtered.csv`                  | 39,913 | 11      | ARG-MGE pairs, distance-filtered            |
| `ARG_MGE_10_filtered.csv`                 | 31,130 | 11      | ARG-MGE pairs, distance-filtered            |
| `Global_Year.csv`                         | 4      | 17      | ARG class burden by time period             |
| `ST_year_Final.csv`                       | 20     | 28      | ST distribution over years (wide format)    |
| `gene_bar_chart.csv`                      | 37     | 3       | per-gene totals and percentages             |
| `resistance_gene_counts.csv`              | 181    | 3       | gene family counts by type                  |
| `sankey.csv`                              | 37     | 5       | ARG -> drug class -> mechanism flows        |
| `acquired_ARGs_classify.csv`              | 33     | 3       | ARG class / mechanism / gene membership     |

Notes on individual tables:

- `metadata.csv` is the final accession-level dataset and the entry point for the whole package: its `Assembly` set is identical to the `Genome_ID` set of `acquired_resistance_genes.csv`, `merged_ARG_metadata.csv`, and `acquired_ARGs_sum_with_metadata.csv`. Its columns and the differences from the older metadata tables are documented in Section 3, and `Length` / `Contigs` in that file satisfy the QC thresholds of Section 5 (3,601,717-4,299,988 bp, 1-300 contigs).
- The ten `ARG_MGE_*_filtered.csv` files share one schema: `Genome_ID, SEQUENCE, ARG_Gene, STRAND, ARG_Start, ARG_End, mge_id, MGE_ID, MGE_Start, MGE_End, Distance_bp`. The maximum absolute `Distance_bp` in these files is 5,000 bp, confirming the +/- 5 kb pairing window; `SEQUENCE` is the contig identifier that enforces the same-contig requirement. The `_1` ... `_10` suffixes index ARG/MGE annotation partitions and are recombined by `Fig.4_RoseChart_and_pairs_proportion.ipynb`.
- `ARGs_MGEs_merged.csv` and `MGE_count_AB.csv` cover 19,843 genomes, i.e. 18 fewer than the 19,861-genome final set. Those 18 genomes have no ARG-MGE pair record, which is expected but should be stated in the manuscript when reporting ARG-MGE network and Spearman sample sizes.
- `acquired_ARGs_classify.csv` is **GB18030-encoded**, not UTF-8, and contains Chinese `Mechanism` labels. It must be read with `encoding='gb18030'` (or `'gbk'`).
- `metadata before mash.csv` has one duplicated `#BioSample` value (39,932 rows, 39,916 unique `#BioSample`; all 39,932 `Assembly` accessions are unique). The de-replication script resolves this deterministically.
- `ST_year_Final.csv` contains raw Chinese text in its earliest-period column header (`≤2000`) and should be read with an explicit encoding.

One non-CSV data file is also deposited, inside the annotation directory:

| File                                                      | Content                                                        | Role                            |
| --------------------------------------------------------- | -------------------------------------------------------------- | ------------------------------- |
| `scripts/04_MGEs_annotation/MGEs_FINAL_99perc_trim.fasta` | 2,714 MGE reference sequences (Pärnänen database, 99% trimmed) | BLASTn query for MGE annotation |

## 5. Assembly QC and species verification

Scripts:

```text
scripts/01_QC/qc_assemblies.sh          generic driver, configured through environment variables
scripts/01_QC/xaa_qc_assemblies.sh      SLURM wrapper, records the production invocation
```

The QC script is a self-contained Bash driver. It enumerates assemblies (`.fna`/`.fa`/`.fasta`, optionally `.gz`) from `ASM_DIR`, verifies gzip integrity, computes assembly statistics with `seqkit stats -a`, applies the length/contig/N thresholds, and optionally adds a `fastANI` screen against a reference genome. The final PASS list is the intersection of the gzip, assembly-statistics, and (when enabled) ANI filters.

Settings and their deposited defaults:

| Variable      | Default            | Role                                          |
| ------------- | ------------------ | --------------------------------------------- |
| `ASM_DIR`     | `wsh`              | directory holding the assemblies              |
| `OUTDIR`      | `O`                | QC output directory                           |
| `THREADS`     | `32`               | GNU `parallel` concurrency                    |
| `MIN_LEN`     | `3600000` (3.6 Mb) | minimum assembly size                         |
| `MAX_LEN`     | `4300000` (4.3 Mb) | maximum assembly size                         |
| `MAX_CONTIGS` | `300`              | maximum number of contigs (`num_seqs`)        |
| `MAX_N_FRAC`  | `0.05`             | maximum fraction of `N` (`sum_gap / sum_len`) |
| `DO_ANI`      | `1`                | enable the `fastANI` screen                   |
| `REF_FOR_ANI` | `ref.fa`           | ANI reference genome                          |
| `ANI_MIN`     | `95.0`             | ANI acceptance threshold (%)                  |

Thresholds as applied:

| Criterion                                          | Threshold  |
| -------------------------------------------------- | ----------:|
| Assembly size                                      | 3.6-4.3 Mb |
| Number of contigs                                  | <300       |
| Ambiguous bases (N)                                | <5%        |
| ANI to *A. baumannii* ATCC 19606 (GCF_009035845.1) | >=95%      |

The production wrapper `xaa_qc_assemblies.sh` contains the exact study invocation:

```text
ASM_DIR=/public/home/wangjm01/zzf/Ab_project/xaa_dir
OUTDIR=qc_xaa
DO_ANI=1
REF_FOR_ANI=/public/home/wangjm01/zzf/Reference_genome/GCF_009035845.1_ASM903584v1_genomic.fna
ANI_MIN=95.0
```

confirming that the ANI reference is *A. baumannii* ATCC 19606 (GCF_009035845.1 / ASM903584v1). The wrapper carries SLURM directives (`--nodes=2`, `--ntasks-per-node=16`, `--mem-per-cpu=6gb`, `--partition=vip`) and site-specific absolute paths, which must be edited before reuse elsewhere.

Software used by the QC step:

```text
SeqKit   2.12.0   (the script requires seqkit >= 2.5.1 column layout)
fastANI  1.34
```

Requirements and output files:

- `awk`, GNU `parallel`, `seqkit` are required; `fastANI` is required only when `DO_ANI=1`.
- Outputs written to `OUTDIR`: `asm_all.txt`, `gz_list.txt`, `bad_gzip.txt`, `asm_good_gzip.txt`, `seqkit_stats.tsv`, `seqkit_stats.flag.tsv` (per-assembly PASS/FAIL plus the failing reason), `pass_by_seqkit.txt`, `fail_by_seqkit.txt`, `ani/` and `ani_summary.tsv` (when `DO_ANI=1`), `pass_by_ani.txt`, `pass_final.txt`, `pass_list.txt`, and `summary.txt`.

Three details in the script need to be read carefully when reproducing or reusing it:

1. Two comments are stale leftovers from a different organism: `MIN_LEN` is annotated "金葡 2.6 Mb" and `MAX_LEN` "金葡 3.2 Mb" (i.e. *S. aureus* examples). The actual values, 3,600,000 and 4,300,000 bp, are the *A. baumannii* thresholds listed above.
2. The `MAX_N_FRAC` comment says "N 比例上限（1%）" but the value is `0.05`, i.e. 5%. The code uses the value, so the effective threshold is 5%.
3. A comment near `ANI_MIN` says "按阈值 98.0%（可改）", but the default and the production wrapper both use 95.0%.

The script header still calls the final list the input "供后续 SNP 流水线使用" (for a downstream SNP pipeline); in this study it is the QC PASS list that feeds de-replication and annotation.

Cross-check against the deposited dataset: the `Length` and `Contigs` columns of `metadata.csv` reproduce the QC outcome. All 19,861 retained genomes have `Length` between 3,601,717 and 4,299,988 bp (inside the 3.6-4.3 Mb window) and `Contigs` between 1 and 300 (none exceed the 300-contig cap). The `N`-fraction and ANI criteria cannot be re-checked from the deposited tables, because neither quantity is recorded in any deposited file.

## 6. Genomic de-replication

Script:

```text
scripts/02_dereplication/country_hostgroup_year_dereplication.sh
```

Invocation:

```text
bash country_hostgroup_year_dereplication.sh <metadata.csv> <mash_candidate_db.sqlite> <output_dir> [dnadiff_cache.sqlite]
```

Required metadata columns:

```text
Assembly
#BioSample
country
host_1
year
```

Thresholds and rule, as implemented:

```text
Mash distance <= 0.001                     (candidate-pair identification only)
dnadiff 1-to-1 AvgIdentity >= 99.99%
reference aligned coverage >= 99%
query aligned coverage >= 99%
Clustering: representative-anchored (no single-linkage propagation) --
            every removed genome must directly satisfy the identity and
            coverage thresholds against its retained representative
```

Versions:

```text
Mash     2.3
dnadiff  1.3
nucmer   4.0.1
```

The script is a SLURM batch script (`#SBATCH` directives are retained, with `--cpus-per-task=32`, `--mem=64G`) and runs its own embedded Python via `python3`; `python3` and `dnadiff` must be on `PATH`.

The `host_1` field is a **pre-harmonised** host-group variable grouped into:

```text
Human / Animal / Other-or-unspecified
```

De-replication is stratified by:

```text
country + host_1 + year
```

Note that the deposited `data/metadata.csv` is the **final** 19,861-genome table and has no `host_1` column. The de-replication script's `<metadata.csv>` argument expects the *pre*-de-replication input table, of which `data/metadata before mash.csv` (39,932 rows with `host_1` and `Location`) is the deposited example — the placeholder name in the script's usage line is unrelated to the deposited file of that name.

The harmonisation itself was performed **before** this script ran and is not represented by any deposited script. The values observed in `data/metadata before mash.csv` are still partially granular (`human`, `others`, plus animal terms such as `bird`, `cattle`, `cat`, `duck`, `wolf`, `horse`, `pig`, `earthworm`, `chicken`, `dog`, `Lizard`, `sheep`, `snake`, `Soricidae`, `mice`, `goose`, `rabbit`, `fish`, `turtle`, `ferret`, `parrot`, `Equus caballus`), so anyone reproducing the stratification must apply the harmonisation explicitly. Applying the grouping `human` -> Human, `others`/missing -> Other, all remaining animal terms -> Animal reproduces the three-group counts that are actually used downstream:

```text
Human  17,719
Animal    581
Other   1,561   (unknown/unspecified host; retained in the dataset, excluded from host-adjusted regression)
```

### 6.1 Intermediate input required

The script requires a precomputed Mash candidate SQLite database (`mash_candidates_le_001.sqlite` in the usage example). **This database is not deposited**, and neither are the commands that generated it. Before permanent deposition, provide either the database or the exact script/commands that produced it. The optional fourth argument is a previously generated dnadiff results database, which can be reused read-only; reusing it changes computation time only, not the thresholds.

## 7. ARG, MGE, and MLST annotation

The annotation scripts for all three modules are deposited in this package. Also deposited is the MGE reference database used by the MGE annotation step. What is **not** deposited are the input assemblies, the CARD database snapshot, and the MLST scheme files, so the steps are documented and re-runnable but the exact upstream corpora must be supplied from their own sources.

### 7.1 ARG annotation (ABRicate / CARD)

Script:

```text
scripts/03_ARGs_annotation/card.sh
```

Command as implemented:

```text
abricate --db card --mincov 90 --minid 80 ./*.fna > card.result.tab
abricate --summary card.result.tab > card.summary.tab
```

Settings and versions:

```text
ABRicate                 1.0.1
CARD database            2,631 nucleotide sequences
ABRicate database date   2026-01-30
Minimum coverage         90%   (--mincov 90)
Minimum identity         80%   (--minid 80)
```

This 90% coverage / 80% identity pair is the ARG detection threshold applied upstream of the binary presence/absence matrices, and should be quoted together with the CARD version in the manuscript. The script runs on `./*.fna` in the working directory and carries SLURM directives (`--ntasks=32`, `--nodes=2`, `--partition=vip`). It does not set an explicit ABRicate database version, so the CARD version above must be pinned separately for exact reproduction.

### 7.2 MGE annotation (BLASTn)

Script:

```text
scripts/04_MGEs_annotation/MGEs.sh
Reference database: scripts/04_MGEs_annotation/MGEs_FINAL_99perc_trim.fasta
```

Command as implemented, executed once per genome:

```text
makeblastdb -in <sample>.fna -dbtype nucl -out ./<sample>/<sample>_db
blastn -query MGEs_FINAL_99perc_trim.fasta -db ./<sample>/<sample>_db \
       -outfmt 6 -perc_identity 90 -qcov_hsp_perc 85 \
       -out ./<sample>/<sample>_mge.txt
```

Settings and versions:

```text
BLAST                    2.16.0
Output format            -outfmt 6 (tabular)
Minimum identity         90%   (-perc_identity 90)
Minimum query coverage   85%   (-qcov_hsp_perc 85)
Query                    MGEs_FINAL_99perc_trim.fasta
```

Deposited MGE reference database:

```text
File            scripts/04_MGEs_annotation/MGEs_FINAL_99perc_trim.fasta
Sequences       2,714
Header format   ><index>_<element>_<GenBank accession>   e.g. >414_tnpA_AL513383.1
Element types   tnpA / tnpA1 / IS / transposase and plasmid-replicon families
Source          Pärnänen MobileGeneticElementDatabase, clustered/trimmed at 99% identity
```

The script expects each sample as an individual FASTA in `xaa/` and writes one output directory per sample. It carries SLURM directives (`--nodes=1`, `--ntasks-per-node=16`, `--mem-per-cpu=4gb`, `--partition=vip`, `--bb=1000`).

### 7.3 MLST annotation

Script:

```text
scripts/05_MLST/MLST.sh
```

Command as implemented:

```text
mlst ./*.fna --scheme abaumannii --csv > mlst_Ab.csv
```

Settings and versions:

```text
mlst                     2.23.0
Scheme                   A. baumannii Pasteur (--scheme abaumannii)
Output                   CSV
GrapeTree                1.5.0   (tree visualisation, run outside this script)
```

The typed sequence types are the `ST` column of `antibiotic_classes_matched_metadata.csv` and `acquired_ARGs_sum_with_metadata.csv` (598 named STs plus the `-` unknown code; see Section 3.1). The `ST` column of `merged_ARG_metadata.csv` is a different, pre-computed grouping variable and must not be confused with this output.

## 8. Multivariable analyses

### 8.1 ARG carriage

Script:

```text
scripts/07_statistics/ARG_logistic_firth_regression_clean.py
```

Input: `merged_ARG_metadata.csv` (resolved relative to the working directory).

Model:

```text
ARG carriage ~ ST_group + host_group + continent + year
```

Analysed set, as implemented:

```text
ST_group    in {ST2, Others}                 (Unknown excluded)
host_group  in {Human, Animal_Environment}   (Unknown excluded)
continent   in {North America, Asia, Europe, Africa, Oceania, South America}
year        non-missing
```

Note the coding: the host factor in the regression dataset is `Animal_Environment`, and `Unknown` host records (1,561 genomes) are **not** dropped from the master dataset but are excluded from the regression. `Unknown` ST records (980 genomes) are excluded in the same way.

Reference groups:

```text
ST_group:   Others
host_group: Human
continent:  North America
```

The year effect is reported per 5-year increase (`Year_per_5_years = (year - median(year)) / 5`; the centring affects only numerical stability and the intercept, not the reported effect). Twelve target ARGs are analysed:

```text
Standard multivariable logistic regression (8):
  blaOXA-23-like, tet(B), armA, sul2, sul1, TEM-12, APH(3')-VIa, AAC(3)-Ia

Firth penalised logistic regression (4, sparse positives):
  blaOXA-24-like, blaOXA-58-like, NDM-1, AAC(6')-Ib7
```

Benjamini-Hochberg FDR correction (`statsmodels`, `method="fdr_bh"`, `alpha=0.05`) is applied to the global tests for these targets. The script also runs a Firth global penalised likelihood-ratio test for the Firth genes. ARG columns are validated to contain only 0/1 before fitting.

### 8.2 Acquired ARG burden

Script:

```text
scripts/07_statistics/ARG_burden_negative_binomial_public.py
```

Inputs: `acquired_resistance_genes.csv` (for `NUM_FOUND`) merged with `merged_ARG_metadata.csv`.

Primary model:

```text
Negative Binomial NB2
NUM_FOUND ~ ST_group + host_group + continent + year
```

Settings, as implemented:

```text
Reference groups:  ST_group = Others, host_group = Human, continent = North America
Year scaling:      YEAR_CENTER = 2015, YEAR_SCALE = 5  -> IRR per 5-year increase
Overdispersion:    Poisson Pearson chi-square diagnostic, threshold 1.5
                   (descriptive only; does not change the primary model)
```

The Poisson model is retained only as an overdispersion diagnostic; the NB2 model is primary. `NUM_FOUND` is taken from the original `acquired_resistance_genes.csv` rather than from the merged metadata, missing `NUM_FOUND` values are not converted to zero, and records with invalid or non-integer counts are excluded before fitting. The estimated NB2 `alpha` is reported. Input reading falls back through `utf-8-sig`, `utf-8`, `gb18030`, and `gbk` encodings.

## 9. Temporal equal-n down-sampling

Script:

```text
scripts/06_temporal_analysis/global_temporal_downsampling.py
```

Required input columns:

```text
year
NUM_FOUND
```

Settings, as implemented (the period boundaries are hard-coded constants and should not be changed when reproducing the analysis):

```text
Periods:        <=2000 | 2001-2009 | 2010-2019 | >=2020
N_RESAMPLES =   1000
RANDOM_SEED =   20260916
Sampling =      without replacement, equal n per period per iteration
Equal n =       size of the smallest temporal period
```

The script writes `period_sample_sizes.tsv`, the per-iteration down-sampled means, the summary statistics (including a Kruskal-Wallis test across periods), and a run-provenance text file recording the seed, the number of resamples, and the per-period sample sizes.

`INPUT_FILE` is shipped as the placeholder `"/path/to/your/final_dataset.csv"` and the script raises an explicit error until it is set. Point it at a table carrying `year` and `NUM_FOUND`; `data/acquired_ARGs_sum_with_metadata.csv` is the natural choice because it contains both columns for all 19,861 genomes.

This is a sensitivity analysis for unequal temporal sample size; it should not be described as correcting geographic, surveillance, lineage, or host-composition bias.

## 10. Representative figure and visualization scripts

The manuscript repeats the same plotting or stratification logic across continents, countries, sequence types, hosts, or time periods. To avoid depositing many near-identical notebooks, one representative script is deposited per visualization pattern.

The 16 deposited figure files are:

| Analysis/figure type               | Deposited file                                   | Notes                                                                            |
| ---------------------------------- | ------------------------------------------------ | -------------------------------------------------------------------------------- |
| Geographic distribution            | `Fig.1_global_distribution.ipynb`                | World-map workflow; needs the bundled shapefile (Section 11)                     |
| Country/continent distribution     | `Fig.1_continent_and country_distribution.ipynb` | Grouped-distribution workflow; note the space in the filename                    |
| Temporal distribution              | `Fig.1_time_distribution.py`                     | Log-scale genome-count-per-year scatter                                          |
| Host distribution                  | `Fig1._host_distribution.ipynb`                  | Host composition; groups into Human / Animal / Others                            |
| Carbapenem-gene profiles           | `Fig.2_carbapenem.ipynb`                         | Five ST/country/host/year-stratified profile variants in one notebook            |
| ARG carriage prevalence            | `Fig.3_ARG_carriage_prevalence.ipynb`            | Prevalence by drug class over time periods                                       |
| ARG-MGE pair composition           | `Fig.4_RoseChart_and_pairs_proportion.ipynb`     | Recombines the ten `ARG_MGE_*_filtered.csv` partitions; rose chart + pair shares |
| ARG-MGE network analysis           | `Fig.4_network_analysis.ipynb`                   | Spearman correlation heatmap + network construction/graph metrics                |
| ST distribution                    | `Fig.S1_ST_distribution.ipynb`                   | ST composition over years from `ST_year_Final.csv`                               |
| Resistance-gene pie charts         | `Fig.S1_pie_chart.ipynb`                         | Gene-family composition pie charts from `resistance_gene_counts.csv`             |
| Gene bar chart                     | `Fig.S2_bar_chart.ipynb`                         | Per-gene totals from `gene_bar_chart.csv`                                        |
| ARG distribution by metadata group | `Fig.S3_ARGs_distribution.ipynb`                 | Top-N gene presence heatmap across countries                                     |
| Sankey visualization               | `Fig.S4_sankey.ipynb`                            | ARG -> drug class -> mechanism Sankey from `sankey.csv`                          |
| ST-specific ARG bubble plot        | `Fig.S5_bubble_plot.ipynb`                       | Three bubble-plot variants by ST                                                 |
| Temporal ARG burden                | `Fig.S6_acquired_ARGs_avg_year.ipynb`            | Mean acquired ARGs per year with Spearman correlation                            |
| ARG-MGE Spearman correlation       | `mge_ARG_spearman.ipynb`                         | Correlation heatmap from `ARGs_MGEs_merged.csv` (ARGs vs. MGEs)                  |

For manuscript panels that differ only by the stratifying variable (continent vs. country vs. host vs. ST) or by the displayed subset, the corresponding representative script documents the common analytical and plotting logic. Near-duplicate copies are intentionally not included. Figure-specific manual layout, panel assembly, labelling, or cosmetic editing performed outside the analytical workflow is not treated as a separate analysis script.

### 10.1 Absolute paths and unresolved dependencies

Every notebook and the temporal script open their inputs through **hard-coded absolute paths**, e.g.:

```text
D:\A.baumannii\data\metadata.csv
D:\python\敏感性分析\...            (in Fig.2_carbapenem.ipynb, Fig1._host_distribution.ipynb)
```

Two notebooks import a helper module that is **not** deposited:

```text
from src.datapreprocess import generate_summary, process_mge, decode_latlon
```

used by `Fig.4_network_analysis.ipynb` and `Fig.S6_acquired_ARGs_avg_year.ipynb`. `Fig.4_network_analysis.ipynb` additionally expects Graphviz/`pygraphviz` (`nx.nx_agraph.graphviz_layout`), and `Fig.S6_acquired_ARGs_avg_year.ipynb` imports `semopy`. None of these are declared in a dependency file. Before re-running these notebooks, replace the absolute paths with relative ones and supply `src/datapreprocess.py`.

`Fig.S6_acquired_ARGs_avg_year.ipynb` was written against an **earlier schema** of the burden table. Its stored outputs list the columns `Assembly` and `Acquired_ARGs`, and the code groups by `Acquired_ARGs`, whereas the deposited `acquired_ARGs_sum_with_metadata.csv` uses `Genome_ID` and `NUM_FOUND`. The column names must be reconciled before this notebook will run against the deposited table.

Earlier drafts of this README also listed representative files that do not exist in this package: `Fig.1_country_continent_distribution.ipynb`, `Fig.1_host_distribution.ipynb`, `Fig.4_MGE_pairs_and_composition.ipynb`, `Fig.S2_ARG_prevalence_bar.ipynb`, `Fig.S2_subtype_composition.ipynb`, and `Fig.S13_ARG_MGE_spearman.ipynb`. They have been replaced above by the deposited filenames that actually carry that logic.

## 11. Bundled world map

`scripts/standard_world_map-main/` contains a standard world map in shapefile/GeoJSON form (审图号 GS(2016)1666), distributed under the MIT License (Copyright (c) 2026 Wu Lei). It is the map source used by `Fig.1_global_distribution.ipynb`. That notebook reads the shapefile from `D:\A.baumannii\scripts\worldmap\standard_world_map-main\std_worldmap.shp`; the deposited location is `scripts/standard_world_map-main/std_worldmap.shp`, so the path must be adjusted.

## 12. Software and package versions

### Bioinformatics software

```text
SeqKit          2.12.0
fastANI         1.34
Mash            2.3
dnadiff         1.3
nucmer          4.0.1
ABRicate        1.0.1
BLAST           2.16.0
mlst            2.23.0
GrapeTree       1.5.0
```

### Reference databases

```text
CARD            2,631 nucleotide sequences; ABRicate DB date 2026-01-30
                (query database for scripts/03_ARGs_annotation/card.sh)
MGE database    Pärnänen MobileGeneticElementDatabase
                deposited as scripts/04_MGEs_annotation/MGEs_FINAL_99perc_trim.fasta
                2,714 sequences, clustered/trimmed at 99% identity
MLST scheme     A. baumannii Pasteur scheme (--scheme abaumannii)
                scheme files not deposited
```

### Shell pipeline tools

The QC, annotation, and MLST scripts are shell scripts that call external binaries. Only `scripts/01_QC/qc_assemblies.sh` states a minimum version in its source:

```text
seqkit          >= 2.5.1   (explicit requirement in qc_assemblies.sh; version used: 2.12.0)
GNU parallel    required by qc_assemblies.sh (no version pinned)
gzip            required by qc_assemblies.sh for .gz integrity checks (no version pinned)
makeblastdb     required by scripts/04_MGEs_annotation/MGEs.sh (BLAST 2.16.0)
```

### Python environment for the core statistical scripts

Versions recorded for the environment in which the core statistics were run:

```text
Python                3.13.9
pandas                3.0.6
numpy                 2.5.3
scipy                 1.18.1
statsmodels           0.15.0
patsy                 1.0.2
firthmodels           0.8.2
seaborn               0.13.2
networkx              3.6.1
geopandas             1.1.3
country_converter     1.3.2
plotly                6.6.0
tqdm                  4.67.3
```




