#!/bin/bash
#SBATCH --job-name=AB_CHY_derep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=AB_CHY_derep_%j.log
#SBATCH --partition=vip
#SBATCH --bb=1000

set -euo pipefail

# ============================================================
# A. baumannii country-host-group-year stratified de-replication
#
# Purpose:
#   1) Read genome accessions from metadata column "Assembly".
#   2) If #BioSample is duplicated, keep one eligible genome
#      deterministically (first occurrence in the input metadata).
#   3) Require country and collection year for stratified de-replication.
#      The input column "host_1" is a PRE-HARMONIZED host-group variable:
#        Human / Animal / Other-or-unspecified.
#      In the study dataset, environmental records and records with
#      missing/unspecified original host metadata were encoded in the
#      Other-or-unspecified group before this script was run; therefore
#      missing original host metadata is not itself an exclusion criterion.
#   4) Restrict genomic de-replication to identical
#      country + harmonized host group + collection year strata.
#   5) Use Mash distance <= 0.001 only to identify candidate pairs.
#   6) Confirm redundancy with dnadiff using:
#        1-to-1 AvgIdentity >= 99.99%
#        REF coverage >= 99%
#        QRY coverage >= 99%
#   7) Use representative-anchored clustering (no single-linkage
#      propagation): every removed genome must directly satisfy the
#      identity and coverage thresholds against its retained representative.
#
# IMPORTANT:
#   - "host_1" must already be harmonized before running this script.
#   - Candidate pairs may be supplied through a precomputed Mash SQLite DB.
#   - Existing dnadiff results can optionally be reused read-only.
#   - Reusing caches changes computation time only, not the thresholds.
# ============================================================

# Usage:
#   bash 01_country_hostgroup_year_dereplication.sh \
#       metadata.csv mash_candidates_le_001.sqlite output_dir [dnadiff_cache.sqlite]
#
# Required metadata columns:
#   Assembly, #BioSample, country, host_1, year
#
# NOTE:
#   host_1 must be the harmonized host-group variable used in the study.
#   Environmental and missing/unspecified original-host records should
#   already have been assigned to the Other/unspecified group.

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <metadata.csv> <mash_candidate_db.sqlite> <output_dir> [dnadiff_cache.sqlite]"
    exit 2
fi

META="$1"
CANDIDATE_DB="$2"
OUT_WORK="$3"
OLD_DNADIFF_DB="${4:-}"

mkdir -p "${OUT_WORK}"

SUPP_DB="${OUT_WORK}/dnadiff_supplemental_country_hostgroup_year.sqlite"

SELECTED_TSV="${OUT_WORK}/selected_complete_unique_biosample.tsv"
MISSING_TSV="${OUT_WORK}/excluded_missing_country_hostgroup_year.tsv"
DUP_BIOSAMPLE_TSV="${OUT_WORK}/duplicate_biosample_removed.tsv"
DUP_ASSEMBLY_TSV="${OUT_WORK}/duplicate_assembly_removed.tsv"
UNMATCHED_TSV="${OUT_WORK}/metadata_assembly_not_found.tsv"

CLUSTER_TSV="${OUT_WORK}/clusters_country_hostgroup_year.tsv"
REP_TXT="${OUT_WORK}/representatives_country_hostgroup_year.txt"
REP_TSV="${OUT_WORK}/representatives_country_hostgroup_year.tsv"
REMOVED_TXT="${OUT_WORK}/removed_country_hostgroup_year.txt"
SUMMARY_TXT="${OUT_WORK}/summary_country_hostgroup_year.txt"
ERROR_LOG="${OUT_WORK}/dnadiff_errors_country_hostgroup_year.log"

THREADS="${SLURM_CPUS_PER_TASK:-32}"
DNADIFF_WORKERS="${THREADS}"
DNADIFF_BATCH=512

if [ -n "${SLURM_TMPDIR:-}" ]; then
    DNADIFF_TMP="${SLURM_TMPDIR}/AB_CHY_dnadiff_${SLURM_JOB_ID:-manual}"
else
    DNADIFF_TMP="/tmp/${USER}_AB_CHY_dnadiff_${SLURM_JOB_ID:-manual}"
fi
mkdir -p "${DNADIFF_TMP}"

cleanup() {
    rm -rf "${DNADIFF_TMP}" || true
}
trap cleanup EXIT

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "============================================================"
echo "A. baumannii country-host-group-year stratified de-replication"
echo "============================================================"
echo "Start time:          $(date)"
echo "Host:                $(hostname)"
echo "Metadata:            ${META}"
echo "Candidate DB:        ${CANDIDATE_DB}"
echo "Old dnadiff cache:   ${OLD_DNADIFF_DB:-none}"
echo "Output workdir:      ${OUT_WORK}"
echo "CPUs/workers:        ${DNADIFF_WORKERS}"
echo "Final NI cutoff:     99.99%"
echo "Coverage cutoff:     99% BOTH"
echo "Stratification:      country + harmonized host group (host_1) + year"
echo

for F in "${META}" "${CANDIDATE_DB}"; do
    if [ ! -s "${F}" ]; then
        echo "ERROR: required file not found or empty:"
        echo "${F}"
        exit 1
    fi
done

if [ -n "${OLD_DNADIFF_DB}" ] && [ ! -s "${OLD_DNADIFF_DB}" ]; then
    echo "ERROR: optional dnadiff cache was specified but not found or empty:"
    echo "${OLD_DNADIFF_DB}"
    exit 1
fi

for CMD in python3 dnadiff; do
    if ! command -v "${CMD}" >/dev/null 2>&1; then
        echo "ERROR: ${CMD} not found"
        exit 1
    fi
done

python3 - \
    "${META}" \
    "${CANDIDATE_DB}" \
    "${OLD_DNADIFF_DB}" \
    "${SUPP_DB}" \
    "${SELECTED_TSV}" \
    "${MISSING_TSV}" \
    "${DUP_BIOSAMPLE_TSV}" \
    "${DUP_ASSEMBLY_TSV}" \
    "${UNMATCHED_TSV}" \
    "${CLUSTER_TSV}" \
    "${REP_TXT}" \
    "${REP_TSV}" \
    "${REMOVED_TXT}" \
    "${SUMMARY_TXT}" \
    "${ERROR_LOG}" \
    "${DNADIFF_TMP}" \
    "${DNADIFF_WORKERS}" \
    "${DNADIFF_BATCH}" \
<<'PY'
import csv
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

(
    metadata_csv,
    candidate_database,
    old_cache_database,
    supplemental_database,
    selected_tsv,
    missing_tsv,
    dup_biosample_tsv,
    dup_assembly_tsv,
    unmatched_tsv,
    cluster_tsv,
    representative_txt,
    representative_tsv,
    removed_txt,
    summary_txt,
    error_log,
    temp_root,
    max_workers,
    dnadiff_batch,
) = sys.argv[1:]

MAX_WORKERS = int(max_workers)
DNADIFF_BATCH = int(dnadiff_batch)

MASH_CUTOFF = 0.001
NI_CUTOFF = 99.99
COV_CUTOFF = 99.0

START = time.time()

MISSING_TOKENS = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "missing",
    "unknown",
    "not available",
    "not applicable",
    "-",
    ".",
}

ACCESSION_RE = re.compile(r"(GC[AF]_[0-9]+\.[0-9]+)", re.I)
YEAR_RE = re.compile(r"(?<!\d)((?:18|19|20)\d{2})(?!\d)")


def clean(v):
    if v is None:
        return ""
    return " ".join(str(v).strip().split())


def is_missing(v):
    return clean(v).casefold() in MISSING_TOKENS


def norm_text(v):
    return clean(v).casefold()


def norm_year(v):
    s = clean(v)
    if is_missing(s):
        return None
    m = YEAR_RE.search(s)
    if not m:
        return None
    y = int(m.group(1))
    # Fixed validation bounds keep historical reruns deterministic.
    if y < 1800 or y > 2100:
        return None
    return str(y)


def extract_accession(text):
    m = ACCESSION_RE.search(os.path.basename(clean(text)))
    if not m:
        return None
    return m.group(1).upper()


def canonical_pair(a, b):
    return (a, b) if a <= b else (b, a)


def extract_percent(text):
    m = re.search(r"\(([0-9.]+)%\)", text)
    if not m:
        raise ValueError("Cannot parse percentage from: " + text)
    return float(m.group(1))


def parse_report(report_file):
    cov_ref = None
    cov_query = None
    ni = None
    in_one_to_one = False

    with open(report_file) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            fields = line.split()

            if fields[0] == "AlignedBases":
                if len(fields) < 3:
                    raise ValueError("Malformed AlignedBases line")
                cov_ref = extract_percent(fields[1])
                cov_query = extract_percent(fields[2])

            elif fields[0] == "1-to-1":
                in_one_to_one = True

            elif in_one_to_one and fields[0] == "AvgIdentity":
                ni = float(fields[1])
                in_one_to_one = False

    if ni is None:
        raise ValueError("1-to-1 AvgIdentity not found")
    if cov_ref is None or cov_query is None:
        raise ValueError("AlignedBases coverage not found")

    return ni, cov_ref, cov_query


def passes_threshold(ni, cov1, cov2):
    return (
        ni is not None
        and cov1 is not None
        and cov2 is not None
        and ni >= NI_CUTOFF
        and cov1 >= COV_CUTOFF
        and cov2 >= COV_CUTOFF
    )


def run_dnadiff(genome1, genome2, mash_distance):
    g1, g2 = canonical_pair(genome1, genome2)
    workdir = tempfile.mkdtemp(prefix="dnadiff_", dir=temp_root)
    try:
        prefix = os.path.join(workdir, "out")
        cmd = ["dnadiff", "-p", prefix, g1, g2]

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        env["OPENBLAS_NUM_THREADS"] = "1"
        env["MKL_NUM_THREADS"] = "1"
        env["NUMEXPR_NUM_THREADS"] = "1"

        result = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )

        report = prefix + ".report"

        if result.returncode != 0:
            return {
                "g1": g1, "g2": g2, "mash": mash_distance,
                "ni": None, "cov1": None, "cov2": None,
                "status": "ERROR",
                "error": f"dnadiff exit {result.returncode}: {result.stderr[-2000:]}",
            }

        if not os.path.isfile(report):
            return {
                "g1": g1, "g2": g2, "mash": mash_distance,
                "ni": None, "cov1": None, "cov2": None,
                "status": "ERROR",
                "error": "report not generated",
            }

        try:
            ni, cov1, cov2 = parse_report(report)
        except Exception as exc:
            return {
                "g1": g1, "g2": g2, "mash": mash_distance,
                "ni": None, "cov1": None, "cov2": None,
                "status": "ERROR",
                "error": "Report parsing failed: " + str(exc),
            }

        return {
            "g1": g1, "g2": g2, "mash": mash_distance,
            "ni": ni, "cov1": cov1, "cov2": cov2,
            "status": "OK", "error": "",
        }

    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ============================================================
# 1. Open candidate DB and build Assembly -> genome mapping
# ============================================================

candidate_conn = sqlite3.connect(
    "file:" + candidate_database + "?mode=ro",
    uri=True,
    timeout=60,
)
candidate_conn.execute("PRAGMA busy_timeout=60000")
candidate_conn.execute("PRAGMA cache_size=-1000000")
candidate_conn.execute("PRAGMA mmap_size=4294967296")

genome_rows = candidate_conn.execute(
    """
    SELECT id, path, degree
    FROM genomes
    ORDER BY id
    """
).fetchall()

id_to_path = {}
id_to_degree = {}
accession_to_entries = {}

for genome_id, path, degree in genome_rows:
    id_to_path[genome_id] = path
    id_to_degree[genome_id] = degree

    acc = extract_accession(path)
    if acc is None:
        continue
    accession_to_entries.setdefault(acc, []).append((genome_id, path, degree))

print(f"Candidate DB genomes: {len(genome_rows):,}", flush=True)


# ============================================================
# 2. Read metadata
# ============================================================

with open(metadata_csv, "r", encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle)

    if reader.fieldnames is None:
        raise RuntimeError("Metadata CSV has no header")

    fields = [f.strip() for f in reader.fieldnames]
    required = ["Assembly", "#BioSample", "country", "host_1", "year"]

    missing_cols = [c for c in required if c not in fields]
    if missing_cols:
        raise RuntimeError(
            "Missing required metadata columns: " + ", ".join(missing_cols)
        )

    # DictReader keys may contain surrounding whitespace; normalize them.
    rows = []
    for row_no, raw in enumerate(reader, start=2):
        row = {}
        for k, v in raw.items():
            if k is not None:
                row[k.strip()] = clean(v)

        row["_row_no"] = row_no
        row["_assembly"] = extract_accession(row.get("Assembly", ""))
        row["_biosample"] = clean(row.get("#BioSample", ""))
        row["_country"] = clean(row.get("country", ""))
        row["_host"] = clean(row.get("host_1", ""))
        row["_year"] = norm_year(row.get("year", ""))

        row["_country_key"] = norm_text(row["_country"])
        row["_host_key"] = norm_text(row["_host"])

        rows.append(row)

print(f"Metadata rows: {len(rows):,}", flush=True)


# ============================================================
# 3. Match metadata Assembly to genome files from candidate DB
# ============================================================

matched = []
unmatched = []

for row in rows:
    acc = row["_assembly"]
    entries = accession_to_entries.get(acc or "", [])

    if not acc or not entries:
        unmatched.append(row)
        continue

    # Normally one path per Assembly. If multiple are present, choose the
    # lexicographically first path deterministically and report via stdout.
    entries = sorted(entries, key=lambda x: x[1])
    genome_id, genome_path, degree = entries[0]

    row["_genome_id"] = genome_id
    row["_genome_path"] = genome_path
    row["_degree"] = degree
    matched.append(row)

with open(unmatched_tsv, "w", newline="", encoding="utf-8") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow(["metadata_row", "Assembly", "#BioSample", "country", "host_1", "year"])
    for r in unmatched:
        w.writerow([
            r["_row_no"], r.get("Assembly", ""), r["_biosample"],
            r["_country"], r["_host"], r.get("year", "")
        ])


# ============================================================
# 4. Validate required stratification metadata
#
# In the study dataset, host_1 had already been harmonized, so no genome
# was excluded because of missing original host metadata. This defensive
# check remains to prevent malformed future inputs.
# ============================================================

complete = []
missing_meta = []

for row in matched:
    reasons = []

    if is_missing(row["_country"]):
        reasons.append("country")
    if is_missing(row["_host"]):
        reasons.append("host_1")
    if row["_year"] is None:
        reasons.append("year")

    if reasons:
        row["_missing_reason"] = ",".join(reasons)
        missing_meta.append(row)
    else:
        row["_stratum"] = (
            row["_country_key"],
            row["_host_key"],
            row["_year"],
        )
        complete.append(row)

with open(missing_tsv, "w", newline="", encoding="utf-8") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow([
        "metadata_row", "Assembly", "#BioSample",
        "country", "host_1", "year", "missing_fields", "genome_path"
    ])
    for r in missing_meta:
        w.writerow([
            r["_row_no"], r["_assembly"], r["_biosample"],
            r["_country"], r["_host"], r.get("year", ""),
            r["_missing_reason"], r["_genome_path"]
        ])


# ============================================================
# 5. Remove duplicate Assembly rows, keeping first CSV occurrence
# ============================================================

seen_assembly = set()
assembly_unique = []
dup_assembly = []

for row in complete:
    acc = row["_assembly"]
    if acc in seen_assembly:
        dup_assembly.append(row)
        continue
    seen_assembly.add(acc)
    assembly_unique.append(row)

with open(dup_assembly_tsv, "w", newline="", encoding="utf-8") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow([
        "metadata_row", "Assembly", "#BioSample",
        "country", "host_1", "year", "genome_path"
    ])
    for r in dup_assembly:
        w.writerow([
            r["_row_no"], r["_assembly"], r["_biosample"],
            r["_country"], r["_host"], r["_year"], r["_genome_path"]
        ])


# ============================================================
# 6. Remove duplicate #BioSample, keeping first CSV occurrence
#
# Blank BioSample values are NOT collapsed together.
# ============================================================

seen_biosample = {}
selected_rows = []
dup_biosample = []

for row in assembly_unique:
    bs = row["_biosample"]

    if is_missing(bs):
        selected_rows.append(row)
        continue

    key = bs.casefold()

    if key in seen_biosample:
        row["_kept_assembly"] = seen_biosample[key]["_assembly"]
        dup_biosample.append(row)
        continue

    seen_biosample[key] = row
    selected_rows.append(row)

with open(dup_biosample_tsv, "w", newline="", encoding="utf-8") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow([
        "metadata_row", "Assembly", "#BioSample",
        "country", "host_1", "year", "kept_Assembly", "genome_path"
    ])
    for r in dup_biosample:
        w.writerow([
            r["_row_no"], r["_assembly"], r["_biosample"],
            r["_country"], r["_host"], r["_year"],
            r["_kept_assembly"], r["_genome_path"]
        ])


# ============================================================
# 7. Build selected genome lookup
# ============================================================

selected = {r["_genome_id"]: r for r in selected_rows}
selected_ids = set(selected)

with open(selected_tsv, "w", newline="", encoding="utf-8") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow([
        "candidate_db_id", "Assembly", "#BioSample",
        "country", "host_1", "year", "genome_path", "global_mash_degree"
    ])
    for r in sorted(selected_rows, key=lambda x: x["_genome_id"]):
        w.writerow([
            r["_genome_id"], r["_assembly"], r["_biosample"],
            r["_country"], r["_host"], r["_year"],
            r["_genome_path"], r["_degree"]
        ])

print(f"Matched metadata genomes:              {len(matched):,}", flush=True)
print(f"Excluded missing required stratum fields: {len(missing_meta):,}", flush=True)
print(f"Duplicate Assembly rows removed:      {len(dup_assembly):,}", flush=True)
print(f"Duplicate BioSamples removed:         {len(dup_biosample):,}", flush=True)
print(f"Assemblies not found in candidate DB: {len(unmatched):,}", flush=True)
print(f"Selected genomes for analysis:        {len(selected_rows):,}", flush=True)


# ============================================================
# 8. Open existing dnadiff cache READ-ONLY
# ============================================================

old_conn = None
if old_cache_database:
    old_conn = sqlite3.connect(
        "file:" + old_cache_database + "?mode=ro",
        uri=True,
        timeout=60,
    )
    old_conn.execute("PRAGMA busy_timeout=60000")
    old_conn.execute("PRAGMA cache_size=-500000")
    old_conn.execute("PRAGMA mmap_size=2147483648")


# ============================================================
# 9. Supplemental cache for only NEW same-stratum comparisons
# ============================================================

supp_conn = sqlite3.connect(supplemental_database, timeout=60)
supp_conn.execute("PRAGMA journal_mode=WAL")
supp_conn.execute("PRAGMA synchronous=NORMAL")
supp_conn.execute("PRAGMA cache_size=-250000")
supp_conn.execute(
    """
    CREATE TABLE IF NOT EXISTS comparisons
    (
        g1 TEXT NOT NULL,
        g2 TEXT NOT NULL,
        mash REAL NOT NULL,
        ni REAL,
        cov1 REAL,
        cov2 REAL,
        status TEXT NOT NULL,
        error TEXT,
        PRIMARY KEY (g1, g2)
    )
    """
)
supp_conn.commit()


def fetch_cache(conn, genome1, genome2):
    g1, g2 = canonical_pair(genome1, genome2)
    row = conn.execute(
        """
        SELECT mash, ni, cov1, cov2, status, error
        FROM comparisons
        WHERE g1=? AND g2=?
        """,
        (g1, g2),
    ).fetchone()

    if row is None:
        return None

    mash, ni, cov1, cov2, status, error = row

    if status != "OK":
        return None

    return {
        "g1": g1, "g2": g2,
        "mash": mash,
        "ni": ni,
        "cov1": cov1,
        "cov2": cov2,
        "status": status,
        "error": error or "",
    }


def get_cached(genome1, genome2):
    # Supplemental cache first, then the already completed global cache.
    r = fetch_cache(supp_conn, genome1, genome2)
    if r is not None:
        return r, "supp_cache"

    if old_conn is not None:
        r = fetch_cache(old_conn, genome1, genome2)
        if r is not None:
            return r, "old_cache"

    return None, None


def save_supp(result):
    supp_conn.execute(
        """
        INSERT OR REPLACE INTO comparisons
        (g1, g2, mash, ni, cov1, cov2, status, error)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            result["g1"], result["g2"], result["mash"],
            result["ni"], result["cov1"], result["cov2"],
            result["status"], result["error"],
        ),
    )


# ============================================================
# 10. Candidate lookup
#
# IMPORTANT:
#   Candidate DB already contains Mash <= 0.001 pairs.
#   Cross-country / cross-host / cross-year neighbors are
#   discarded BEFORE any dnadiff lookup or new dnadiff run.
# ============================================================

def get_same_stratum_candidates(rep_idx, assigned):
    rep_meta = selected[rep_idx]
    rep_stratum = rep_meta["_stratum"]

    cursor = candidate_conn.execute(
        """
        SELECT neighbor, mash
        FROM
        (
            SELECT b AS neighbor, mash
            FROM pairs
            WHERE a=?

            UNION ALL

            SELECT a AS neighbor, mash
            FROM pairs
            WHERE b=?
        )
        ORDER BY mash ASC
        """,
        (rep_idx, rep_idx),
    )

    for neighbor, mash in cursor:
        if mash > MASH_CUTOFF:
            continue
        if neighbor not in selected_ids:
            continue
        if assigned.get(neighbor, -1) != -1:
            continue
        if selected[neighbor]["_stratum"] != rep_stratum:
            continue
        yield neighbor, mash


# ============================================================
# 11. Conservative representative-anchored clustering
#
# Uses the same low global Mash-degree first order as the
# existing global script, but clustering is restricted to
# country + harmonized host group (host_1) + year strata.
# ============================================================

representative_order = sorted(
    selected_ids,
    key=lambda i: (id_to_degree[i], id_to_path[i])
)

assigned = {i: -1 for i in selected_ids}
representatives = []
member_metrics = {}

old_cache_hits = 0
supp_cache_hits = 0
new_dnadiff_count = 0
error_count = 0
cross_stratum_skipped = 0

with open(error_log, "w", encoding="utf-8") as handle:
    handle.write("genome1\tgenome2\terror\n")

executor = ThreadPoolExecutor(max_workers=MAX_WORKERS)

try:
    for rep_idx in representative_order:

        if assigned[rep_idx] != -1:
            continue

        assigned[rep_idx] = rep_idx
        representatives.append(rep_idx)

        rep_path = id_to_path[rep_idx]

        candidates = list(get_same_stratum_candidates(rep_idx, assigned))

        if not candidates:
            if len(representatives) % 100 == 0:
                assigned_n = sum(v != -1 for v in assigned.values())
                print(
                    f"Representatives: {len(representatives):,} | "
                    f"assigned: {assigned_n:,}/{len(selected_ids):,} | "
                    f"new dnadiff: {new_dnadiff_count:,} | "
                    f"old cache: {old_cache_hits:,} | "
                    f"supp cache: {supp_cache_hits:,}",
                    flush=True,
                )
            continue

        missing = []

        # Reuse cached comparisons first.
        for candidate_idx, mash_distance in candidates:
            if assigned[candidate_idx] != -1:
                continue

            candidate_path = id_to_path[candidate_idx]
            cached, source = get_cached(rep_path, candidate_path)

            if cached is None:
                missing.append((candidate_idx, mash_distance))
                continue

            if source == "old_cache":
                old_cache_hits += 1
            else:
                supp_cache_hits += 1

            if passes_threshold(cached["ni"], cached["cov1"], cached["cov2"]):
                if cached["g1"] == rep_path:
                    cov_rep = cached["cov1"]
                    cov_member = cached["cov2"]
                else:
                    cov_rep = cached["cov2"]
                    cov_member = cached["cov1"]

                assigned[candidate_idx] = rep_idx
                member_metrics[candidate_idx] = (
                    mash_distance,
                    cached["ni"],
                    cov_rep,
                    cov_member,
                    source,
                )

        # Run only missing same-country + same-host-group + same-year pairs.
        for start in range(0, len(missing), DNADIFF_BATCH):
            chunk = missing[start:start + DNADIFF_BATCH]
            future_map = {}

            for candidate_idx, mash_distance in chunk:
                if assigned[candidate_idx] != -1:
                    continue

                candidate_path = id_to_path[candidate_idx]
                future = executor.submit(
                    run_dnadiff,
                    rep_path,
                    candidate_path,
                    mash_distance,
                )
                future_map[future] = (candidate_idx, mash_distance)

            for future in as_completed(future_map):
                candidate_idx, mash_distance = future_map[future]
                result = future.result()
                new_dnadiff_count += 1

                save_supp(result)

                if new_dnadiff_count % 100 == 0:
                    supp_conn.commit()

                if result["status"] != "OK":
                    error_count += 1
                    with open(error_log, "a", encoding="utf-8") as handle:
                        err = (
                            result["error"]
                            .replace("\t", " ")
                            .replace("\n", " ")
                        )
                        handle.write(
                            f'{result["g1"]}\t{result["g2"]}\t{err}\n'
                        )
                    continue

                if passes_threshold(
                    result["ni"],
                    result["cov1"],
                    result["cov2"],
                ):
                    if result["g1"] == rep_path:
                        cov_rep = result["cov1"]
                        cov_member = result["cov2"]
                    else:
                        cov_rep = result["cov2"]
                        cov_member = result["cov1"]

                    assigned[candidate_idx] = rep_idx
                    member_metrics[candidate_idx] = (
                        mash_distance,
                        result["ni"],
                        cov_rep,
                        cov_member,
                        "new_dnadiff",
                    )

                if new_dnadiff_count % 5000 == 0:
                    assigned_n = sum(v != -1 for v in assigned.values())
                    print(
                        f"dnadiff progress | new: {new_dnadiff_count:,} | "
                        f"assigned: {assigned_n:,}/{len(selected_ids):,} | "
                        f"representatives: {len(representatives):,} | "
                        f"old cache: {old_cache_hits:,} | "
                        f"supp cache: {supp_cache_hits:,}",
                        flush=True,
                    )

            supp_conn.commit()

        if len(representatives) % 100 == 0:
            assigned_n = sum(v != -1 for v in assigned.values())
            print(
                f"Representatives: {len(representatives):,} | "
                f"assigned: {assigned_n:,}/{len(selected_ids):,} | "
                f"new dnadiff: {new_dnadiff_count:,} | "
                f"old cache: {old_cache_hits:,} | "
                f"supp cache: {supp_cache_hits:,}",
                flush=True,
            )

finally:
    executor.shutdown(wait=True)
    supp_conn.commit()


# Safety: any still-unassigned selected genome becomes its own representative.
for i in representative_order:
    if assigned[i] == -1:
        assigned[i] = i
        representatives.append(i)


# ============================================================
# 12. Final cluster outputs
# ============================================================

cluster_id = {
    rep_idx: n
    for n, rep_idx in enumerate(representatives, start=1)
}

cluster_size = {rep_idx: 0 for rep_idx in representatives}
for i in selected_ids:
    cluster_size[assigned[i]] += 1

with open(representative_txt, "w", encoding="utf-8") as handle:
    for rep_idx in representatives:
        handle.write(id_to_path[rep_idx] + "\n")

with open(removed_txt, "w", encoding="utf-8") as handle:
    for i in sorted(selected_ids):
        if assigned[i] != i:
            handle.write(id_to_path[i] + "\n")

with open(representative_tsv, "w", newline="", encoding="utf-8") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow([
        "cluster_id", "Assembly", "#BioSample",
        "country", "host_1", "year",
        "cluster_size", "genome_path"
    ])
    for rep_idx in representatives:
        r = selected[rep_idx]
        w.writerow([
            f"Cluster_{cluster_id[rep_idx]}",
            r["_assembly"],
            r["_biosample"],
            r["_country"],
            r["_host"],
            r["_year"],
            cluster_size[rep_idx],
            id_to_path[rep_idx],
        ])

with open(cluster_tsv, "w", newline="", encoding="utf-8") as handle:
    w = csv.writer(handle, delimiter="\t", lineterminator="\n")
    w.writerow([
        "cluster_id",
        "representative",
        "representative_Assembly",
        "genome",
        "Assembly",
        "#BioSample",
        "country",
        "host_1",
        "year",
        "cluster_size",
        "mash_distance_to_rep",
        "NI_percent",
        "coverage_rep_percent",
        "coverage_genome_percent",
        "comparison_source",
    ])

    for i in sorted(selected_ids):
        rep_idx = assigned[i]
        r = selected[i]
        rep_meta = selected[rep_idx]

        if i == rep_idx:
            w.writerow([
                f"Cluster_{cluster_id[rep_idx]}",
                id_to_path[rep_idx],
                rep_meta["_assembly"],
                id_to_path[i],
                r["_assembly"],
                r["_biosample"],
                r["_country"],
                r["_host"],
                r["_year"],
                cluster_size[rep_idx],
                "NA", "NA", "NA", "NA", "representative",
            ])
        else:
            mash, ni, cov_rep, cov_member, source = member_metrics[i]
            w.writerow([
                f"Cluster_{cluster_id[rep_idx]}",
                id_to_path[rep_idx],
                rep_meta["_assembly"],
                id_to_path[i],
                r["_assembly"],
                r["_biosample"],
                r["_country"],
                r["_host"],
                r["_year"],
                cluster_size[rep_idx],
                f"{mash:.8f}",
                f"{ni:.5f}",
                f"{cov_rep:.5f}",
                f"{cov_member:.5f}",
                source,
            ])


# ============================================================
# 13. Summary
# ============================================================

singleton_count = sum(v == 1 for v in cluster_size.values())
multi_count = sum(v > 1 for v in cluster_size.values())
largest_cluster = max(cluster_size.values()) if cluster_size else 0
removed_count = len(selected_ids) - len(representatives)
elapsed_h = (time.time() - START) / 3600.0

with open(summary_txt, "w", encoding="utf-8") as handle:
    handle.write("A. baumannii country-host-group-year stratified genomic de-replication\n")
    handle.write("=========================================================\n")
    handle.write(f"Metadata rows:\t{len(rows)}\n")
    handle.write(f"Matched metadata genomes:\t{len(matched)}\n")
    handle.write(f"Assemblies not found:\t{len(unmatched)}\n")
    missing_country_n = sum("country" in r.get("_missing_reason", "").split(",") for r in missing_meta)
    missing_hostgroup_n = sum("host_1" in r.get("_missing_reason", "").split(",") for r in missing_meta)
    missing_year_n = sum("year" in r.get("_missing_reason", "").split(",") for r in missing_meta)
    handle.write(f"Excluded missing required stratum fields:\t{len(missing_meta)}\n")
    handle.write(f"  missing country (including combinations):\t{missing_country_n}\n")
    handle.write(f"  missing harmonized host group (including combinations):\t{missing_hostgroup_n}\n")
    handle.write(f"  missing year (including combinations):\t{missing_year_n}\n")
    handle.write(f"Duplicate Assembly rows removed:\t{len(dup_assembly)}\n")
    handle.write(f"Duplicate BioSamples removed:\t{len(dup_biosample)}\n")
    handle.write(f"Selected genomes analyzed:\t{len(selected_ids)}\n")
    handle.write("Stratification:\tcountry + harmonized host group (host_1) + year\n")
    handle.write("Clustering:\trepresentative-anchored; no single-linkage propagation\n")
    handle.write("Representative order:\tascending global Mash degree, then genome path\n")
    handle.write(f"Mash candidate cutoff:\t{MASH_CUTOFF}\n")
    handle.write(f"NI cutoff:\t{NI_CUTOFF}%\n")
    handle.write(f"Coverage cutoff:\t{COV_CUTOFF}% BOTH genomes\n")
    handle.write(f"Representatives retained:\t{len(representatives)}\n")
    handle.write(f"Genomes collapsed/removed:\t{removed_count}\n")
    handle.write(f"Singleton clusters:\t{singleton_count}\n")
    handle.write(f"Clusters with >1 genome:\t{multi_count}\n")
    handle.write(f"Largest cluster:\t{largest_cluster}\n")
    handle.write(f"Old cache hits:\t{old_cache_hits}\n")
    handle.write(f"Supplemental cache hits:\t{supp_cache_hits}\n")
    handle.write(f"New dnadiff comparisons:\t{new_dnadiff_count}\n")
    handle.write(f"dnadiff errors:\t{error_count}\n")
    handle.write(f"Elapsed hours:\t{elapsed_h:.3f}\n")

print()
print("============================================================")
print("FINISHED")
print("============================================================")
print(f"Selected genomes:        {len(selected_ids):,}")
print(f"Representatives retained:{len(representatives):,}")
print(f"Collapsed/removed:       {removed_count:,}")
print(f"Old cache hits:          {old_cache_hits:,}")
print(f"Supplemental cache hits: {supp_cache_hits:,}")
print(f"New dnadiff comparisons: {new_dnadiff_count:,}")
print(f"dnadiff errors:          {error_count:,}")
print(f"Elapsed hours:           {elapsed_h:.2f}")
print()
print("Outputs:")
print(cluster_tsv)
print(representative_tsv)
print(summary_txt)

supp_conn.close()
if old_conn is not None:
    old_conn.close()
candidate_conn.close()
PY

echo
echo "============================================================"
echo "Output files"
echo "============================================================"
echo "Selected genomes:"
echo "${SELECTED_TSV}"
echo
echo "Excluded missing metadata:"
echo "${MISSING_TSV}"
echo
echo "Duplicate BioSample rows removed:"
echo "${DUP_BIOSAMPLE_TSV}"
echo
echo "Metadata Assembly not found:"
echo "${UNMATCHED_TSV}"
echo
echo "Clusters:"
echo "${CLUSTER_TSV}"
echo
echo "Representatives:"
echo "${REP_TSV}"
echo
echo "Summary:"
echo "${SUMMARY_TXT}"
echo
echo "Finished: $(date)"
