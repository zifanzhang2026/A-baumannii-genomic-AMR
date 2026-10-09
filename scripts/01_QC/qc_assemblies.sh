#!/usr/bin/env bash
set -euo pipefail

# ========= 配置 =========
ASM_DIR="${ASM_DIR:-wsh}"       # assemblies 目录（支持 .fna/.fa/.fasta 及 .gz）
OUTDIR="${OUTDIR:-O}"              # QC 输出目录
THREADS="${THREADS:-32}"                 # 并发
MIN_LEN="${MIN_LEN:-3600000}"            # 总长度下限（示例：金葡 2.6 Mb）
MAX_LEN="${MAX_LEN:-4300000}"            # 总长度上限（示例：金葡 3.2 Mb）
MAX_CONTIGS="${MAX_CONTIGS:-300}"        # contig 数上限
MAX_N_FRAC="${MAX_N_FRAC:-0.05}"         # N 比例上限（1%）
DO_ANI="${DO_ANI:-1}"                    # 是否做 fastANI（0/1）
REF_FOR_ANI="${REF_FOR_ANI:-ref.fa}"     # ANI 参考（DO_ANI=1 时使用）
# ========================

mkdir -p "$OUTDIR" "$OUTDIR/logs"

# 0) 依赖自检
need() { command -v "$1" >/dev/null 2>&1 || { echo "[ERR] 缺少 $1，请安装后再跑"; exit 1; }; }
need awk; need parallel
if ! command -v seqkit >/dev/null 2>&1; then
  echo "[ERR] 缺少 seqkit（用于装配统计）。conda 安装：conda install -c bioconda seqkit"
  exit 1
fi
if [[ "$DO_ANI" == "1" ]] && ! command -v fastANI >/dev/null 2>&1; then
  echo "[ERR] DO_ANI=1 但未安装 fastANI。conda 安装：conda install -c bioconda fastani"; exit 1
fi

# 1) 列出所有候选文件（绝对路径）
find "$(realpath "$ASM_DIR")" -maxdepth 1 -type f \
  \( -iname "*.fna" -o -iname "*.fna.gz" -o -iname "*.fa" -o -iname "*.fa.gz" -o -iname "*.fasta" -o -iname "*.fasta.gz" \) \
  | sort > "$OUTDIR/asm_all.txt"
echo "[INFO] 总样本数：$(wc -l < "$OUTDIR/asm_all.txt")"

# 2) gzip 完整性检查（仅 .gz）
echo "[INFO] 开始 gzip 完整性检查（仅 .gz）..."
awk '/\.gz$/ {print}' "$OUTDIR/asm_all.txt" > "$OUTDIR/gz_list.txt" || true
: > "$OUTDIR/bad_gzip.txt"
if [[ -s "$OUTDIR/gz_list.txt" ]]; then
  parallel -j "$THREADS" 'gzip -t {} || echo {}' :::: "$OUTDIR/gz_list.txt" \
    | tee "$OUTDIR/bad_gzip.txt" >/dev/null
fi
GOOD_LIST="$OUTDIR/asm_good_gzip.txt"
if [[ -s "$OUTDIR/bad_gzip.txt" ]]; then
  grep -Fxv -f "$OUTDIR/bad_gzip.txt" "$OUTDIR/asm_all.txt" > "$GOOD_LIST" || true
  echo "[WARN] 损坏 gzip 文件数：$(wc -l < "$OUTDIR/bad_gzip.txt")"
else
  cp "$OUTDIR/asm_all.txt" "$GOOD_LIST"
fi
echo "[INFO] gzip 通过的样本数：$(wc -l < "$GOOD_LIST")"

# 3) seqkit 统计（长度/contig/N 比例）
echo "[INFO] 运行 seqkit stats ..."
seqkit stats -a $(cat "$GOOD_LIST") > "$OUTDIR/seqkit_stats.tsv"

# 4) 基于门槛筛选 PASS/FAIL
# seqkit 列（>=2.5.1）：file format type num_seqs sum_len min len avg len max len Q1 Q2 Q3 sum_gap N50 Q20(%) Q30(%) GC(%)
# 其中 sum_gap 是 N 的总量（大写 N），N 比例 = sum_gap / sum_len
awk -v MIN_LEN="$MIN_LEN" -v MAX_LEN="$MAX_LEN" -v MAX_CONTIGS="$MAX_CONTIGS" -v MAX_N_FRAC="$MAX_N_FRAC" '
BEGIN{OFS="\t"}
NR==1{print $0,"\tPASS_FAIL","\tREASON"; next}
{
  file=$1; num=$4; 
  
  # 关键修复：去掉数字中的逗号
  len_str = $5; gsub(/,/, "", len_str); len = len_str + 0
  gaps_str = $(NF-5); gsub(/,/, "", gaps_str); gaps = gaps_str + 0
  
  pass=1; reason=""
  if (len<MIN_LEN || len>MAX_LEN){pass=0; reason=reason"len_out_of_range; "}
  if (num>MAX_CONTIGS){pass=0; reason=reason"too_many_contigs; "}
  if (gaps/len>MAX_N_FRAC){pass=0; reason=reason"high_N_fraction; "}
  if (pass==1){pf="PASS"; reason="."} else {pf="FAIL"}
  print $0, pf, reason
}' "$OUTDIR/seqkit_stats.tsv" > "$OUTDIR/seqkit_stats.flag.tsv"

awk '$0 ~ /^#/ {next} NR>1 && $0 !~ /FAIL/ {print $1}' "$OUTDIR/seqkit_stats.flag.tsv" > "$OUTDIR/pass_by_seqkit.txt"
awk '$0 ~ /^#/ {next} NR>1 && $0 ~ /FAIL/ {print $1}' "$OUTDIR/seqkit_stats.flag.tsv" > "$OUTDIR/fail_by_seqkit.txt"

# 5) （可选）fastANI 物种/相近度快速筛（默认关闭）
if [[ "$DO_ANI" == "1" ]]; then
  echo "[INFO] 运行 fastANI（对 seqkit PASS 的样本）..."
  mkdir -p "$OUTDIR/ani"
  parallel -j "$THREADS" 'fastANI -q {} -r "'"$REF_FOR_ANI"'" -o "'"$OUTDIR"'/ani/{/}.tsv" 2> /dev/null || true' :::: "$OUTDIR/pass_by_seqkit.txt"
  # 聚合 ANI 结果（fastANI 输出列：query ref ANI fragments matchedFragments）
  echo -e "file\tANI" > "$OUTDIR/ani_summary.tsv"
  for f in "$OUTDIR"/ani/*.tsv; do
    [[ -s "$f" ]] || continue
    ani=$(awk 'NR==1{print $3}' "$f")
    q=$(awk 'NR==1{print $1}' "$f")
    echo -e "$q\t$ani"
  done >> "$OUTDIR/ani_summary.tsv"
  # 按阈值 98.0%（可改）筛
  ANI_MIN="${ANI_MIN:-95.0}"
  awk -v T="$ANI_MIN" 'NR==1{next} {if($2+0 >= T) print $1}' "$OUTDIR/ani_summary.tsv" > "$OUTDIR/pass_by_ani.txt"
  # 交集：seqkit PASS ∩ ANI PASS
  grep -Fx -f "$OUTDIR/pass_by_ani.txt" "$OUTDIR/pass_by_seqkit.txt" > "$OUTDIR/pass_final.txt" || true
else
  cp "$OUTDIR/pass_by_seqkit.txt" "$OUTDIR/pass_final.txt"
fi

# 6) 生成最终清单（供后续 SNP 流水线使用）
# 去掉 gzip 坏包 + 基础装配筛选（+ 可选 ANI）的最终 PASS 列表
echo "[INFO] 最终 PASS 样本数：$(wc -l < "$OUTDIR/pass_final.txt")"
cp "$OUTDIR/pass_final.txt" "$OUTDIR/pass_list.txt"

# 7) 汇总报告
echo "====== QC SUMMARY ======" | tee "$OUTDIR/summary.txt"
echo "Total files:            $(wc -l < "$OUTDIR/asm_all.txt")" | tee -a "$OUTDIR/summary.txt"
echo "Bad gzip files:         $( [[ -s $OUTDIR/bad_gzip.txt ]] && wc -l < $OUTDIR/bad_gzip.txt || echo 0 )" | tee -a "$OUTDIR/summary.txt"
echo "Seqkit PASS:            $(wc -l < "$OUTDIR/pass_by_seqkit.txt")" | tee -a "$OUTDIR/summary.txt"
if [[ "$DO_ANI" == "1" ]]; then
  echo "ANI PASS (>=${ANI_MIN}%): $(wc -l < "$OUTDIR/pass_by_ani.txt")" | tee -a "$OUTDIR/summary.txt"
fi
echo "PASS final (for SNP):   $(wc -l < "$OUTDIR/pass_list.txt")" | tee -a "$OUTDIR/summary.txt"
echo "Outputs in:             $OUTDIR" | tee -a "$OUTDIR/summary.txt"

