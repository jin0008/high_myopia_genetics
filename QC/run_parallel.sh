#!/usr/bin/env bash
# 오류 발생 시 중단하되, 파이프 오류도 체크하도록 설정
set -uo pipefail

###############################################################################
#  run_parallel.sh — GNU Parallel을 이용한 샘플별 병렬 QC 실행 (수정본)
###############################################################################

PROJECT_DIR=""
BED_DIR=""
JOBS=10            # 동시 실행 샘플 수
THREADS=10         # 샘플당 thread 수
REF_GENOME="/media/hanjinu/PM883/db/refs/hg38_broad/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR=""
SAMPLE_LIST=""    # 지정하지 않으면 자동 생성

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  -p, --project-dir  DIR    high_myopia 디렉토리 경로
  -b, --bed-dir      DIR    capture kit BED 파일 디렉토리

Optional:
  -j, --jobs         INT    동시 실행 샘플 수 (default: 4)
  -t, --threads      INT    샘플당 thread 수 (default: 4)
  -r, --reference    FILE   Reference FASTA
  -o, --output-dir   DIR    결과 디렉토리 (default: project_dir/qc_results)
  -s, --sample-list  FILE   sample_list.txt (없으면 자동 생성)
  -h, --help                도움말

EOF
    exit 0
}

# 인자 파싱
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--project-dir)  PROJECT_DIR="${2%/}"; shift 2;;
        -b|--bed-dir)      BED_DIR="${2%/}";     shift 2;;
        -j|--jobs)         JOBS="$2";         shift 2;;
        -t|--threads)      THREADS="$2";      shift 2;;
        -r|--reference)    REF_GENOME="$2";   shift 2;;
        -o|--output-dir)   OUTPUT_DIR="${2%/}";   shift 2;;
        -s|--sample-list)  SAMPLE_LIST="$2";  shift 2;;
        -h|--help)         usage;;
        *) echo "Unknown: $1"; usage;;
    esac
done

# 필수 인자 확인
[[ -z "$PROJECT_DIR" ]] && { echo "ERROR: --project-dir required"; usage; }
[[ -z "$BED_DIR" ]]     && { echo "ERROR: --bed-dir required"; usage; }
[[ -z "$OUTPUT_DIR" ]]  && OUTPUT_DIR="${PROJECT_DIR}/qc_results"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─────────────────── tool check ───────────────────
command -v parallel &>/dev/null || { echo "ERROR: GNU parallel not found."; exit 1; }
for tool in samtools mosdepth bedtools python3; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not found"; exit 1; }
done

# ─────────────────── sample_list.txt 생성 ───────────────────
mkdir -p "$OUTPUT_DIR"

if [[ -z "$SAMPLE_LIST" ]]; then
    SAMPLE_LIST="${OUTPUT_DIR}/sample_list.txt"
    echo "=== Generating sample_list.txt ==="
    > "$SAMPLE_LIST"
    
    count=0
    shopt -s nullglob
    for d in "${PROJECT_DIR}"/*; do
        [[ -d "$d" ]] || continue
        
        # 경로 정규화 및 샘플명 추출
        current_dir="${d%/}"
        sample=$(basename "$current_dir")
        bam="${current_dir}/${sample}.dedup.bam"
        
        if [[ -f "$bam" ]]; then
            echo "$sample" >> "$SAMPLE_LIST"
            ((count++))
        fi
    done
    shopt -u nullglob
    echo "Found ${count} samples → ${SAMPLE_LIST}"
else
    [[ ! -f "$SAMPLE_LIST" ]] && { echo "ERROR: Sample list not found: ${SAMPLE_LIST}"; exit 1; }
    count=$(wc -l < "$SAMPLE_LIST")
    echo "Using provided sample list: ${SAMPLE_LIST} (${count} samples)"
fi

[[ $count -eq 0 ]] && { echo "ERROR: No samples to process."; exit 1; }

# 환경 변수 export (parallel에서 사용)
export PROJECT_DIR BED_DIR OUTPUT_DIR THREADS REF_GENOME SCRIPT_DIR

# ─────────────────── 샘플별 실행 함수 ───────────────────
run_sample_qc() {
    local sample="$1"
    # 경로 재구성 시 오타 방지를 위해 명시적 선언
    local sample_dir="${PROJECT_DIR}/${sample}"
    local bam="${sample_dir}/${sample}.dedup.bam"
    
    local p1_dir="${OUTPUT_DIR}/phase1_bam_qc/${sample}"
    local p2_dir="${OUTPUT_DIR}/phase2_kit_inference/${sample}"
    local p3_dir="${OUTPUT_DIR}/phase3_coverage/${sample}"

    mkdir -p "$p1_dir" "$p2_dir" "$p3_dir"

    echo "[$(date '+%H:%M:%S')] START ${sample}"

    # BAM 파일 존재 여부 재확인
    if [[ ! -f "$bam" ]]; then
        echo "  ✗ ${sample}: BAM not found at ${bam}"
        return 1
    fi

    # ── BAM index ──
    if [[ ! -f "${bam}.bai" ]] && [[ ! -f "${bam%.bam}.bai" ]]; then
        samtools index -@ "$THREADS" "$bam"
    fi

    # ── Phase 1: BAM QC ──
    [[ ! -f "${p1_dir}/flagstat.txt" ]] && samtools flagstat -@ "$THREADS" "$bam" > "${p1_dir}/flagstat.txt"
    [[ ! -f "${p1_dir}/samtools_stats.txt" ]] && samtools stats -@ "$THREADS" "$bam" > "${p1_dir}/samtools_stats.txt"
    [[ ! -f "${p1_dir}/idxstats.txt" ]] && samtools idxstats "$bam" > "${p1_dir}/idxstats.txt"
    [[ ! -f "${p1_dir}/${sample}.mosdepth.summary.txt" ]] && mosdepth --no-per-base --threads "$THREADS" --fast-mode "${p1_dir}/${sample}" "$bam"

    # ── Phase 2: Kit inference ──
    local total_reads
    total_reads=$(samtools view -c -F 0x904 -@ "$THREADS" "$bam")

    local best_kit="UNKNOWN"
    local best_pct=0
    > "${p2_dir}/kit_comparison.tsv"

    for bedfile in "${BED_DIR}"/*.bed; do
        [[ -f "$bedfile" ]] || continue
        local kit_name=$(basename "$bedfile" .bed)
        local on_target_reads=$(samtools view -c -F 0x904 -@ "$THREADS" -L "$bedfile" "$bam")

        local on_target_pct="0"
        if [[ "$total_reads" -gt 0 ]]; then
            on_target_pct=$(awk "BEGIN{printf \"%.4f\", ${on_target_reads}/${total_reads}*100}")
        fi

        echo -e "${kit_name}\t${on_target_reads}\t${total_reads}\t${on_target_pct}%" >> "${p2_dir}/kit_comparison.tsv"

        if awk "BEGIN{exit !(${on_target_pct} > ${best_pct})}"; then
            best_pct="$on_target_pct"
            best_kit="$kit_name"
        fi
    done

    echo -e "${sample}\t${best_kit}\t${best_pct}%" > "${p2_dir}/kit_result.tsv"

    # ── Phase 3: Coverage QC ──
    local target_bed="${BED_DIR}/${best_kit}.bed"
    [[ ! -f "$target_bed" ]] && target_bed=$(find "$BED_DIR" -name "*.bed" | head -1)

    if [[ -n "$target_bed" ]]; then
        mosdepth --by "$target_bed" --no-per-base --threads "$THREADS" \
                 --thresholds 1,5,10,15,20,30,50,100 \
                 "${p3_dir}/${sample}_target" "$bam"
    fi

    echo "[$(date '+%H:%M:%S')] DONE  ${sample} (Kit: ${best_kit})"
}
export -f run_sample_qc

# ─────────────────── 실행 ───────────────────
echo "=== Starting Parallel Execution ==="
mkdir -p "${OUTPUT_DIR}/logs"

parallel --jobs "$JOBS" \
         --bar \
         --joblog "${OUTPUT_DIR}/logs/parallel_joblog.txt" \
         run_sample_qc \
         :::: "$SAMPLE_LIST"

# 결과 병합 및 리포트 (Phase 4)
echo "=== Summary & Finalizing ==="
echo -e "sample\tinferred_kit\ton_target_pct" > "${OUTPUT_DIR}/phase2_kit_inference/kit_inference_summary.tsv"
while IFS= read -r s; do
    [[ -f "${OUTPUT_DIR}/phase2_kit_inference/${s}/kit_result.tsv" ]] && cat "${OUTPUT_DIR}/phase2_kit_inference/${s}/kit_result.tsv" >> "${OUTPUT_DIR}/phase2_kit_inference/kit_inference_summary.tsv"
done < "$SAMPLE_LIST"

# Python 리포트 생성 (파일이 존재할 경우에만 실행)
if [[ -f "${SCRIPT_DIR}/summarize_qc.py" ]]; then
    python3 "${SCRIPT_DIR}/summarize_qc.py" --project-dir "$PROJECT_DIR" --output-dir "$OUTPUT_DIR"
fi

echo "Done. Results in ${OUTPUT_DIR}"
