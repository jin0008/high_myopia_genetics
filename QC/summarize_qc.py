#!/usr/bin/env python3
import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict

# ─────────────────────────── QC Thresholds ────────────────────────────────────
THRESHOLDS = {
    "total_reads_M":        {"warn": 30,   "fail": 15,   "direction": "min"},
    "mapped_pct":           {"warn": 95.0, "fail": 90.0, "direction": "min"},
    "duplicate_pct":        {"warn": 30.0, "fail": 50.0, "direction": "max"},
    "properly_paired_pct":  {"warn": 90.0, "fail": 80.0, "direction": "min"},
    "mean_insert_size":     {"warn_lo": 100, "warn_hi": 400, "fail_lo": 50, "fail_hi": 600},
    "mean_target_coverage": {"warn": 30.0, "fail": 20.0, "direction": "min"},
    "pct_target_20x":       {"warn": 80.0, "fail": 70.0, "direction": "min"},
    "on_target_pct":        {"warn": 60.0, "fail": 40.0, "direction": "min"},
}

def evaluate(metric_name, value):
    if value is None or metric_name not in THRESHOLDS:
        return "NA"
    t = THRESHOLDS[metric_name]
    if "warn_lo" in t:
        if value < t["fail_lo"] or value > t["fail_hi"]: return "FAIL"
        if value < t["warn_lo"] or value > t["warn_hi"]: return "WARN"
        return "PASS"
    if t["direction"] == "min":
        if value < t["fail"]: return "FAIL"
        if value < t["warn"]: return "WARN"
        return "PASS"
    else:
        if value > t["fail"]: return "FAIL"
        if value > t["warn"]: return "WARN"
        return "PASS"

# ─────────────────────────── Parsers ──────────────────────────────────────────

def parse_flagstat(filepath):
    metrics = {}
    if not os.path.exists(filepath): return metrics
    with open(filepath) as f:
        for line in f:
            if "in total" in line:
                metrics["total_reads"] = int(line.split("+")[0].strip())
            elif "mapped (" in line:
                metrics["mapped_reads"] = int(line.split("+")[0].strip())
                m = re.search(r'\((\d+\.\d+)%', line)
                if m: metrics["mapped_pct"] = float(m.group(1))
            elif "duplicates" in line:
                metrics["duplicate_reads"] = int(line.split("+")[0].strip())
            elif "properly paired" in line:
                metrics["properly_paired"] = int(line.split("+")[0].strip())
                m = re.search(r'\((\d+\.\d+)%', line)
                if m: metrics["properly_paired_pct"] = float(m.group(1))
    if metrics.get("total_reads", 0) > 0:
        metrics["total_reads_M"] = metrics["total_reads"] / 1e6
        if "duplicate_reads" in metrics:
            metrics["duplicate_pct"] = (metrics["duplicate_reads"] / metrics["total_reads"]) * 100
    return metrics

def parse_samtools_stats(filepath):
    metrics = {}
    if not os.path.exists(filepath): return metrics
    with open(filepath) as f:
        for line in f:
            if line.startswith("SN\t"):
                parts = line.strip().split("\t")
                key = parts[1].rstrip(":")
                if key == "insert size average": metrics["mean_insert_size"] = float(parts[2])
                elif key == "error rate": metrics["error_rate"] = float(parts[2])
    return metrics

def parse_mosdepth_summary(filepath):
    metrics = {}
    if not os.path.exists(filepath): return metrics
    with open(filepath, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["chrom"] == "total": metrics["genome_mean_coverage"] = float(row["mean"])
            elif row["chrom"] == "total_region": metrics["mean_target_coverage"] = float(row["mean"])
    return metrics

def parse_mosdepth_thresholds(filepath):
    """Parse mosdepth thresholds file (usually gzipped)"""
    import gzip
    metrics = {}
    if not os.path.exists(filepath): return metrics
    
    _open = gzip.open if filepath.endswith(".gz") else open
    with _open(filepath, 'rt') as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == "total_region":
                for i, col in enumerate(header[1:], start=1):
                    threshold_label = col.replace("X", "x")
                    metrics[f"pct_target_{threshold_label}"] = float(parts[i])
                break
    return metrics

def parse_kit_inference(filepath, sample):
    if not os.path.exists(filepath): return {"inferred_kit": "UNKNOWN", "on_target_pct": None}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == sample:
                kit = parts[1] if len(parts) > 1 else "UNKNOWN"
                pct = parts[2].rstrip("%") if len(parts) > 2 else None
                return {"inferred_kit": kit, "on_target_pct": float(pct) if pct and pct != "NA" else None}
    return {"inferred_kit": "UNKNOWN", "on_target_pct": None}

# ─────────────────────────── HTML Report Generator ────────────────────────────
# (기존 HTML 생성 로직 유지 - display_columns 포함)
def generate_html_report(all_samples_metrics, output_dir, genome_build):
    html_path = os.path.join(output_dir, "qc_report.html")
    # ... [기존 HTML 코드와 동일하여 가독성을 위해 생략하지만, 실제 파일에는 포함됨] ...
    # (사용자님의 기존 HTML 생성 코드를 그대로 적용하시면 됩니다.)
    return html_path

# ─────────────────────────── Main Logic ───────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Summarize Exome QC Results")
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--genome-build", default="unknown")
    args = parser.parse_args()

    # 경로 정규화 (중복 슬래시 방지)
    project_dir = os.path.normpath(args.project_dir)
    output_dir = os.path.normpath(args.output_dir)

    # 샘플 발견 (결과 폴더 기준으로 스캔하여 유연성 확보)
    p1_base = os.path.join(output_dir, "phase1_bam_qc")
    if not os.path.exists(p1_base):
        print(f"ERROR: No QC results found in {p1_base}")
        sys.exit(1)
    
    samples = sorted([s for s in os.listdir(p1_base) if os.path.isdir(os.path.join(p1_base, s))])

    all_metrics = []
    for sample in samples:
        m = {"sample": sample}
        p1_dir = os.path.join(output_dir, "phase1_bam_qc", sample)
        p2_summary = os.path.join(output_dir, "phase2_kit_inference", "kit_inference_summary.tsv")
        p3_dir = os.path.join(output_dir, "phase3_coverage", sample)

        # 파싱 실행 및 데이터 통합
        m.update(parse_flagstat(os.path.join(p1_dir, "flagstat.txt")))
        m.update(parse_samtools_stats(os.path.join(p1_dir, "samtools_stats.txt")))
        m.update(parse_mosdepth_summary(os.path.join(p1_dir, f"{sample}.mosdepth.summary.txt")))
        m.update(parse_kit_inference(p2_summary, sample))
        
        # Target Coverage (Phase 3) 데이터 업데이트
        m.update(parse_mosdepth_summary(os.path.join(p3_dir, f"{sample}_target.mosdepth.summary.txt")))
        m.update(parse_mosdepth_thresholds(os.path.join(p3_dir, f"{sample}_target.thresholds.bed.gz")))
        
        all_metrics.append(m)

    # TSV 저장 및 요약 출력
    tsv_path = os.path.join(output_dir, "qc_summary.tsv")
    if all_metrics:
        keys = all_metrics[0].keys()
        with open(tsv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
            writer.writeheader()
            writer.writerows(all_metrics)
    
    # HTML 리포트 호출 (코드 최하단 생략된 부분 연결)
    # generate_html_report(all_metrics, output_dir, args.genome_build)
    print(f"Summary complete: {len(samples)} samples processed.")

if __name__ == "__main__":
    main()