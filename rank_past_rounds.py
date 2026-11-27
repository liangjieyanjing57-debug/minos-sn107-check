#!/usr/bin/env python3
"""Minos (SN107): 過去の本番ラウンドを手元で解き、「自分なら何位だったか」を出す。

使い方（minos_subnet のフォルダで、.venv を有効にしてから）:
    python rank_past_rounds.py                     # 直近の公開済みラウンド3回を標準設定(configs/gatk.conf)で
    python rank_past_rounds.py --rounds 5          # 5回分
    python rank_past_rounds.py --conf my.conf      # 別の設定ファイルで

チェーン・ウォレット・登録は一切使わない。採点は Minos 公式コードの
_run_and_score（バリデーターと同じ処理）を v2 方式で呼ぶ。
同点のときは「一番遅く提出した」とみなして順位を付ける（速さに頼らない前提）。
"""
import argparse
import csv
import json
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

PLATFORM = "https://api.theminos.ai"
SDF_FILES = ["done", "mainIndex", "nameIndex0", "namedata0", "namepointer0",
             "progress", "seqdata0", "seqpointer0", "sequenceIndex0", "summary.txt"]


def get_json(path):
    with urllib.request.urlopen(PLATFORM + path, timeout=60) as r:
        return json.load(r)


def fetch(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    subprocess.run(["curl", "-fsSL", "--max-time", "900", "-o", str(dest), url], check=True)


def ensure_reference(chrom, base):
    d = base / "datasets" / "reference" / chrom
    for ext in ("fa", "fa.fai", "dict"):
        fetch(f"{PLATFORM}/reference/{chrom}/{chrom}.{ext}", d / f"{chrom}.{ext}")
    sdf = d / f"{chrom}.sdf"
    for f in SDF_FILES:
        fetch(f"{PLATFORM}/reference/{chrom}/{chrom}.sdf/{f}", sdf / f)
    return d / f"{chrom}.fa", sdf


def dust_weights(cfg):
    """network-config の配分式: 1位=winner_weight、2位〜dust_top_n位で残りを dust_decay の等比で分ける。"""
    top, n, decay = cfg["winner_weight"], cfg["dust_top_n"], cfg["dust_decay"]
    rest = 1.0 - top - cfg.get("burn_rate", 0.0)
    terms = [decay ** i for i in range(n - 1)]
    s = sum(terms)
    return [top] + [rest * t / s for t in terms]


def score_of(e):
    for k in ("current_score", "combined_final"):
        v = e.get(k)
        if isinstance(v, (int, float)) and v > 0:
            return float(v)
    return None


def estimate_rank(my_score, entries):
    """同点は自分を最後尾に置く（提出が一番遅かった場合）。"""
    others = [s for s in (score_of(e) for e in entries) if s is not None]
    better_or_equal = sum(1 for s in others if s >= my_score - 1e-12)
    return better_or_equal + 1, sorted(others, reverse=True)


def revealed_round_ids(limit, skip=0):
    """確定済みの最新ラウンドから72分ずつさかのぼり、課題ファイルが公開済みのものを集める。"""
    from datetime import datetime, timedelta
    latest = get_json("/scoring/rounds/latest-finalized/leaderboard")["round"]["round_id"]
    t0 = datetime.fromisoformat(latest)
    out = []
    for i in range(skip, skip + limit * 4):
        rid = (t0 - timedelta(minutes=72 * i)).isoformat()
        q = urllib.parse.quote(rid, safe="")
        try:
            lb = get_json("/scoring/rounds/" + q + "/leaderboard")
            v = get_json("/verification/round/" + q)
        except Exception as e:
            print(f"  {rid}: 取得できず ({e})")
            continue
        if not lb.get("round", {}).get("is_finalized"):
            continue
        valid = sum(1 for e in lb.get("entries", []) if score_of(e) is not None)
        if valid < 20:
            print(f"  {rid}: 有効な点数が{valid}件しかないラウンドなので除外")
            continue
        if not v.get("reveal"):
            print(f"  {rid}: 課題ファイル未公開")
            continue
        lb["_reveal_files"] = v["reveal"].get("files", [])
        out.append((rid, lb))
        if len(out) >= limit:
            break
    print(f"対象ラウンド: {[r for r, _ in out]}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--conf", default="configs/gatk.conf")
    ap.add_argument("--workdir", default="past_rounds")
    ap.add_argument("--skip", type=int, default=0, help="新しい方から何ラウンド飛ばすか")
    args = ap.parse_args()

    base = Path.cwd()
    if not (base / "neurons" / "miner.py").exists():
        sys.exit("minos_subnet のフォルダで実行してください")
    sys.path.insert(0, str(base))
    from neurons import miner  # 公式の採点処理をそのまま使う
    from utils import scoring_version

    net = get_json("/scoring/network-config")
    version = net.get("scoring_version", "v2")
    scoring_version.record_used(version)
    weights = dust_weights(net)
    print(f"採点方式: {version} / 1位の取り分 {weights[0]:.3f} / 20位 {weights[-1]:.4f}")

    tool_config = miner._load_config_for_scoring(args.conf, "gatk")
    work = base / args.workdir
    results = []
    for rid, lb in revealed_round_ids(args.rounds, args.skip):
        region = lb["round"]["region"]
        chrom = region.split(":")[0]
        tag = rid.replace(":", "").replace("+0000", "")
        rdir = work / tag
        print(f"\n=== ラウンド {rid}  領域 {region} ===")
        allowed = {"input.bam", "input.bam.bai", "truth.vcf.gz", "truth.vcf.gz.tbi",
                   "mutations.vcf.gz", "mutations.vcf.gz.tbi"}
        for f in lb["_reveal_files"]:
            if f.get("name") in allowed and str(f.get("url", "")).startswith("https://"):
                fetch(f["url"], rdir / f["name"])
        if not (rdir / "input.bam").exists():
            print("  課題ファイルを取得できませんでした")
            results.append([rid, region, "", "", "", "", ""])
            continue
        try:
            ref, sdf = ensure_reference(chrom, base)
            score = miner._run_and_score(
                tool="gatk", tool_config=tool_config,
                bam_path=(rdir / "input.bam").resolve(), ref_path=ref.resolve(),
                truth_path=(rdir / "truth.vcf.gz").resolve(),
                mutations_path=(rdir / "mutations.vcf.gz").resolve(),
                region=region, reference_sdf=str(sdf.resolve()), scoring_version=version,
            )
        except Exception as e:  # 1ラウンドの失敗で全体を止めない
            print(f"  エラー: {e}")
            score = None
        if score is None:
            print("  採点できませんでした（上のエラーを確認）")
            results.append([rid, region, "", "", "", "", ""])
            continue
        rank, others = estimate_rank(score, lb["entries"])
        w = weights[rank - 1] if rank <= len(weights) else 0.0
        r20 = others[19] if len(others) >= 20 else None
        print(f"  自分の点数 {score:.6f} / 1位 {others[0]:.6f} / 20位 {r20 if r20 is None else round(r20, 6)}")
        print(f"  推定順位 {rank}位（同点は最後尾扱い） / 推定取り分 {w:.4f}")
        results.append([rid, region, f"{score:.6f}", f"{others[0]:.6f}",
                        "" if r20 is None else f"{r20:.6f}", rank, f"{w:.4f}"])

    out = work / "results.csv"
    work.mkdir(exist_ok=True)
    with open(out, "w", newline="") as f:
        csv.writer(f).writerows([["round", "region", "my_score", "top", "rank20", "est_rank", "est_weight"]] + results)
    scored = [r for r in results if r[6] != ""]
    if scored:
        avg_w = sum(float(r[6]) for r in scored) / len(scored)
        print(f"\n平均の推定取り分 {avg_w:.4f}（1日の鉱夫配分 約62 TAO なら 約{avg_w * 62:.2f} TAO/日）")
    print(f"結果: {out}")


if __name__ == "__main__":
    main()
