#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_scene_split.py —— 按“场景（scene）”固定划分 训练/验证/测试 的生成器

设计原则（对应第 2 周任务“按场景固定训练、验证和测试划分”）：
  1. 以“场景（一次独立采集）”为最小划分单元，同一景的所有裁剪块/1D 剖面
     必须落在同一集合，杜绝跨集泄漏（patch-level 随机划分会泄漏）。
  2. 确定性 + 可追加（append-only）：
     - 已在 scene_inventory.csv 的 assigned_split 列显式指定的，原样采用；
     - 未指定的，用 scene_id 的 SHA-256 哈希按固定比例桶分配，结果完全可复现；
     - 新增场景时，旧场景的分配永远不会被改动（只在尾部追加新 scene）。
  3. 纯标准库实现，不依赖 torch，可在任意 Python 3.8+ 运行。

用法：
  python scripts/make_scene_split.py                 # 由 inventory 生成 scene_split.json
  python scripts/make_scene_split.py --check         # 仅做完整性/泄漏检查，不写文件
  python scripts/make_scene_split.py --ratio 70,15,15 --seed 2020

输出：
  data/splits/scene_split.json  —— 权威划分文件（被 verify_io_and_split.py 读取）
"""
import argparse
import csv
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

# ---- 路径（相对本脚本定位项目根） ----
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DEFAULT_INVENTORY = os.path.join(PROJECT_ROOT, "data", "splits", "scene_inventory.csv")
DEFAULT_OUT = os.path.join(PROJECT_ROOT, "data", "splits", "scene_split.json")
SEED = 2020


def parse_ratio(s: str):
    parts = [int(x) for x in s.split(",")]
    if len(parts) != 3 or sum(parts) != 100:
        raise ValueError("ratio 必须是 3 个整数且和为 100，例如 70,15,15")
    return parts


def hash_bucket(scene_id: str, seed: int, ratio):
    """确定性哈希到 train/val/test 桶。"""
    h = hashlib.sha256(f"{seed}:{scene_id}".encode("utf-8")).hexdigest()
    v = int(h[:8], 16) % 100  # 0..99
    tr, va, _ = ratio
    if v < tr:
        return "train"
    if v < tr + va:
        return "val"
    return "test"


def load_inventory(path: str):
    """读取登记表：先滤过以 # 开头的注释行，再交给 DictReader（首行才是表头）。"""
    import io
    with open(path, "r", encoding="utf-8-sig") as f:  # utf-8-sig 自动去 BOM
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    rows = []
    reader = csv.DictReader(io.StringIO("".join(lines)))
    for row in reader:
        sid = (row.get("scene_id") or "").strip()
        if not sid:
            continue
        rows.append(row)
    return rows


def build_split(rows, ratio, seed, allow_hash_fill=True):
    splits = {"train": [], "val": [], "test": []}
    scene_meta = {}
    problems = []
    for r in rows:
        sid = r["scene_id"].strip()
        forced = (r.get("assigned_split") or "").strip().lower()
        if forced in ("train", "val", "test"):
            split = forced
        elif allow_hash_fill:
            split = hash_bucket(sid, seed, ratio)
        else:
            problems.append(f"{sid}: 未指定 assigned_split 且禁止哈希补全")
            split = None
        if split is None:
            continue
        splits[split].append(sid)
        scene_meta[sid] = {
            "platform": r.get("platform", ""),
            "sensor_mode": r.get("sensor_mode", ""),
            "acquisition_date": r.get("acquisition_date", ""),
            "products": r.get("products", ""),
            "polarization": r.get("polarization", ""),
            "source": r.get("source", ""),
            "status": r.get("status", ""),
            "split": split,
            "note": r.get("note", ""),
        }
    return splits, scene_meta, problems


def check_integrity(rows, splits):
    """检查：每个 scene 唯一分配、无重复、无跨集。"""
    errors = []
    seen = {}
    for split, ids in splits.items():
        for sid in ids:
            if sid in seen:
                errors.append(f"泄漏/重复：{sid} 同时出现在 {seen[sid]} 与 {split}")
            seen[sid] = split
    # inventory 与 split 一致性
    inv_ids = {r["scene_id"].strip() for r in rows}
    split_ids = set(seen.keys())
    missing = inv_ids - split_ids
    if missing:
        errors.append(f"以下 scene 未进入任何集合：{sorted(missing)}")
    extra = split_ids - inv_ids
    if extra:
        errors.append(f"以下 scene 不在登记表中却出现在划分：{sorted(extra)}")
    return errors


def main():
    ap = argparse.ArgumentParser(description="按场景固定划分 训练/验证/测试")
    ap.add_argument("--inventory", default=DEFAULT_INVENTORY)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--ratio", default="70,15,15", help="train,val,test 百分比，和为100")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--check", action="store_true", help="仅做完整性/泄漏检查，不写文件")
    ap.add_argument("--no-hash-fill", action="store_true", help="禁止对未指定场景哈希补全")
    args = ap.parse_args()

    ratio = parse_ratio(args.ratio)
    rows = load_inventory(args.inventory)
    splits, scene_meta, problems = build_split(rows, ratio, args.seed, allow_hash_fill=not args.no_hash_fill)

    errors = check_integrity(rows, splits)
    errors += problems

    n_train, n_val, n_test = len(splits["train"]), len(splits["val"]), len(splits["test"])
    print(f"[make_scene_split] 场景总数={len(rows)}  "
          f"train={n_train} val={n_val} test={n_test}  ratio={args.ratio}")

    if errors:
        print("[make_scene_split] ❌ 检查未通过：")
        for e in errors:
            print("   -", e)
        sys.exit(1)
    print("[make_scene_split] ✅ 完整性/泄漏检查通过：每景唯一分配，无跨集。")

    if args.check:
        return

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    doc = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "strategy": "scene-level fixed split; append-only; deterministic hash fallback",
        "ratio": {"train": ratio[0], "val": ratio[1], "test": ratio[2]},
        "seed": args.seed,
        "source_inventory": os.path.relpath(args.inventory, PROJECT_ROOT),
        "splits": splits,
        "scene_meta": scene_meta,
    }
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=2)
    print(f"[make_scene_split] ✅ 已写出 {args.out}")


if __name__ == "__main__":
    main()
