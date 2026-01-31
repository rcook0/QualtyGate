from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from typing import List, Tuple, Dict, Any

import numpy as np
import yaml

from qgate_v12 import (
    QGateV12Config, ClassifyCfg, BandsCfg, ADXCfg,
    compute_qgate_v12, align_htf_to_ltf, GateState, GateBand
)

def _parse_time(s: str) -> int:
    """
    Accepts:
      - ISO8601 like 2026-01-05T09:00:00Z or with offset
      - MT5-like 'YYYY.MM.DD HH:MI' or 'YYYY.MM.DD HH:MI:SS'
    Returns epoch seconds UTC (best effort).
    """
    s = s.strip()
    if s.endswith("Z"):
        dt = datetime.fromisoformat(s[:-1]).replace(tzinfo=timezone.utc)
        return int(dt.timestamp())
    if "T" in s:
        dt = datetime.fromisoformat(s)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return int(dt.timestamp())
    if len(s) == 16:
        dt = datetime.strptime(s, "%Y.%m.%d %H:%M").replace(tzinfo=timezone.utc)
    else:
        dt = datetime.strptime(s, "%Y.%m.%d %H:%M:%S").replace(tzinfo=timezone.utc)
    return int(dt.timestamp())


def read_ohlc_csv(path: str) -> Tuple[List[int], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    times: List[int] = []
    o, h, l, c = [], [], [], []
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            times.append(_parse_time(row["time"]))
            o.append(float(row["open"]))
            h.append(float(row["high"]))
            l.append(float(row["low"]))
            c.append(float(row["close"]))

    idx = np.argsort(np.array(times, dtype=np.int64))
    times = [times[i] for i in idx]
    o = np.array([o[i] for i in idx], dtype=float)
    h = np.array([h[i] for i in idx], dtype=float)
    l = np.array([l[i] for i in idx], dtype=float)
    c = np.array([c[i] for i in idx], dtype=float)
    return times, o, h, l, c


def epoch_to_iso(t: int) -> str:
    if t == 0:
        return ""
    return datetime.fromtimestamp(t, tz=timezone.utc).isoformat().replace("+00:00", "Z")


def state_str(x: int) -> str:
    return "PASS" if x == int(GateState.PASS) else "FAIL"


def band_str(x: int) -> str:
    if x == int(GateBand.POOR): return "POOR"
    if x == int(GateBand.OK): return "OK"
    if x == int(GateBand.GOOD): return "GOOD"
    return "GREAT"


def _cfg_from_profile(doc: Dict[str, Any], profile: str, side: str) -> QGateV12Config:
    qg = doc["qgate"]
    defaults = qg["defaults"]
    prof = qg["profiles"][profile][side]

    atr = defaults["atr"]
    er = defaults["er"]
    tr = defaults["tr"]
    adx = defaults["adx"]
    bands = defaults["bands"]

    cls = prof["classification"]
    weights = prof["weights"]

    return QGateV12Config(
        atr_period=int(atr["period"]),
        atr_baseline_len=int(atr["baseline_len"]),
        atr_ratio_ref=float(atr["ratio_ref"]),
        er_lookback=int(er["lookback"]),
        er_ref=float(er["ref"]),
        tratr_ref=float(tr["ref_tr_atr"]),
        adx=ADXCfg(enabled=bool(adx["enabled"]), period=int(adx["period"]), ref=float(adx["ref"])),
        w_atr=float(weights["atr"]),
        w_er=float(weights["er"]),
        w_tr=float(weights["tr"]),
        w_adx=float(weights["adx"]),
        classification=ClassifyCfg(
            hysteresis=bool(cls["hysteresis"]),
            q_pass=float(cls["q_pass"]),
            q_fail=float(cls["q_fail"]),
            bands=BandsCfg(float(bands["poor"]), float(bands["ok"]), float(bands["good"]))
        ),
        vetoes_enabled=False
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True, choices=["m5_m30","m15_h1"])
    ap.add_argument("--ltf", required=True, help="LTF OHLC CSV (time,open,high,low,close)")
    ap.add_argument("--htf", required=True, help="HTF OHLC CSV (time,open,high,low,close)")
    ap.add_argument("--out", required=True, help="Output golden CSV")
    ap.add_argument("--config", default="../config/qgate_v12_profiles.yaml", help="Profiles YAML")
    args = ap.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        doc = yaml.safe_load(f)

    cfg_ltf = _cfg_from_profile(doc, args.profile, "ltf")
    cfg_htf = _cfg_from_profile(doc, args.profile, "htf")

    ltf_times, lo, lh, ll, lc = read_ohlc_csv(args.ltf)
    htf_times, ho, hh, hl, hc = read_ohlc_csv(args.htf)

    ltf = compute_qgate_v12(lo, lh, ll, lc, cfg_ltf)
    htf = compute_qgate_v12(ho, hh, hl, hc, cfg_htf)

    mapped_htf_time, mapped_q_htf = align_htf_to_ltf(ltf_times, htf_times, htf.q)
    _, mapped_state_htf = align_htf_to_ltf(ltf_times, htf_times, htf.state.astype(float))
    _, mapped_band_htf  = align_htf_to_ltf(ltf_times, htf_times, htf.band.astype(float))
    _, mapped_allow_htf = align_htf_to_ltf(ltf_times, htf_times, htf.allow.astype(float))

    mapped_state_htf = mapped_state_htf.astype(int)
    mapped_band_htf  = mapped_band_htf.astype(int)
    mapped_allow_htf = mapped_allow_htf.astype(float) >= 0.5

    allow_combined = ltf.allow & mapped_allow_htf

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "time_ltf",
            "q_ltf","state_ltf","band_ltf","allow_ltf",
            "time_htf_mapped",
            "q_htf_mapped","state_htf_mapped","band_htf_mapped","allow_htf_mapped",
            "allow_combined"
        ])
        for i in range(len(ltf_times)):
            w.writerow([
                epoch_to_iso(ltf_times[i]),
                f"{ltf.q[i]:.10f}", state_str(int(ltf.state[i])), band_str(int(ltf.band[i])), int(ltf.allow[i]),
                epoch_to_iso(mapped_htf_time[i]),
                f"{mapped_q_htf[i]:.10f}", state_str(int(mapped_state_htf[i])), band_str(int(mapped_band_htf[i])), int(mapped_allow_htf[i]),
                int(allow_combined[i])
            ])

    print(f"Wrote {args.out} with {len(ltf_times)} rows for profile={args.profile}.")


if __name__ == "__main__":
    main()
