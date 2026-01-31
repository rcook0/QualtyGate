from __future__ import annotations

import csv
from datetime import datetime, timezone
from typing import List, Tuple

import numpy as np

from qgate_v11 import (
    QGateV11Config, ClassifyCfg, BandsCfg,
    compute_qgate_v11, align_htf_to_ltf, GateState, GateBand
)

def _parse_time(s: str) -> int:
    """
    Accepts:
      - ISO8601 like 2026-01-05T09:00:00Z or with offset
      - MT5-like 'YYYY.MM.DD HH:MI' or 'YYYY.MM.DD HH:MI:SS'
    Returns epoch seconds UTC (best effort).
    """
    s = s.strip()
    # ISO 'Z'
    if s.endswith("Z"):
        s2 = s[:-1]
        dt = datetime.fromisoformat(s2)
        dt = dt.replace(tzinfo=timezone.utc)
        return int(dt.timestamp())

    # ISO with T
    if "T" in s:
        dt = datetime.fromisoformat(s)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return int(dt.timestamp())

    # MT5 style: YYYY.MM.DD HH:MI[:SS]
    # Treat as UTC if tz absent; for parity, keep both files in the same time basis.
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
            t = _parse_time(row["time"])
            times.append(t)
            o.append(float(row["open"]))
            h.append(float(row["high"]))
            l.append(float(row["low"]))
            c.append(float(row["close"]))

    # sort ascending by time
    idx = np.argsort(np.array(times, dtype=np.int64))
    times = [times[i] for i in idx]
    o = np.array([o[i] for i in idx], dtype=float)
    h = np.array([h[i] for i in idx], dtype=float)
    l = np.array([l[i] for i in idx], dtype=float)
    c = np.array([c[i] for i in idx], dtype=float)
    return times, o, h, l, c


def state_str(x: int) -> str:
    return "PASS" if x == int(GateState.PASS) else "FAIL"


def band_str(x: int) -> str:
    if x == int(GateBand.POOR): return "POOR"
    if x == int(GateBand.OK): return "OK"
    if x == int(GateBand.GOOD): return "GOOD"
    return "GREAT"


def epoch_to_iso(t: int) -> str:
    if t == 0:
        return ""
    return datetime.fromtimestamp(t, tz=timezone.utc).isoformat().replace("+00:00", "Z")


def main():
    # Place ltf.csv/htf.csv next to this script, or edit paths.
    ltf_path = "ltf.csv"
    htf_path = "htf.csv"
    out_path = "golden_out.csv"

    # Choose one profile and keep it identical in MT5 EA.
    # Profile B default: M15 -> H1
    cfg_ltf = QGateV11Config(
        classification=ClassifyCfg(hysteresis=True, q_pass=0.72, q_fail=0.66, bands=BandsCfg(0.55, 0.70, 0.85)),
        vetoes_enabled=False
    )
    cfg_htf = QGateV11Config(
        classification=ClassifyCfg(hysteresis=True, q_pass=0.75, q_fail=0.70, bands=BandsCfg(0.55, 0.70, 0.85)),
        vetoes_enabled=False
    )

    ltf_times, lo, lh, ll, lc = read_ohlc_csv(ltf_path)
    htf_times, ho, hh, hl, hc = read_ohlc_csv(htf_path)

    ltf = compute_qgate_v11(lo, lh, ll, lc, cfg_ltf)
    htf = compute_qgate_v11(ho, hh, hl, hc, cfg_htf)

    mapped_htf_time, mapped_q_htf = align_htf_to_ltf(ltf_times, htf_times, htf.q)
    _, mapped_state_htf = align_htf_to_ltf(ltf_times, htf_times, htf.state.astype(float))
    _, mapped_band_htf  = align_htf_to_ltf(ltf_times, htf_times, htf.band.astype(float))
    _, mapped_allow_htf = align_htf_to_ltf(ltf_times, htf_times, htf.allow.astype(float))

    mapped_state_htf = mapped_state_htf.astype(int)
    mapped_band_htf  = mapped_band_htf.astype(int)
    mapped_allow_htf = mapped_allow_htf.astype(float) >= 0.5

    allow_combined = ltf.allow & mapped_allow_htf

    with open(out_path, "w", newline="") as f:
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

    print(f"Wrote {out_path} with {len(ltf_times)} rows.")


if __name__ == "__main__":
    main()
