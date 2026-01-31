from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Optional, List, Tuple
import numpy as np


class GateState(IntEnum):
    FAIL = 0
    PASS = 1


class GateBand(IntEnum):
    POOR = 0
    OK = 1
    GOOD = 2
    GREAT = 3


def clamp01(x: np.ndarray) -> np.ndarray:
    return np.clip(x, 0.0, 1.0)


def safe_div(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    out = np.zeros_like(num, dtype=float)
    mask = den != 0.0
    out[mask] = num[mask] / den[mask]
    return out


def true_range(high: np.ndarray, low: np.ndarray, close: np.ndarray) -> np.ndarray:
    """
    TR[t] = max(high-low, abs(high-prev_close), abs(low-prev_close))
    For t=0, prev_close is close[0] (conventional).
    """
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    close = np.asarray(close, dtype=float)

    prev_close = np.empty_like(close)
    prev_close[0] = close[0]
    prev_close[1:] = close[:-1]

    a = high - low
    b = np.abs(high - prev_close)
    c = np.abs(low - prev_close)
    return np.maximum(a, np.maximum(b, c))


def wilder_atr(high: np.ndarray, low: np.ndarray, close: np.ndarray, period: int) -> np.ndarray:
    """
    Wilder ATR (RMA of TR):
      ATR[period-1] = mean(TR[0..period-1])
      ATR[t] = (ATR[t-1]*(period-1) + TR[t]) / period
    """
    tr = true_range(high, low, close)
    n = len(tr)
    atr = np.full(n, np.nan, dtype=float)
    if period <= 0 or n < period:
        return atr

    atr[period - 1] = np.mean(tr[:period])
    alpha_num = period - 1.0
    for t in range(period, n):
        atr[t] = (atr[t - 1] * alpha_num + tr[t]) / period
    return atr


def sma_strict(x: np.ndarray, length: int) -> np.ndarray:
    """
    Strict SMA: only produces a value when all window values are finite.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    out = np.full(n, np.nan, dtype=float)
    if length <= 0 or n < length:
        return out
    for i in range(length - 1, n):
        w = x[i - length + 1 : i + 1]
        if np.any(~np.isfinite(w)):
            continue
        out[i] = float(np.mean(w))
    return out


def efficiency_ratio(close: np.ndarray, lookback: int) -> np.ndarray:
    """
    Kaufman ER:
      net = abs(close[t] - close[t-L])
      den = sum_{k=t-L+1..t} abs(close[k] - close[k-1])
      ER = net/den (if den==0 -> 0)
    """
    close = np.asarray(close, dtype=float)
    n = len(close)
    out = np.full(n, np.nan, dtype=float)
    if lookback <= 0 or n <= lookback:
        return out

    diffs = np.abs(np.diff(close))
    for t in range(lookback, n):
        net = abs(close[t] - close[t - lookback])
        den = float(np.sum(diffs[t - lookback : t]))
        out[t] = (net / den) if den != 0.0 else 0.0
    return out


def adx_wilder(high: np.ndarray, low: np.ndarray, close: np.ndarray, period: int) -> np.ndarray:
    """
    Wilder ADX with strict seeding for parity.
    Returns NaN until first ADX seed at index (2*period - 2).
    """
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    close = np.asarray(close, dtype=float)

    n = len(close)
    out = np.full(n, np.nan, dtype=float)
    if period <= 0 or n < (2*period - 1):
        return out

    tr = true_range(high, low, close)

    dm_p = np.zeros(n, dtype=float)
    dm_m = np.zeros(n, dtype=float)
    for t in range(1, n):
        up = high[t] - high[t-1]
        down = low[t-1] - low[t]
        if up > down and up > 0:
            dm_p[t] = up
        if down > up and down > 0:
            dm_m[t] = down

    tr_sm = np.full(n, np.nan, dtype=float)
    dmp_sm = np.full(n, np.nan, dtype=float)
    dmm_sm = np.full(n, np.nan, dtype=float)

    seed_idx = period - 1
    tr_sm[seed_idx] = np.sum(tr[:period])
    dmp_sm[seed_idx] = np.sum(dm_p[:period])
    dmm_sm[seed_idx] = np.sum(dm_m[:period])

    for t in range(seed_idx + 1, n):
        tr_sm[t] = tr_sm[t-1] - (tr_sm[t-1] / period) + tr[t]
        dmp_sm[t] = dmp_sm[t-1] - (dmp_sm[t-1] / period) + dm_p[t]
        dmm_sm[t] = dmm_sm[t-1] - (dmm_sm[t-1] / period) + dm_m[t]

    di_p = np.zeros(n, dtype=float)
    di_m = np.zeros(n, dtype=float)
    valid = np.isfinite(tr_sm) & (tr_sm != 0.0)
    di_p[valid] = 100.0 * (dmp_sm[valid] / tr_sm[valid])
    di_m[valid] = 100.0 * (dmm_sm[valid] / tr_sm[valid])

    dx = np.zeros(n, dtype=float)
    denom = di_p + di_m
    m = denom != 0.0
    dx[m] = 100.0 * (np.abs(di_p[m] - di_m[m]) / denom[m])

    adx_seed_idx = 2*period - 2
    out[adx_seed_idx] = np.mean(dx[seed_idx:adx_seed_idx+1])

    for t in range(adx_seed_idx + 1, n):
        out[t] = (out[t-1] * (period - 1.0) + dx[t]) / period

    return out


@dataclass(frozen=True)
class BandsCfg:
    poor: float = 0.55
    ok: float = 0.70
    good: float = 0.85


@dataclass(frozen=True)
class ClassifyCfg:
    hysteresis: bool = True
    q_pass: float = 0.70
    q_fail: float = 0.65
    bands: BandsCfg = BandsCfg()


@dataclass(frozen=True)
class ADXCfg:
    enabled: bool = True
    period: int = 14
    ref: float = 25.0


@dataclass(frozen=True)
class QGateV12Config:
    atr_period: int = 14
    atr_baseline_len: int = 50
    atr_ratio_ref: float = 1.25

    er_lookback: int = 20
    er_ref: float = 0.35

    tratr_ref: float = 1.00

    adx: ADXCfg = ADXCfg()

    w_atr: float = 0.28
    w_er: float = 0.28
    w_tr: float = 0.14
    w_adx: float = 0.30

    classification: ClassifyCfg = ClassifyCfg()

    vetoes_enabled: bool = False
    atr_min: float = 0.0
    spread_max_points: float = 0.0


@dataclass(frozen=True)
class QGateV12Arrays:
    q: np.ndarray
    allow: np.ndarray
    state: np.ndarray
    band: np.ndarray

    s_atr: np.ndarray
    s_er: np.ndarray
    s_tr: np.ndarray
    s_adx: np.ndarray

    atr: np.ndarray
    atr_sma: np.ndarray
    atr_ratio: np.ndarray
    er: np.ndarray
    tr: np.ndarray
    tr_atr: np.ndarray
    adx: np.ndarray


def _band_from_q(q: float, bands: BandsCfg) -> GateBand:
    if q < bands.poor:
        return GateBand.POOR
    if q < bands.ok:
        return GateBand.OK
    if q < bands.good:
        return GateBand.GOOD
    return GateBand.GREAT


def compute_qgate_v12(
    open_: np.ndarray,
    high: np.ndarray,
    low: np.ndarray,
    close: np.ndarray,
    cfg: QGateV12Config,
    spread_points: Optional[np.ndarray] = None
) -> QGateV12Arrays:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    close = np.asarray(close, dtype=float)

    tr = true_range(high, low, close)
    atr = wilder_atr(high, low, close, cfg.atr_period)
    atr_sma = sma_strict(atr, cfg.atr_baseline_len)
    atr_ratio = safe_div(atr, atr_sma)

    er = efficiency_ratio(close, cfg.er_lookback)
    tr_atr = safe_div(tr, atr)

    adx = adx_wilder(high, low, close, cfg.adx.period) if cfg.adx.enabled else np.full_like(close, np.nan, dtype=float)

    s_atr = clamp01(safe_div(atr_ratio, np.array(cfg.atr_ratio_ref)))
    s_er  = clamp01(safe_div(er,       np.array(cfg.er_ref)))
    s_tr  = clamp01(safe_div(tr_atr,   np.array(cfg.tratr_ref)))
    s_adx = clamp01(safe_div(adx,      np.array(cfg.adx.ref))) if cfg.adx.enabled else np.zeros_like(close, dtype=float)

    invalid_base = (~np.isfinite(atr)) | (~np.isfinite(atr_sma)) | (~np.isfinite(er)) | (atr <= 0.0) | (atr_sma <= 0.0)
    invalid_adx = cfg.adx.enabled & (~np.isfinite(adx))
    invalid = invalid_base | invalid_adx

    w_atr = cfg.w_atr
    w_er = cfg.w_er
    w_tr = cfg.w_tr
    w_adx = cfg.w_adx if cfg.adx.enabled else 0.0
    wsum = w_atr + w_er + w_tr + w_adx

    q = np.zeros_like(close, dtype=float)
    if wsum > 0.0:
        q = (w_atr*s_atr + w_er*s_er + w_tr*s_tr + w_adx*s_adx) / wsum
    q = np.where(invalid, 0.0, q)

    n = len(q)
    band = np.zeros(n, dtype=np.int32)
    state = np.zeros(n, dtype=np.int32)
    allow = np.zeros(n, dtype=bool)

    cls = cfg.classification
    q_pass = float(cls.q_pass)
    q_fail = float(cls.q_fail) if cls.hysteresis else float(cls.q_pass)

    prev_state = GateState.FAIL
    for i in range(n):
        qi = float(q[i])
        band[i] = int(_band_from_q(qi, cls.bands))

        if invalid[i]:
            st = GateState.FAIL
        else:
            if not cls.hysteresis:
                st = GateState.PASS if qi >= q_pass else GateState.FAIL
            else:
                if prev_state == GateState.FAIL:
                    st = GateState.PASS if qi >= q_pass else GateState.FAIL
                else:
                    st = GateState.FAIL if qi <= q_fail else GateState.PASS

        state[i] = int(st)
        prev_state = st

        ok = (st == GateState.PASS)

        if cfg.vetoes_enabled:
            if cfg.atr_min > 0.0 and float(atr[i]) < cfg.atr_min:
                ok = False
            if cfg.spread_max_points > 0.0:
                if spread_points is None:
                    raise ValueError("spread_points required when spread veto enabled")
                if float(spread_points[i]) > cfg.spread_max_points:
                    ok = False

        allow[i] = ok

    return QGateV12Arrays(
        q=q, allow=allow, state=state, band=band,
        s_atr=s_atr, s_er=s_er, s_tr=s_tr, s_adx=s_adx,
        atr=atr, atr_sma=atr_sma, atr_ratio=atr_ratio,
        er=er, tr=tr, tr_atr=tr_atr, adx=adx
    )


def align_htf_to_ltf(
    ltf_times: List[int],
    htf_times: List[int],
    htf_values: np.ndarray
) -> Tuple[List[int], np.ndarray]:
    """
    No-lookahead asof alignment:
      for each ltf_time, pick the last htf_time <= ltf_time.
    times are epoch seconds (int). Must be ascending.
    Returns mapped_htf_time per ltf bar and mapped values array.
    """
    mapped_t: List[int] = []
    mapped_v = np.zeros(len(ltf_times), dtype=float)

    j = 0
    last_idx = -1
    for i, t in enumerate(ltf_times):
        while j < len(htf_times) and htf_times[j] <= t:
            last_idx = j
            j += 1
        if last_idx < 0:
            mapped_t.append(0)
            mapped_v[i] = 0.0
        else:
            mapped_t.append(htf_times[last_idx])
            mapped_v[i] = float(htf_values[last_idx])
    return mapped_t, mapped_v
