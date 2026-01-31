# QGate v1.2 Specification (Delta from v1.1)

## 1. Purpose
QGate evaluates market regime quality to decide whether **new risk may be added**:
- It does not generate signals
- It does not manage trades
- It provides:
  - `q ∈ [0,1]`
  - `state ∈ {FAIL,PASS}` with optional hysteresis
  - `band ∈ {POOR,OK,GOOD,GREAT}` for interpretability
  - optional HTF overlay alignment

v1.2 extends v1.1 with an additional regime component: **ADX**.

---

## 2. Inputs
- Required: OHLC aligned bars
- Optional: spread series (for vetoes; not part of q)
- Config:
  - ATR: period, baseline length, ratio_ref
  - ER: lookback, ref
  - TR/ATR: ref
  - ADX: period, ref, enabled
  - weights per component
  - classification thresholds & hysteresis
  - HTF profile (optional)

---

## 3. Core metrics (v1.2 addition: ADX)

### 3.1 Wilder Directional Movement (DM) and ADX
Let:
- `upMove = high[t] - high[t-1]`
- `downMove = low[t-1] - low[t]`

Define:
- `DM+ = upMove` if `upMove > downMove` and `upMove > 0` else `0`
- `DM- = downMove` if `downMove > upMove` and `downMove > 0` else `0`

Compute True Range `TR[t]` as in v1.0.

Wilder smoothing over period `N`:
- `TR_sm[t]  = TR_sm[t-1]  - TR_sm[t-1]/N  + TR[t]`
- `DM+_sm[t] = DM+_sm[t-1] - DM+_sm[t-1]/N + DM+[t]`
- `DM-_sm[t] = DM-_sm[t-1] - DM-_sm[t-1]/N + DM-[t]`

Seed at `t = N-1`:
- `TR_sm[N-1]  = sum(TR[0..N-1])`
- `DM+_sm[N-1] = sum(DM+[0..N-1])`
- `DM-_sm[N-1] = sum(DM-[0..N-1])`

Directional Indicators:
- `DI+ = 100 * DM+_sm / TR_sm`
- `DI- = 100 * DM-_sm / TR_sm`

Directional Index:
- `DX = 100 * abs(DI+ - DI-) / (DI+ + DI-)` (if denom=0 -> 0)

ADX seeding:
- First ADX at `t = 2N - 2`:
  - `ADX[2N-2] = mean(DX[N-1 .. 2N-2])`

Then Wilder smoothing:
- `ADX[t] = (ADX[t-1]*(N-1) + DX[t]) / N`

---

## 4. Normalization and scoring
Existing components:
- `s_atr = clamp( (ATR/ATR_SMA) / atr_ratio_ref, 0, 1 )`
- `s_er  = clamp( ER / er_ref, 0, 1 )`
- `s_tr  = clamp( (TR/ATR) / tratr_ref, 0, 1 )`

New component:
- `s_adx = clamp( ADX / adx_ref, 0, 1 )`

Composite:
- `q = weighted_average( enabled component scores )`

v1.2 default `adx_ref = 25` (typical “trend strength” threshold).

---

## 5. Classification
Same as v1.1:
- bands by q cutpoints (default: 0.55 / 0.70 / 0.85)
- PASS/FAIL based on `q_pass` and optional hysteresis (`q_fail`)

---

## 6. HTF overlay
Same as v1.1:
- align HTF to LTF using **last closed HTF bar** (no-lookahead)
- default combine: AND (`allow_ltf && allow_htf`)

---

## 7. Backward compatibility
- If `adx.enabled=false`, v1.2 reduces to v1.1 behavior (aside from config defaults).
