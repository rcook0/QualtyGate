# QGate v1.2 (Core + HTF Combiner + Golden Parity Tests)

This workspace delivers **QGate v1.2**:
- v1.1 features preserved (q-score, bands, PASS/FAIL with hysteresis, HTF combiner semantics)
- **New regime component:** **ADX (Average Directional Index)**, computed with Wilder smoothing, normalized into `s_adx ∈ [0,1]`
- Updated YAML + JSON Schema contracts
- **Golden parity workflow** (Python generates `golden_out.csv`; MT5 EA compares q/state/band + HTF alignment)

## What’s new in v1.2
- New component:
  - `adx[t]` (Wilder ADX, length `adx.period`)
  - `s_adx = clamp(adx / adx.ref, 0, 1)`
- Composite:
  - `q = weighted_average(enabled_components)`
  - Default weights updated to include ADX (see `config/qgate_v12_profiles.yaml`)

## Directory layout
- `config/` – YAML config + JSON Schema
- `docs/` – v1.2 specification
- `python/` – reference implementation + golden generator
- `mql5/` – MQL5 include modules + golden test EA

---

## 1) Export data (MT5) for XAUUSD (recommended)
Export *bar OPEN times* (MT5 default) and columns:
`time,open,high,low,close`

### Pair A (M5 → M30)
- LTF: **M5**: 10,000 bars → save as `ltf.csv`
- HTF: **M30**: 2,000 bars → save as `htf.csv`

### Pair B (M15 → H1)
- LTF: **M15**: 8,000 bars → `ltf.csv`
- HTF: **H1**: 2,000 bars → `htf.csv`

Place the two CSVs next to `python/golden_generate_v12.py` for the Python step.

---

## 2) Generate golden outputs in Python
Requirements:
- Python 3.10+
- `numpy`
- `pyyaml`

Install:
```bash
pip install -r python/requirements.txt
```

Run (example: M15→H1 profile):
```bash
cd python
python golden_generate_v12.py --profile m15_h1 --ltf ltf.csv --htf htf.csv --out golden_out.csv
```

This writes `golden_out.csv` containing:
- LTF q/state/band/allow
- mapped HTF q/state/band/allow (no-lookahead)
- combined allow (AND)

---

## 3) Run MT5 Golden Test EA
1. Copy files into MT5 data folder:
   - `mql5/Include/*.mqh` → `MQL5/Include/`
   - `mql5/Experts/QGateGoldenTestEA_v12.mq5` → `MQL5/Experts/`
2. Put CSV files into `MQL5/Files/`:
   - `ltf.csv`, `htf.csv`, `golden_out.csv`
3. Compile `QGateGoldenTestEA_v12.mq5` in MetaEditor.
4. Attach EA to any chart.
5. Review **Experts** log output summary.

Pass criteria:
- `Max |q_ltf diff| <= EPS_Q` and `Max |q_htf diff| <= EPS_Q`
- mismatch counts all zero:
  - LTF: q/state/band
  - HTF: mapped_time/q/state/band
  - combined allow

---

## 4) Notes on determinism and parity
- Both Python and MQL5 compute ATR/ER/TR/ADX **purely from CSV OHLC** (no indicator handles). This avoids platform ATR implementation drift.
- Time alignment uses bar **open time** and no-lookahead: `htf_time <= ltf_time`.

---

## 5) Next steps (production integration)
- Integrate `QGateV12` as a reusable component in Algo:
  - Gate **entries and scale-ins only**
  - Never gate exits/stops
- Keep vetoes outside q (spread/news/session filters remain separate gates).
