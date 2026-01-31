#pragma once

enum GateState { GS_FAIL=0, GS_PASS=1 };
enum GateBand  { GB_POOR=0, GB_OK=1, GB_GOOD=2, GB_GREAT=3 };

struct BandsCfg {
   double poor;
   double ok;
   double good;
};

struct ClassifyCfg {
   bool   hysteresis;
   double q_pass;
   double q_fail;
   BandsCfg bands;
};

struct QGateV11Cfg {
   int    atrPeriod;
   int    atrBaselineLen;
   double atrRatioRef;

   int    erLookback;
   double erRef;

   double tratrRef;

   double wAtr, wEr, wTr;

   ClassifyCfg cls;

   bool   vetoesEnabled;
   double atrMin;
   double spreadMaxPoints;
};

struct QGateV11Out {
   double q;
   int    state; // GateState
   int    band;  // GateBand
   bool   allow;

   // debug
   double atr, atrSma, atrRatio;
   double er, tr, trAtr;
   double s_atr, s_er, s_tr;
};

static double _Clamp01(const double x) {
   if(x < 0.0) return 0.0;
   if(x > 1.0) return 1.0;
   return x;
}

static int _BandFromQ(const double q, const BandsCfg &b) {
   if(q < b.poor) return (int)GB_POOR;
   if(q < b.ok)   return (int)GB_OK;
   if(q < b.good) return (int)GB_GOOD;
   return (int)GB_GREAT;
}

static double _TrueRange(const double high, const double low, const double prevClose) {
   const double a = high - low;
   const double b = MathAbs(high - prevClose);
   const double c = MathAbs(low  - prevClose);
   return MathMax(a, MathMax(b, c));
}

// Arrays must be ascending by time: 0 oldest .. n-1 newest.
static void ComputeTR(const double &high[], const double &low[], const double &close[], double &tr[]) {
   const int n = ArraySize(close);
   ArrayResize(tr, n);
   for(int i=0;i<n;i++){
      double prevClose = (i==0 ? close[0] : close[i-1]);
      tr[i] = _TrueRange(high[i], low[i], prevClose);
   }
}

static void ComputeWilderATR(const double &tr[], const int period, double &atr[]) {
   const int n = ArraySize(tr);
   ArrayResize(atr, n);
   for(int i=0;i<n;i++) atr[i] = EMPTY_VALUE;

   if(period <= 0 || n < period) return;

   double sum = 0.0;
   for(int i=0;i<period;i++) sum += tr[i];
   atr[period-1] = sum / (double)period;

   const double alpha_num = (double)(period - 1);
   for(int i=period;i<n;i++){
      atr[i] = (atr[i-1]*alpha_num + tr[i]) / (double)period;
   }
}

static void ComputeSMA_Strict(const double &x[], const int len, double &out[]) {
   const int n = ArraySize(x);
   ArrayResize(out, n);
   for(int i=0;i<n;i++) out[i] = EMPTY_VALUE;

   if(len <= 0 || n < len) return;

   for(int i=len-1;i<n;i++){
      bool ok = true;
      double sum = 0.0;
      for(int j=i-len+1; j<=i; j++){
         if(x[j] == EMPTY_VALUE) { ok=false; break; }
         sum += x[j];
      }
      if(ok) out[i] = sum / (double)len;
   }
}

static void ComputeER(const double &close[], const int lookback, double &er[]) {
   const int n = ArraySize(close);
   ArrayResize(er, n);
   for(int i=0;i<n;i++) er[i] = EMPTY_VALUE;

   if(lookback <= 0 || n <= lookback) return;

   for(int t=lookback; t<n; t++){
      const double net = MathAbs(close[t] - close[t-lookback]);
      double den = 0.0;
      for(int k=t-lookback+1; k<=t; k++){
         den += MathAbs(close[k] - close[k-1]);
      }
      er[t] = (den != 0.0 ? net/den : 0.0);
   }
}

static void ComputeQGateV11_All(
   const datetime &time[],
   const double &open[],
   const double &high[],
   const double &low[],
   const double &close[],
   const QGateV11Cfg &cfg,
   QGateV11Out &out[]
){
   const int n = ArraySize(close);
   ArrayResize(out, n);

   double tr[], atr[], atrSma[], er[];
   ComputeTR(high, low, close, tr);
   ComputeWilderATR(tr, cfg.atrPeriod, atr);
   ComputeSMA_Strict(atr, cfg.atrBaselineLen, atrSma);
   ComputeER(close, cfg.erLookback, er);

   int prevState = (int)GS_FAIL;
   const double wsum = cfg.wAtr + cfg.wEr + cfg.wTr;
   const bool useHys = cfg.cls.hysteresis;

   for(int i=0;i<n;i++){
      QGateV11Out o;
      o.q = 0.0;
      o.state = (int)GS_FAIL;
      o.band = (int)GB_POOR;
      o.allow = false;

      o.tr = tr[i];
      o.atr = atr[i];
      o.atrSma = atrSma[i];
      o.er = er[i];

      bool invalid = (o.atr == EMPTY_VALUE) || (o.atrSma == EMPTY_VALUE) || (o.er == EMPTY_VALUE) ||
                     (o.atr <= 0.0) || (o.atrSma <= 0.0) || (wsum <= 0.0);

      o.atrRatio = invalid ? 0.0 : (o.atr / o.atrSma);
      o.trAtr    = (o.atr > 0.0 ? o.tr / o.atr : 0.0);

      o.s_atr = invalid ? 0.0 : _Clamp01(o.atrRatio / cfg.atrRatioRef);
      o.s_er  = invalid ? 0.0 : _Clamp01(o.er / cfg.erRef);
      o.s_tr  = invalid ? 0.0 : _Clamp01(o.trAtr / cfg.tratrRef);

      o.q = invalid ? 0.0 : (cfg.wAtr*o.s_atr + cfg.wEr*o.s_er + cfg.wTr*o.s_tr) / wsum;
      o.band = _BandFromQ(o.q, cfg.cls.bands);

      if(invalid){
         o.state = (int)GS_FAIL;
      } else {
         if(!useHys){
            o.state = (o.q >= cfg.cls.q_pass ? (int)GS_PASS : (int)GS_FAIL);
         } else {
            if(prevState == (int)GS_FAIL){
               o.state = (o.q >= cfg.cls.q_pass ? (int)GS_PASS : (int)GS_FAIL);
            } else {
               o.state = (o.q <= cfg.cls.q_fail ? (int)GS_FAIL : (int)GS_PASS);
            }
         }
      }

      prevState = o.state;

      bool allow = (o.state == (int)GS_PASS);

      if(cfg.vetoesEnabled){
         if(cfg.atrMin > 0.0 && o.atr < cfg.atrMin) allow = false;
         // spread veto requires external input; excluded in golden tests
      }

      o.allow = allow;
      out[i] = o;
   }
}
