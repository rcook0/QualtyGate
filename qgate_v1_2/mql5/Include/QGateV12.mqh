#pragma once

enum GateState { GS_FAIL=0, GS_PASS=1 };
enum GateBand  { GB_POOR=0, GB_OK=1, GB_GOOD=2, GB_GREAT=3 };

struct BandsCfg { double poor; double ok; double good; };

struct ClassifyCfg {
   bool   hysteresis;
   double q_pass;
   double q_fail;
   BandsCfg bands;
};

struct ADXCfg { bool enabled; int period; double ref; };

struct QGateV12Cfg {
   int    atrPeriod;
   int    atrBaselineLen;
   double atrRatioRef;

   int    erLookback;
   double erRef;

   double tratrRef;

   ADXCfg adx;

   double wAtr, wEr, wTr, wAdx;

   ClassifyCfg cls;

   bool   vetoesEnabled;
   double atrMin;
   double spreadMaxPoints;
};

struct QGateV12Out {
   double q;
   int    state;
   int    band;
   bool   allow;

   double s_atr, s_er, s_tr, s_adx;

   double atr, atrSma, atrRatio;
   double er, tr, trAtr;
   double adx;
};

static double _Clamp01(const double x) { if(x < 0.0) return 0.0; if(x > 1.0) return 1.0; return x; }

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

static void ComputeTR(const double &high[], const double &low[], const double &close[], double &tr[]) {
   const int n = ArraySize(close);
   ArrayResize(tr, n);
   for(int i=0;i<n;i++){
      double prevClose = (i==0 ? close[0] : close[i-1]);
      tr[i] = _TrueRange(high[i], low[i], prevClose);
   }
}

static void ComputeWilderATR_fromTR(const double &tr[], const int period, double &atr[]) {
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

static void ComputeADX_Wilder(const double &high[], const double &low[], const double &close[], const int period, double &adx[]) {
   const int n = ArraySize(close);
   ArrayResize(adx, n);
   for(int i=0;i<n;i++) adx[i] = EMPTY_VALUE;
   if(period <= 0 || n < (2*period - 1)) return;

   double tr[];
   ComputeTR(high, low, close, tr);

   double dmp[], dmm[];
   ArrayResize(dmp, n);
   ArrayResize(dmm, n);
   for(int i=0;i<n;i++){ dmp[i]=0.0; dmm[i]=0.0; }

   for(int t=1; t<n; t++){
      double up = high[t] - high[t-1];
      double down = low[t-1] - low[t];
      if(up > down && up > 0.0) dmp[t] = up;
      if(down > up && down > 0.0) dmm[t] = down;
   }

   const int seed = period - 1;

   double tr_sm[], dmp_sm[], dmm_sm[];
   ArrayResize(tr_sm, n);
   ArrayResize(dmp_sm, n);
   ArrayResize(dmm_sm, n);
   for(int i=0;i<n;i++){ tr_sm[i]=EMPTY_VALUE; dmp_sm[i]=EMPTY_VALUE; dmm_sm[i]=EMPTY_VALUE; }

   double sumTR=0.0, sumP=0.0, sumM=0.0;
   for(int i=0;i<period;i++){ sumTR += tr[i]; sumP += dmp[i]; sumM += dmm[i]; }
   tr_sm[seed]=sumTR; dmp_sm[seed]=sumP; dmm_sm[seed]=sumM;

   for(int t=seed+1; t<n; t++){
      tr_sm[t]  = tr_sm[t-1]  - (tr_sm[t-1]  / period) + tr[t];
      dmp_sm[t] = dmp_sm[t-1] - (dmp_sm[t-1] / period) + dmp[t];
      dmm_sm[t] = dmm_sm[t-1] - (dmm_sm[t-1] / period) + dmm[t];
   }

   double dip[], dim[];
   ArrayResize(dip, n);
   ArrayResize(dim, n);
   for(int i=0;i<n;i++){ dip[i]=0.0; dim[i]=0.0; }

   for(int t=seed; t<n; t++){
      if(tr_sm[t] == 0.0 || tr_sm[t] == EMPTY_VALUE) continue;
      dip[t] = 100.0 * (dmp_sm[t] / tr_sm[t]);
      dim[t] = 100.0 * (dmm_sm[t] / tr_sm[t]);
   }

   double dx[];
   ArrayResize(dx, n);
   for(int i=0;i<n;i++) dx[i]=0.0;

   for(int t=seed; t<n; t++){
      double denom = dip[t] + dim[t];
      if(denom == 0.0) { dx[t]=0.0; continue; }
      dx[t] = 100.0 * (MathAbs(dip[t] - dim[t]) / denom);
   }

   const int adx_seed = 2*period - 2;
   double sumDX=0.0;
   for(int t=seed; t<=adx_seed; t++) sumDX += dx[t];
   adx[adx_seed] = sumDX / (double)period;

   for(int t=adx_seed+1; t<n; t++){
      adx[t] = (adx[t-1] * (period - 1.0) + dx[t]) / (double)period;
   }
}

static void ComputeQGateV12_All(
   const datetime &time[],
   const double &open[],
   const double &high[],
   const double &low[],
   const double &close[],
   const QGateV12Cfg &cfg,
   QGateV12Out &out[]
){
   const int n = ArraySize(close);
   ArrayResize(out, n);

   double tr[], atr[], atrSma[], er[];
   ComputeTR(high, low, close, tr);
   ComputeWilderATR_fromTR(tr, cfg.atrPeriod, atr);
   ComputeSMA_Strict(atr, cfg.atrBaselineLen, atrSma);
   ComputeER(close, cfg.erLookback, er);

   double adx[];
   if(cfg.adx.enabled) ComputeADX_Wilder(high, low, close, cfg.adx.period, adx);
   else { ArrayResize(adx, n); for(int i=0;i<n;i++) adx[i]=EMPTY_VALUE; }

   int prevState = (int)GS_FAIL;
   const bool useHys = cfg.cls.hysteresis;

   const double wAdx = (cfg.adx.enabled ? cfg.wAdx : 0.0);
   const double wsum = cfg.wAtr + cfg.wEr + cfg.wTr + wAdx;

   for(int i=0;i<n;i++){
      QGateV12Out o;
      o.q=0.0; o.state=(int)GS_FAIL; o.band=(int)GB_POOR; o.allow=false;

      o.tr=tr[i]; o.atr=atr[i]; o.atrSma=atrSma[i]; o.er=er[i]; o.adx=adx[i];

      bool invalid_base = (o.atr==EMPTY_VALUE) || (o.atrSma==EMPTY_VALUE) || (o.er==EMPTY_VALUE) || (o.atr<=0.0) || (o.atrSma<=0.0);
      bool invalid_adx = cfg.adx.enabled && (o.adx==EMPTY_VALUE);
      bool invalid = invalid_base || invalid_adx || (wsum <= 0.0);

      o.atrRatio = invalid ? 0.0 : (o.atr / o.atrSma);
      o.trAtr    = (o.atr>0.0 ? o.tr / o.atr : 0.0);

      o.s_atr = invalid ? 0.0 : _Clamp01(o.atrRatio / cfg.atrRatioRef);
      o.s_er  = invalid ? 0.0 : _Clamp01(o.er / cfg.erRef);
      o.s_tr  = invalid ? 0.0 : _Clamp01(o.trAtr / cfg.tratrRef);
      o.s_adx = (cfg.adx.enabled && !invalid ? _Clamp01(o.adx / cfg.adx.ref) : 0.0);

      if(!invalid){
         o.q = (cfg.wAtr*o.s_atr + cfg.wEr*o.s_er + cfg.wTr*o.s_tr + wAdx*o.s_adx) / wsum;
      } else o.q = 0.0;

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
      }
      o.allow = allow;
      out[i]=o;
   }
}
