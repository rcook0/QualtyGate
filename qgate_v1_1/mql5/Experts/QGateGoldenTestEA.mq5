#property strict
#property description "QGate v1.1 Golden Test EA (no trading). Put CSVs in MQL5/Files."

#include <QGateV11.mqh>
#include <QGateHTFCombiner.mqh>

input string LTF_CSV        = "ltf.csv";
input string HTF_CSV        = "htf.csv";
input string GOLDEN_OUT_CSV = "golden_out.csv";

input double EPS_Q = 1e-6;

static string _NormTimeStr(string s){
   s = StringTrim(s);
   s = StringReplace(s, "T", " ");
   s = StringReplace(s, "Z", "");
   s = StringReplace(s, "-", ".");
   s = StringReplace(s, "/", ".");
   return s;
}

static datetime _ParseTime(string s){
   s = _NormTimeStr(s);
   return (datetime)StringToTime(s);
}

struct BarRow {
   datetime t;
   double o,h,l,c;
};

static bool ReadBarsCSV(const string filename, BarRow &rows[]){
   int fh = FileOpen(filename, FILE_READ|FILE_CSV|FILE_ANSI);
   if(fh == INVALID_HANDLE){
      Print("Failed to open ", filename, " err=", GetLastError());
      return false;
   }

   // header
   if(!FileIsEnding(fh)){
      FileReadString(fh); // time
      FileReadString(fh); // open
      FileReadString(fh); // high
      FileReadString(fh); // low
      FileReadString(fh); // close
   }

   BarRow tmp[];
   int n = 0;
   while(!FileIsEnding(fh)){
      string ts = FileReadString(fh);
      if(ts == "") break;
      BarRow r;
      r.t = _ParseTime(ts);
      r.o = FileReadNumber(fh);
      r.h = FileReadNumber(fh);
      r.l = FileReadNumber(fh);
      r.c = FileReadNumber(fh);
      ArrayResize(tmp, n+1);
      tmp[n] = r;
      n++;
   }
   FileClose(fh);

   // sort ascending by time
   for(int i=1;i<n;i++){
      BarRow key = tmp[i];
      int j=i-1;
      while(j>=0 && tmp[j].t > key.t){
         tmp[j+1] = tmp[j];
         j--;
      }
      tmp[j+1] = key;
   }

   ArrayResize(rows, n);
   for(int i=0;i<n;i++) rows[i] = tmp[i];
   Print(filename, ": loaded ", n, " bars");
   return (n > 0);
}

struct GoldenRow {
   datetime t_ltf;
   double q_ltf;
   string state_ltf;
   string band_ltf;
   int allow_ltf;

   datetime t_htf_map;
   double q_htf;
   string state_htf;
   string band_htf;
   int allow_htf;

   int allow_combined;
};

static bool ReadGoldenOut(const string filename, GoldenRow &rows[]){
   int fh = FileOpen(filename, FILE_READ|FILE_CSV|FILE_ANSI);
   if(fh == INVALID_HANDLE){
      Print("Failed to open ", filename, " err=", GetLastError());
      return false;
   }

   // header row: 11 columns
   for(int i=0;i<11;i++) FileReadString(fh);

   GoldenRow tmp[];
   int n=0;
   while(!FileIsEnding(fh)){
      string tltf = FileReadString(fh);
      if(tltf == "") break;

      GoldenRow r;
      r.t_ltf = _ParseTime(tltf);
      r.q_ltf = StringToDouble(FileReadString(fh));
      r.state_ltf = FileReadString(fh);
      r.band_ltf  = FileReadString(fh);
      r.allow_ltf = (int)FileReadNumber(fh);

      string thtf = FileReadString(fh);
      r.t_htf_map = (thtf=="" ? (datetime)0 : _ParseTime(thtf));
      r.q_htf = StringToDouble(FileReadString(fh));
      r.state_htf = FileReadString(fh);
      r.band_htf  = FileReadString(fh);
      r.allow_htf = (int)FileReadNumber(fh);

      r.allow_combined = (int)FileReadNumber(fh);

      ArrayResize(tmp, n+1);
      tmp[n]=r;
      n++;
   }
   FileClose(fh);

   // sort ascending by t_ltf
   for(int i=1;i<n;i++){
      GoldenRow key=tmp[i];
      int j=i-1;
      while(j>=0 && tmp[j].t_ltf > key.t_ltf){
         tmp[j+1]=tmp[j];
         j--;
      }
      tmp[j+1]=key;
   }

   ArrayResize(rows,n);
   for(int i=0;i<n;i++) rows[i]=tmp[i];
   Print(filename, ": loaded ", n, " rows");
   return (n>0);
}

static int _StateFromStr(const string s){
   if(s=="PASS") return (int)GS_PASS;
   return (int)GS_FAIL;
}
static int _BandFromStr(const string s){
   if(s=="POOR") return (int)GB_POOR;
   if(s=="OK") return (int)GB_OK;
   if(s=="GOOD") return (int)GB_GOOD;
   return (int)GB_GREAT;
}

int OnInit(){
   BarRow ltfBars[], htfBars[];
   GoldenRow gold[];

   if(!ReadBarsCSV(LTF_CSV, ltfBars)) return INIT_FAILED;
   if(!ReadBarsCSV(HTF_CSV, htfBars)) return INIT_FAILED;
   if(!ReadGoldenOut(GOLDEN_OUT_CSV, gold)) return INIT_FAILED;

   // Build arrays (ascending)
   int nL = ArraySize(ltfBars);
   datetime tL[]; double oL[], hL[], lL[], cL[];
   ArrayResize(tL,nL); ArrayResize(oL,nL); ArrayResize(hL,nL); ArrayResize(lL,nL); ArrayResize(cL,nL);
   for(int i=0;i<nL;i++){ tL[i]=ltfBars[i].t; oL[i]=ltfBars[i].o; hL[i]=ltfBars[i].h; lL[i]=ltfBars[i].l; cL[i]=ltfBars[i].c; }

   int nH = ArraySize(htfBars);
   datetime tH[]; double oH[], hH[], lH[], cH[];
   ArrayResize(tH,nH); ArrayResize(oH,nH); ArrayResize(hH,nH); ArrayResize(lH,nH); ArrayResize(cH,nH);
   for(int i=0;i<nH;i++){ tH[i]=htfBars[i].t; oH[i]=htfBars[i].o; hH[i]=htfBars[i].h; lH[i]=htfBars[i].l; cH[i]=htfBars[i].c; }

   // CONFIG: must match Python golden_generate.py
   QGateV11Cfg cfgL;
   cfgL.atrPeriod=14; cfgL.atrBaselineLen=50; cfgL.atrRatioRef=1.25;
   cfgL.erLookback=20; cfgL.erRef=0.35;
   cfgL.tratrRef=1.0;
   cfgL.wAtr=0.40; cfgL.wEr=0.40; cfgL.wTr=0.20;
   cfgL.cls.hysteresis=true; cfgL.cls.q_pass=0.72; cfgL.cls.q_fail=0.66;
   cfgL.cls.bands.poor=0.55; cfgL.cls.bands.ok=0.70; cfgL.cls.bands.good=0.85;
   cfgL.vetoesEnabled=false;

   QGateV11Cfg cfgH = cfgL;
   cfgH.cls.q_pass=0.75; cfgH.cls.q_fail=0.70;

   QGateV11Out outL[], outH[];
   ComputeQGateV11_All(tL, oL, hL, lL, cL, cfgL, outL);
   ComputeQGateV11_All(tH, oH, hH, lH, cH, cfgH, outH);

   int mism_q_ltf=0, mism_state_ltf=0, mism_band_ltf=0;
   int mism_htf_time=0, mism_q_htf=0, mism_state_htf=0, mism_band_htf=0;
   int mism_allow_comb=0;

   double max_abs_q_ltf=0.0, max_abs_q_htf=0.0;

   int iL = 0;

   for(int g=0; g<ArraySize(gold); g++){
      datetime t = gold[g].t_ltf;

      while(iL < nL && tL[iL] < t) iL++;
      if(iL >= nL || tL[iL] != t){
         Print("Golden time not found in LTF bars: ", TimeToString(t, TIME_DATE|TIME_MINUTES));
         continue;
      }

      double dq = MathAbs(outL[iL].q - gold[g].q_ltf);
      if(dq > max_abs_q_ltf) max_abs_q_ltf = dq;
      if(dq > EPS_Q) mism_q_ltf++;

      int stg = _StateFromStr(gold[g].state_ltf);
      int bng = _BandFromStr(gold[g].band_ltf);

      if(outL[iL].state != stg) mism_state_ltf++;
      if(outL[iL].band  != bng) mism_band_ltf++;

      int idxH = AlignIndex_NoLookahead(tH, t);
      datetime mappedTime = (idxH < 0 ? (datetime)0 : tH[idxH]);

      if(mappedTime != gold[g].t_htf_map) mism_htf_time++;

      if(idxH >= 0){
         double dqH = MathAbs(outH[idxH].q - gold[g].q_htf);
         if(dqH > max_abs_q_htf) max_abs_q_htf = dqH;
         if(dqH > EPS_Q) mism_q_htf++;

         int stH = _StateFromStr(gold[g].state_htf);
         int bnH = _BandFromStr(gold[g].band_htf);

         if(outH[idxH].state != stH) mism_state_htf++;
         if(outH[idxH].band  != bnH) mism_band_htf++;
      }

      bool allow_comb = (outL[iL].allow && (idxH>=0 ? outH[idxH].allow : false));
      if((int)allow_comb != gold[g].allow_combined) mism_allow_comb++;
   }

   Print("==== QGate v1.1 Golden Test Summary ====");
   Print("Max |q_ltf diff| = ", DoubleToString(max_abs_q_ltf, 10), " (eps=", DoubleToString(EPS_Q,10), ")");
   Print("Max |q_htf diff| = ", DoubleToString(max_abs_q_htf, 10), " (eps=", DoubleToString(EPS_Q,10), ")");
   Print("LTF mismatches: q=", mism_q_ltf, ", state=", mism_state_ltf, ", band=", mism_band_ltf);
   Print("HTF mismatches: mapped_time=", mism_htf_time, ", q=", mism_q_htf, ", state=", mism_state_htf, ", band=", mism_band_htf);
   Print("Combined allow mismatches: ", mism_allow_comb);

   return INIT_SUCCEEDED;
}

void OnDeinit(const int reason){}
void OnTick(){}
