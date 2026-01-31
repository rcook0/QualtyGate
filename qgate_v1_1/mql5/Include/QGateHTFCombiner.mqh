#pragma once

// Arrays must be ascending by time (open times).
// Returns last index j where htfTime[j] <= t, or -1 if none.
static int AlignIndex_NoLookahead(const datetime &htfTime[], const datetime t) {
   int lo = 0;
   int hi = ArraySize(htfTime) - 1;
   int ans = -1;
   while(lo <= hi){
      int mid = (lo + hi) / 2;
      if(htfTime[mid] <= t){
         ans = mid;
         lo = mid + 1;
      } else {
         hi = mid - 1;
      }
   }
   return ans;
}
