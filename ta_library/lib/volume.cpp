#include "volume.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int chaikin_money_flow(int start_idx, int end_idx, int window, const double volume[], const double high[],
                       const double low[], const double close[], double out[]){

    double vol_sum = 0.0;
    double mfv_sum = 0.0;
    int i;
    for (i = start_idx; i < window; i++){
        vol_sum += volume[i];
        mfv_sum += ((close[i] - low[i]) - (high[i] - close[i])) / (high[i] - low[i]) * volume[i];
    }

    i = start_idx + window - 1;
    out[i] = mfv_sum / vol_sum;

    int j;
    for (int i = start_idx + window; i <= end_idx; i++){
        j = i - window;
        vol_sum -= volume[j];
        vol_sum += volume[i];

        if (high[i] > low[i]){
            mfv_sum += ((close[i] - low[i]) - (high[i] - close[i])) / (high[i] - low[i]) * volume[i];
        }

        if (high[j] > low[j]){
            mfv_sum -= ((close[j] - low[j]) - (high[j] - close[j])) / (high[j] - low[j]) * volume[j];
        }
        out[i] = mfv_sum / vol_sum;
    }

    return 0;
}
