#include <stdio.h>
#include "common.h"


int max(int start_idx, int end_idx, int length, const double input[], double out[]){
    // Block taken from TaLib
    int i = start_idx + length;
    int trailing_idx = start_idx;
    int max_idx = -1;


    int f;
    int j;
    double highest;
    double tmp;
    while( i <= end_idx ) {
        tmp = input[i];
        if( max_idx < trailing_idx ) {
            max_idx = trailing_idx;
            highest = input[max_idx];
            j = max_idx;
            while( ++j <= i ) {
                tmp = input[j];
                if( tmp > highest ) {
                    max_idx = j;
                    highest = tmp;
                }
            }
        }else if( tmp >= highest ){
            max_idx = i;
            highest = tmp;
        }
        out[i] = highest;

        i++;
        trailing_idx++;
    }

    return 0;
}

// Variable window max value.
// Windows is an array of periods e.g. dominant cycle period
// Coef is a multiplier of the window. Generally  set to 1.
int vw_max(int start_idx, int end_idx, double coef, const double input[], const double windows[], double out[]){
    // Block taken from TaLib
    int i = start_idx + (int) windows[start_idx] + 1;
    int min_idx = i;
    int trailing_idx = start_idx;
    int max_idx = -1;

    int j, cw;
    double highest;

    while( i <= end_idx ) {
        cw = (int) (coef * windows[i] + 0.5);
        j = i - (cw - 1);

        if (j < min_idx){
            j = min_idx;
        }

        highest = 0;
        while(++j <= i) {
            if(input[j] > highest) {
                highest = input[j];
            }
        }
        out[i] = highest;
        i++;
    }
    return 0;
}

int vw_min(int start_idx, int end_idx, double coef, const double input[], const double windows[], double out[]){
    // Block taken from TaLib
    int i = start_idx + (int) windows[start_idx] + 1;
    int min_idx = i;
    int max_idx = -1;

    int j, cw;
    double lowest;

    while( i <= end_idx ) {
        cw = (int) (coef * windows[i] + 0.5);
        j = i - (cw - 1);

        if (j < min_idx){
            j = min_idx;
        }

        lowest = 99999999;
        while(++j <= i) {
            if(input[j] < lowest) {
                lowest = input[j];
            }
        }

        out[i] = lowest;
        i++;
    }
    return 0;
}

int min(int start_idx, int end_idx, int length, const double input[], double out[]){
    int i = start_idx + length;
    int trailing_idx = start_idx;
    int max_idx = -1;

    int j;
    double lowest;
    double tmp;
    while( i <= end_idx ) {
        /* Set the lowest low */
        tmp = input[i];
        if( max_idx < trailing_idx ) {
            max_idx = trailing_idx;
            lowest = input[max_idx];
            j = max_idx;
            while( ++j <= i ) {
                tmp = input[j];
                if( tmp < lowest ) {
                    max_idx = j;
                    lowest = tmp;
                }
            }
        }else if( tmp <= lowest ){
            max_idx = i;
            lowest = tmp;
        }
        out[i] = lowest;

        i++;
        trailing_idx++;
    }
    return 0;
}

int peaks_and_valleys(int start_idx, int end_idx, const double input[], double peaks[], double valleys[]){
    int i = start_idx + 2;
    peaks[start_idx] = input[start_idx];
    peaks[start_idx + 1] = input[start_idx + 1];
    valleys[start_idx] = input[start_idx];
    valleys[start_idx + 1] = input[start_idx + 1];
    while (i <= end_idx){
        peaks[i] = peaks[i - 1];
        valleys[i] = valleys[i - 1];

        if (input[i - 1] > input[i] && input[i - 1] > input[i - 2]){
            peaks[i] = input[i - 1];
        }

        if (input[i - 1] < input[i] && input[i - 1] < input[i - 2]){
            valleys[i] = input[i - 1];
        }

        i++;
    }

}

