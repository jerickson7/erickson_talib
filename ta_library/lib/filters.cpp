#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "filters.h"
#include "common.h"


void get_super_smoother_constants(int period, double* c1, double* c2, double* c3){
    // Neg Sqrt(2) / PI
    double a1 = exp(-1.414 * PI / period);
    *c2 = 2.0 * a1 * cos(1.414 * PI / period);
    *c3 = -1 * a1 * a1;
    *c1 = 1 - *c2 - *c3;
}

// 2 Pole butterworth, of a 2 period average to filter out nyquist freq
int super_smoother(int start_idx, int end_idx, int length, const double input[], double out[]){
    double c1, c2, c3;
    get_super_smoother_constants(length, &c1, &c2, &c3);
    int i = start_idx + 2;
    out[start_idx] = input[start_idx];
    out[start_idx + 1] =  (input[start_idx] + input[start_idx + 1]) / 2.0;
    while (i <= end_idx){
        out[i] = c1 * (input[i] + input[i-1]) / 2 + c2 * out[i-1] + c3 * out[i-2];
        i++;
    }
    return 0;
}

int high_pass_filter(int start_idx, int end_idx, double length, const double input[], double out[]){
    int i = start_idx + 1;
    double ang = TAU / length;
    double alpha = (cos(ang) + sin(ang) - 1) / cos(ang);

    out[start_idx] = 0;
    while (i <= end_idx){
        out[i] = (1 + alpha / 2) * (input[i] - input[i - 1]) + (1 - alpha) * out[i - 1];
        i++;
    }
    return 0;
}

int two_pole_high_pass_filter(int start_idx, int end_idx, double length, const double input[], double out[]){
    double rad = 0.707 * TAU / length;
    double alpha = (cos(rad) + sin(rad) - 1) / cos(rad);
    double alpha_1 = (1 - alpha / 2) * (1 - alpha / 2);
    double alpha_2 = 2 * ( 1 - alpha);
    double alpha_3 = (1 - alpha) * (1 - alpha);

    int i = start_idx + 2;
    out[start_idx] = 0;
    out[start_idx + 1] = 0;
    while (i <= end_idx){
        out[i] = alpha_1 * (input[i] - 2*input[i - 1] + input[i - 2]) + alpha_2 * out[i - 1] - alpha_3 * out[i - 2];
        i++;
    }

}

int roofing_filter(int start_idx, int end_idx, int hp_length, int lp_length, const double input[], double out[]){
    double* hp = (double*) calloc(end_idx - start_idx + 1, sizeof(double));
    two_pole_high_pass_filter(start_idx, end_idx, hp_length, input, hp);
    super_smoother(start_idx, end_idx, lp_length, hp, out);
    free(hp);
    return 0;
}


int band_pass_filter(int start_idx, int end_idx, int length, double bandwidth, const double input[], double out[], double lead[]){

    int i = start_idx;
    double beta = cos(TAU / length);
    double gamma = 1 / cos(TAU * bandwidth / length);
    double alpha1 = gamma - sqrt(gamma * gamma - 1);
    double alpha2 = (cos(1.5 * bandwidth * 6.2838 / length) + sin(1.5 * bandwidth * TAU / length) - 1) / cos(1.5 * bandwidth * TAU / length);

    while (i <= end_idx){
        if (i > start_idx + 3){
            out[i] = 0.5 * (1 - alpha1) * (input[i] - input[i - 2]) + beta * (1 + alpha1) * out[i - 1] - alpha1 * out[i - 2];
            if (i > start_idx + 4){
                lead[i] = (1 + alpha2 / 2) * (out[i] - out[i - 1]) + (1 - alpha2) * lead[i - 1]; // One pole hp
            } else {
                lead[i] = out[i];
            }
        } else {
            out[i] = 0;
            lead[i] = 0;
        }
        i++;
    }
    return 0;
}

// A single pole high pass subtracted from price
int decycler(int start_idx, int end_idx, int length, const double input[], double out[]){
    int i = start_idx + 1;
    double ang = TAU / length;
    double alpha = (cos(ang) + sin(ang) - 1) / cos(ang);

    out[start_idx] = input[start_idx];
    while (i <= end_idx){
        out[i] = (alpha / 2) * (input[i] + input[i - 1]) + (1 - alpha) * out[i - 1];
        i++;
    }
    return 0;
}


int automatic_gain_control(int start_idx, int end_idx, double decay, const double input[], double out[]){
    int i = start_idx;
    double peak = 0.0;
    while (i <= end_idx){
        peak *= decay;
        if (fabs(input[i]) > peak){
            peak = fabs(input[i]);
        }

        if (peak != 0) {
            out[i] = input[i] / peak;
        } else {
            out[i] = 0;
        }
        i++;
    }
    return 0;
}

