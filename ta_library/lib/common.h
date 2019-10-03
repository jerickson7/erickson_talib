#ifndef COMMON_H
#define COMMON_H

#define TAU 6.2838
#define PI 3.14159

int max(int start_idx, int end_idx, int length, const double input[], double out[]);
int min(int start_idx, int end_idx, int length, const double input[], double out[]);
int vw_max(int start_idx, int end_idx, double coef, const double input[], const double windows[], double out[]);
int vw_min(int start_idx, int end_idx, double coef, const double input[], const double windows[], double out[]);

int peaks_and_valleys(int start_idx, int end_idx, const double input[], double peaks[], double valleys[]);


int get_index(double arr[], long arr_indices[], long current_idx, int last_available);
int get_next_index(double arr[], long arr_indices[]);


//int rsi(int start_idx, int end_idx, int length, const double input[], double out[]);

#endif
