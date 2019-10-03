#ifndef cycle_measurement_H
#define cycle_measurement_H

int autocorr_periodogram(int start_idx, int end_idx, int avg_window, int min_cycle, int max_cycle,
                         double spectrum_alpha, const double input[],
                         double spectrum[], double cycle_period[], double snr[]);

int mem(int start_idx, int end_idx, int window, int poles, int min_cycle, int max_cycle, double spectrum_alpha,
        int pred_horizon, const double input[], double spectrum[], double cycle_period[], double snr[],
        double mem_pred[]);

#endif
