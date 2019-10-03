#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include "cycle_measurement.h"
#include "common.h"

int autocorr_periodogram(int start_idx, int end_idx, int avg_window, int min_cycle, int max_cycle,
                         double spectrum_alpha, const double input[],
                         double spectrum[], double cycle_period[], double snr[]){

    double sx, sy, sxx, syy, sxy, spx, sp;
    double x, y;
    double denom;
    double dc; // Dom cycle
    double spectrum_sum;

    int spectrum_width = max_cycle - min_cycle + 1;
    double* correlations = (double*) calloc(max_cycle + 1, sizeof(double));

    double* cos_part = (double*) calloc(spectrum_width, sizeof(double));
    double* sin_part = (double*) calloc(spectrum_width, sizeof(double));
    double* sq_sum = (double*) calloc(spectrum_width, sizeof(double));

    // The un-normalized spectrum
    double* raw_spectrum = (double*) calloc(spectrum_width, sizeof(double));

    double max_pwr = 0.0;

    // Looping variables
    int lag, j, n, i;
    int idx;
    int first_run = 1;

    for (i = start_idx + avg_window + max_cycle; i <= end_idx; i++){

        // Find the autocorrelation at each lag period
        for (lag = 0; lag <= max_cycle; lag++){
            sx = sy = sxx = syy = sxy = 0;
            for (j = 0; j < avg_window; j++){
                x = input[i - j];
                y = input[i - lag - j];
                sx += x;
                sy += y;
                sxx += x * x;
                sxy += x * y;
                syy += y * y;
            }

            denom = (avg_window * sxx - sx * sx) * (avg_window * syy - sy * sy);
            if (denom > 0){
                correlations[lag] = (avg_window * sxy - sx * sy) / sqrt(denom);
            } else {
                correlations[lag] = 0;
            }
        }

        max_pwr *= 0.995; // Decay last max

        // Compute spectrum.
        for (j = min_cycle; j <= max_cycle; j++){
            idx = j - min_cycle;
            cos_part[idx] = 0;
            sin_part[idx] = 0;
            for (n = avg_window; n < max_cycle; n++){
                cos_part[idx] += correlations[n] * cos(TAU * n / j);
                sin_part[idx] += correlations[n] * sin(TAU * n / j);
            }

            sq_sum[idx] = cos_part[idx] * cos_part[idx] + sin_part[idx] * sin_part[idx];
            if (first_run == 1){
                raw_spectrum[idx] = sq_sum[idx] * sq_sum[idx];
            } else {
                raw_spectrum[idx] = spectrum_alpha * sq_sum[idx] * sq_sum[idx] + (1.0 - spectrum_alpha) * raw_spectrum[idx];
            }

            // Find max power for normalization
            if (raw_spectrum[idx] > max_pwr) {
                max_pwr = raw_spectrum[idx];
            }
            spectrum_sum += raw_spectrum[idx];

        }

        // Normalize spectrum using the max
        for (j = min_cycle; j <= max_cycle; j++){
            idx = j - min_cycle;
            spectrum[i * spectrum_width + idx] = raw_spectrum[idx] / max_pwr;
        }

        // Calculate dominant using center of gravity of spectrum
        spx = 0;
        sp = 0;
        for (j = min_cycle; j <= max_cycle; j++){
            idx = j - min_cycle;
            spx += j * spectrum[i * spectrum_width + idx];
            sp += spectrum[i * spectrum_width + idx];
        }

        double dc;
        if (first_run == 1) { // In case sp is 0 on the first run set cycle period to average.
            dc = ((double) (max_cycle + min_cycle)) / 2.0;
        } else {
            dc = cycle_period[i - 1];
        }

        if (sp != 0) {
            dc = spx / sp;
        }

        cycle_period[i] = dc;
        snr[i] = 1 / (spectrum_sum / spectrum_width);
        first_run = 0;
    }

    free(correlations);
    free(cos_part);
    free(sin_part);
    free(sq_sum);
    free(raw_spectrum);

    return 0;
}


// Calculates the n AR coefficients that minimize MSE on input data.
// Used by mem. Returns MSE and Coefficients.
// From numerical recipes
void mem_coef(const double input[], int i0, int i1, int poles, double* mse, double* coef){
    int i, j, k;
    int window = i1 - i0 + 1;

    double p = 0.0;

    double* wk1 = (double*) calloc(window, sizeof(double));
    double* wk2 = (double*) calloc(window, sizeof(double));
    double* wkm = (double*) calloc(poles, sizeof(double));


    for(i=i0; i <= i1; i++){
        p += input[i] * input[i];
    }
    *mse = p / window;
    wk1[0] = input[i0];
    wk2[window - 2] = input[i1];
    for (j=1; j < window - 1; j++){
        wk1[j] = input[i0 + j];
        wk2[j - 1] = input[i0 + j];
    }


    for(k = 0; k < poles; k++){
        double num = 0.0, denom = 0.0;
        for (j = 0; j < (window - k - 1); j++){
            num += (wk1[j] * wk2[j]);
            denom += (wk1[j] * wk1[j] + wk2[j] * wk2[j]);
        }
        coef[k] = 2.0 * num / denom;
        *mse *= (1.0 - coef[k] * coef[k]);
        for (i = 0; i < k; i++){
            coef[i] = wkm[i] - coef[k] * wkm[k - 1 - i];
        }

        if (k == poles - 1){
            free(wk1);
            free(wk2);
            free(wkm);
            return;
        }

        for (i = 0; i <= k; i++) wkm[i] = coef[i];
        for (j = 0; j < (window - k - 2); j++){
            wk1[j] -= (wkm[k] * wk2[j]);
            wk2[j] = wk2[j + 1] - wkm[k] * wk1[j + 1];
        }

    }
    printf("This should never print \n");
}


int mem(int start_idx, int end_idx, int window, int poles, int min_cycle, int max_cycle, double spectrum_alpha,
        int pred_horizon, const double input[], double spectrum[], double cycle_period[], double snr[],
        double mem_pred[]){


    double* coef = (double*) calloc(poles, sizeof(double));
    double* pred_arr = (double*) calloc(poles + pred_horizon, sizeof(double));

    int spectrum_width = max_cycle - min_cycle + 1;

    // Unnormalized spectrum
    double* raw_spectrum = (double*) calloc(spectrum_width, sizeof(double));


    // Loop variables
    int i, j, n;
    int idx;

    double mse, current_period, current_freq, r_sum, i_sum, wr, wi, wpr, wpi, wtemp, theta, sum;
    double pwr, spectrum_sum, max_pwr, real_part, imag_part;
    max_pwr = 0;
    int first_run = 1;

    double spx, sp;

    for (i = start_idx + (window + 1); i <= end_idx; i++){
        mem_coef(input, i - (window - 1), i, poles, &mse, coef);


        // Decay spectrum max.
        max_pwr = 0;

        // Calculate Spectrum at time step i
        for (j = min_cycle; j <= max_cycle; j++){
            idx = j - min_cycle;
            current_freq = 1.0 / j;
            r_sum = 1.0;
            i_sum = 0.0;
            wr = 1.0;
            wi = 0.0;

            theta = TAU * current_freq;
            wpr = cos(theta);
            wpi = sin(theta);
            for (n = 0; n < poles; n++){
                wr = (wtemp = wr) * wpr - wi * wpi;
                wi = wi * wpr + wtemp * wpi;
                r_sum -= coef[n] * wr;
                i_sum -= coef[n] * wi;
            }


            pwr = mse / (r_sum * r_sum + i_sum * i_sum);
            if (first_run == 1){
                raw_spectrum[idx] = pwr;
            } else {
                // Smooth the spectrum using an ema.
                raw_spectrum[idx] = spectrum_alpha * pwr + (1.0 - spectrum_alpha) * raw_spectrum[idx];
            }

            if (pwr > max_pwr){
                max_pwr = pwr;
            }
        }

        // Normalize spectrum to the max ampitude and get spectrum sum;
        spectrum_sum = 0.0;
        for (j=0; j < spectrum_width; j++){
            spectrum[spectrum_width * i + j] = raw_spectrum[j] / max_pwr;
            spectrum_sum += spectrum[spectrum_width * i + j];
        }

        // Calculate DC using center of gravity of spectrum
        spx = 0;
        sp = 0;
        for (j = 0; j < spectrum_width; j++){
            spx += (j + min_cycle) * spectrum[spectrum_width * i + j];
            sp += spectrum[spectrum_width * i + j];
        }


        if (sp != 0) {
            cycle_period[i] = spx / sp;
        } else {
            if (first_run == 1) { // In case sp is 0 on the first run set cycle period to average.
                cycle_period[i] = ((double) (max_cycle + min_cycle)) / 2.0;
            } else {
                cycle_period[i] = cycle_period[i - 1];
            }
        }

        first_run = 0;

        // Calculate Prediction
        // load the last poles into pred array
        for (n = 1; n <= poles; n++){
            pred_arr[n - 1] = input[i - (poles - n)];
        }

        for (n = 0; n < pred_horizon; n++){
            sum = 0;
            for (j = 0; j < poles; j++){
                sum += pred_arr[j + n] * coef[poles - j - 1];
            }
            pred_arr[poles + n] = sum;
        }
        mem_pred[i] = pred_arr[poles + pred_horizon - 1];
        snr[i] = 1 / (spectrum_sum / spectrum_width);

    }
    free(coef);
    free(pred_arr);
    free(raw_spectrum);
    return 0;
}
