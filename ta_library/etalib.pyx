# cython: language_level=3
# distutils: language = c++

cimport numpy as np
import numpy as np
from cython.view cimport array
from numpy import NaN
from cython import boundscheck, wraparound, binding


INT = np.int64
FLOAT = np.float64

ctypedef np.int64_t INT_t
ctypedef np.float64_t FLOAT_t

cdef extern from "numpy/arrayobject.h":
    int PyArray_TYPE(np.ndarray)
    np.ndarray PyArray_EMPTY(int, np.npy_intp*, int, int)
    np.ndarray PyArray_ZEROS(int, np.npy_intp*, int, int)
    int PyArray_FLAGS(np.ndarray)
    np.ndarray PyArray_GETCONTIGUOUS(np.ndarray)

np.import_array()

cdef np.npy_int get_start_idx(np.npy_intp length, double* a1) except -1:
    cdef:
        double val
    for i from 0 <= i < length:
        val = a1[i]
        if val != val:
            continue
        return i
    else:
        raise Exception("inputs are all NaN")

cdef np.npy_int get_start_idx_int(np.npy_intp length, int* a1) except -1:
    cdef:
        double val
    for i from 0 <= i < length:
        val = a1[i]
        if val != val:
            continue
        return i
    else:
        raise Exception("inputs are all NaN")

cdef np.ndarray make_double_array(np.npy_intp length, int lookback):
    cdef:
        np.ndarray outreal
        double* outreal_data
    outreal = PyArray_EMPTY(1, &length, np.NPY_DOUBLE, np.NPY_DEFAULT)
    outreal_data = <double*>outreal.data
    for i from 0 <= i < min(lookback, length):
        outreal_data[i] = NaN

    return outreal

cdef extern from "filters.h":
    int high_pass_filter(int start_idx, int end_idx, double length, const double input[], double out[])
    int two_pole_high_pass_filter(int start_idx, int end_idx, double length, const double input[], double out[])
    int super_smoother(int start_idx, int end_idx, int length, const double input[], double out[])
    int decycler(int start_idx, int end_idx, int length, const double input[], double out[])
    int automatic_gain_control(int start_idx, int end_idx, double decay, const double input[], double out[])
    int band_pass_filter(int start_idx, int end_idx, int length, double bandwidth, const double input[], double out[], double lead[])
    int roofing_filter(int start_idx, int end_idx, int hp_length, int lp_length, const double input[], double out[])

cdef extern from "cycle_measurement.h":
    int autocorr_periodogram(int start_idx, int end_idx, int avg_window, int min_cycle, int max_cycle,
                             double spectrum_alpha, const double input[],
                             double spectrum[], double cycle_period[], double snr[])

    int mem(int start_idx, int end_idx, int window, int poles, int min_cycle, int max_cycle, double spectrum_alpha,
            int pred_horizon, const double input[], double spectrum[], double cycle_period[], double snr[],
            double mem_pred[]);

cdef extern from "common.h":
    int vw_max(int start_idx, int end_idx, double coef, const double input[], const double windows[], double out[])
    int vw_min(int start_idx, int end_idx, double coef, const double input[], const double windows[], double out[])
    int peaks_and_valleys(int start_idx, int end_idx, const double input[], double peaks[], double valleys[])

cdef extern from "entropy.h":
    int entropy(int start_idx, int end_idx, int window, int word_length, int encoding_states,
                const int input[], double out[])

cdef extern from "volume.h":
    int chaikin_money_flow(int start_idx, int end_idx, int window, const double volume[], const double high[],
                           const double low[], const double close[], double out[])

def eta_chaikin_money_flow(int window, np.ndarray volume not None, np.ndarray high not None,
                           np.ndarray low not None, np.ndarray close not None):

    cdef:
        np.npy_intp out_length

    out_length = volume.shape[0]
    start_idx = 0
    end_idx = <int> out_length - 1
    out = make_double_array(out_length, window + 1)
    ret_code = chaikin_money_flow(start_idx, end_idx, window, <double*> volume.data, <double*> high.data,
                                  <double*> low.data, <double*> close.data, <double*> out.data)

    return out

def eta_entropy(int window, int word_length, int encoding_states, np.ndarray real not None):
    cdef:
        np.npy_intp out_length

    out_length = real.shape[0]
    start_idx = get_start_idx_int(out_length, <int*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + window + 1
    ent_out = make_double_array(out_length, lookback)

    real = real.astype(np.int32) # Make it a 32 bit integer, same as
    ret_code = entropy(start_idx, end_idx, window, word_length, encoding_states, <int*> real.data, <double*>ent_out.data)
    return ent_out

def eta_vw_max(double coef, np.ndarray real not None, np.ndarray window not None):
    cdef:
        np.npy_intp out_length
        int ret_code, length

    if real.shape[0] != window.shape[0]:
        raise Exception("ETA - vw_max - Arrays real and periods have different lengths.")

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    periods_idx = get_start_idx(out_length, <double*>window.data)
    end_idx = <int> out_length - 1
    lookback = periods_idx + window[periods_idx] + 1
    max_out = make_double_array(out_length, lookback)

    ret_code = vw_max(max(start_idx, periods_idx), end_idx, coef, <double*>real.data, <double*> window.data, <double*> max_out.data)
    return max_out

def eta_vw_min(double coef, np.ndarray real not None, np.ndarray window not None):
    cdef:
        np.npy_intp out_length
        int ret_code, length

    if real.shape[0] != window.shape[0]:
        raise Exception("ETA - vw_min - Arrays real and periods have different lengths.")

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    periods_idx = get_start_idx(out_length, <double*>window.data)
    end_idx = <int> out_length - 1
    lookback = periods_idx + window[periods_idx] + 1
    min_out = make_double_array(out_length, lookback)

    ret_code = vw_min(max(start_idx, periods_idx), end_idx, coef, <double*>real.data, <double*> window.data, <double*> min_out.data)
    return min_out


def eta_mem(int window, int poles, int min_cycle, int max_cycle, double spectrum_alpha, int prediction_horizon, np.ndarray real not None):
    cdef:
        np.npy_intp out_length
        int ret_code, length

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + window + 1
    period_out = make_double_array(out_length, lookback)
    snr_out = make_double_array(out_length, lookback)
    pred_out = make_double_array(out_length, lookback)

    # 2D Numpy array's are internally stored as 1d arrays. Pass the stride and rows and columns to C.
    spectrum_width = max_cycle - min_cycle + 1
    spectrum_result = np.zeros((out_length, spectrum_width), dtype=np.float64)
    cdef double[:, :] spectrum_view = spectrum_result

    ret_code = mem(start_idx, end_idx, window, poles, min_cycle, max_cycle, spectrum_alpha, prediction_horizon, <double*>real.data,
                    &spectrum_view[0, 0],  # Pointer to first element in array
                    <double*>period_out.data, # Period
                    <double*>snr_out.data, # SNR
                    <double*>pred_out.data) # Prediction

    return period_out, snr_out, pred_out, spectrum_result

def eta_super_smoother(int window, np.ndarray real not None):
    cdef:
        np.npy_intp out_length

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + window + 1
    ss_out = make_double_array(out_length, lookback)
    ret_code = super_smoother(start_idx, end_idx, window, <double*>real.data, <double*>ss_out.data)
    return ss_out

def eta_autocorr_periodogram(int avg_window, int min_cycle, int max_cycle, double spectrum_alpha, np.ndarray real not None):
    cdef:
        np.npy_intp out_length, spectrum_length


    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + avg_window + max_cycle
    period_out = make_double_array(out_length, lookback)
    snr_out = make_double_array(out_length, lookback)

    spectrum_width = max_cycle - min_cycle + 1
    spectrum_result = np.zeros((out_length, spectrum_width), dtype=np.float64)
    cdef double[:, :] spectrum_view = spectrum_result


    ret_code = autocorr_periodogram(start_idx, end_idx, avg_window, min_cycle, max_cycle, spectrum_alpha, <double*>real.data,
                    &spectrum_view[0, 0],  # Pointer to first element in array
                    <double*>period_out.data, # Period
                    <double*>snr_out.data)

    return period_out, snr_out, spectrum_result

def eta_roofing_filter(int hp_window, int lp_window, np.ndarray real not None):
    cdef:
        np.npy_intp out_length

    if lp_window > hp_window:
        print("Low pass window higher to High pass window")

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + hp_window + 1
    rf_out = make_double_array(out_length, lookback)
    ret_code = roofing_filter(start_idx, end_idx, hp_window, lp_window, <double*>real.data, <double*>rf_out.data)
    return rf_out

def eta_agc(double decay, np.ndarray real not None):
    cdef:
        np.npy_intp out_length

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    agc_out = make_double_array(out_length, start_idx)
    ret_code = automatic_gain_control(start_idx, end_idx, decay, <double*>real.data, <double*>agc_out.data)
    return agc_out

def eta_band_pass_filter(int window, double bandwidth, np.ndarray real not None):
    cdef:
        np.npy_intp out_length

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + window + 1
    bp_out = make_double_array(out_length, lookback)
    lead_out = make_double_array(out_length, lookback)
    ret_code = band_pass_filter(start_idx, end_idx, window, bandwidth, <double*>real.data, <double*>bp_out.data, <double*>lead_out.data)
    return bp_out, lead_out

def eta_decycler(int window, np.ndarray real not None):
    cdef:
        np.npy_intp out_length

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + window + 1
    out = make_double_array(out_length, lookback)
    ret_code = decycler(start_idx, end_idx, window, <double*>real.data, <double*>out.data)
    return out

def eta_high_pass(int window, np.ndarray real not None):
    cdef:
        np.npy_intp out_length

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + window + 1
    out = make_double_array(out_length, lookback)
    ret_code = high_pass_filter(start_idx, end_idx, window, <double*>real.data, <double*>out.data)
    return out

def eta_two_pole_high_pass(int window, np.ndarray real not None):
    cdef:
        np.npy_intp out_length

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    lookback = start_idx + window + 1
    out = make_double_array(out_length, lookback)
    ret_code = two_pole_high_pass_filter(start_idx, end_idx, window, <double*>real.data, <double*>out.data)
    return out


def eta_peaks_and_valleys(np.ndarray real not None, **kwargs):
    cdef:
        np.npy_intp out_length

    out_length = real.shape[0]
    start_idx = get_start_idx(out_length, <double*>real.data)
    end_idx = <int> out_length - 1
    peaks = make_double_array(out_length, start_idx)
    valleys = make_double_array(out_length, start_idx)
    ret_code = peaks_and_valleys(start_idx, end_idx, <double*>real.data, <double*>peaks.data, <double*>valleys.data)
    return peaks, valleys

