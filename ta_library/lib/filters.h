#ifndef filters_H
#define filters_H

int high_pass_filter(int start_idx, int end_idx, double length, const double input[], double out[]);
int two_pole_high_pass_filter(int start_idx, int end_idx, double length, const double input[], double out[]);
int super_smoother(int start_idx, int end_idx, int length, const double input[], double out[]);
int decycler(int start_idx, int end_idx, int length, const double input[], double out[]);
int automatic_gain_control(int start_idx, int end_idx, double decay, const double input[], double out[]);
int band_pass_filter(int start_idx, int end_idx, int length, double bandwidth, const double input[], double out[], double lead[]);
int roofing_filter(int start_idx, int end_idx, int hp_length, int lp_length, const double input[], double out[]);

#endif
