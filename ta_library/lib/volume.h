#ifndef volume_H
#define volume_H

int chaikin_money_flow(int start_idx, int end_idx, int window, const double volume[], const double high[],
                       const double low[], const double close[], double out[]);

#endif
