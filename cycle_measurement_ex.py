from ta_library import etalib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Create sinewave to test cycle measurements on
freq = 30  # Cycle period
amp = 1
sample = 500
x = np.arange(sample)
y = np.sin(2 * np.pi * amp * x / freq)
#plt.plot(x, y)
#plt.show()

# Mem on sinewave
mem_period_s, mem_snr_s, mem_pred_s, mem_spectrum_s = etalib.eta_mem(40, 16, 10, 48, 0.2, 4, y)

# Autocorrelation periodogram on sinewave
acp_period_s, acp_snr_s, acp_spectrum_s = etalib.eta_autocorr_periodogram(3, 10, 48, 0.25, y)


df = pd.read_csv('BTCUSDT_3600.csv')
df.set_index('close_time', inplace=True)

cnp = df['close'].to_numpy()

# Use roofing filter to create "sinewave" with cyclic components with periods between 10 and 48
rf = etalib.eta_roofing_filter(48, 10, cnp)
rf_norm = etalib.eta_agc(0.995, rf)
period_out, snr_out, pred_out, spectrum_result = etalib.eta_mem(48, 10, 10, 48, 0.5, 4, rf)
acp_period_p, acp_snr_p, acp_spectrum_p = etalib.eta_autocorr_periodogram(3, 10, 48, 0.25, rf)
