from ta_library import etalib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv('BTCUSDT_3600.csv')
df.set_index('close_time', inplace=True)
cnp = df['close'].to_numpy()

# Decycler is close minus but smoothed, resulting in smooth low lag trend line
# Compare the decycler to just subtracting highpass from close.
# More info chap 4 of Cycle Analytics for Traders by John Ehlers
df['hp'] = etalib.eta_high_pass(12, cnp)
df['decycler'] = etalib.eta_decycler(12, cnp)

(df['close'] - df['hp']).plot(label='Close - HP')
df['decycler'].plot(label='Decycler')
df['close'].plot(label='CLose')

plt.legend()
plt.show()

# A Low period decycler can be very useful as proxy for price when for trading systems.
# It removes a significant amount of noise with very little lag.


