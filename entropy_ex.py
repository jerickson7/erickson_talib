from ta_library import etalib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_csv('BTCUSDT_3600.csv')
df.set_index('close_time', inplace=True)

# Create binary encoding of closes. 0 = Down Close, 1 = Up close.
diff = df['close'].diff()
encoding = np.where(diff > 0, 1, 0)

# Estimate entropy using word length of 3. Encoding states must be set to 2 because its a binary encoding.
# Be careful of setting the word length too high. Ideally the number of possible words (states ^ word length) should be
# much less than window size.
# Word length 4 -> 2 ^ 4 = 16 < 240. Worth length 10 ->  2 ^ 10 = 1024 > 240.

entropy3 = pd.Series(etalib.eta_entropy(240, 3, 2, encoding))
entropy4 = pd.Series(etalib.eta_entropy(240, 4, 2, encoding))


# The two word lengths are strongly correlated, showing the word length value isn't terribly important.
# Don't optimize it.
print("Correlation: ", entropy3.corr(entropy4))  # Correlation should be close to 1.
entropy3.plot(label='Word Len 3')
entropy4.plot(label='Word Len 4')
plt.legend()
plt.show()

# Entropy of market data will generally be fairly close to 1.  It can be a useful feature.

