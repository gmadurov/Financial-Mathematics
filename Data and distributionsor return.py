import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DF_AMZN = pd.read_csv('AMZN.csv', delimiter=',')

AR_AMZN = np.array(DF_AMZN)

for i, val in enumerate(DF_AMZN.Open):
    print(DF_AMZN.Open.iloc[i])
print(DF_AMZN.head())