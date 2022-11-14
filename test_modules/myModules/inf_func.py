import numpy as np

def ff(x):
	return x**2

f = np.vectorize(ff)