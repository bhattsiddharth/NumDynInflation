import numpy as np

def ff(x):
	return x*3

f = np.vectorize(ff)