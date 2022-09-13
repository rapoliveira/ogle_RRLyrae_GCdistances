import numpy as np

def log_prior(theta):
    m, b, log_f = theta
    if 0 < m < 1 and 0.0 < b < 2 and -3 < log_f < 0:
        return 0.0
    return -np.inf