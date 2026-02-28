from math import *
import numpy as np

def check_harmony(D, L, P_downstream, P_p, P_upstream, mu, q, **kwargs):
    return (q * P_p) - (3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream))
