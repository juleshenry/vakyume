from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_2_17__mu(self, L: float, d: float, delta_P: float, q: float, **kwargs):
    # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
    # Solve for mu by rearranging the equation delta_P = 0.105 * mu * L * q / d**4
    try:
        if 'v' in kwargs:
            v = kwargs['v']
            expected_v = float(kwargs['v'])
            calculated_mu = (delta_P * d**4) / (0.105 * L * q)
            return [calculated_mu, abs(expected_v - v)]  # Return the difference between expected and actual velocity if 'v' is provided
        else:
            mu = (delta_P * d**4) / (0.105 * L * q)
            return [mu]
    except ZeroDivisionError:
        raise ValueError("The denominator in the equation cannot be zero.")
