from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_19__S_Th(self, T_k, **kwargs):
    # [.pyeqn] R = (S_p / (S_Th)) ** 0.6
    return (T_k - 273.15) * 9/5 + 32
