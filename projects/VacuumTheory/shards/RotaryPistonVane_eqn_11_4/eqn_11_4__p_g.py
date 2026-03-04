from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_11_4__p_g(self, p_v, __knots, p_s, **kwargs):
    # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
    return [(p_s * (__knots - p_v)) / p_v] if p_v != 0 else []
