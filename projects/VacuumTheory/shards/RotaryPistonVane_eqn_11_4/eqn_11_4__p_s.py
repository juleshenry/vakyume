from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_11_4__p_s(self, p_g, p_v, **kwargs):
    # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
    return [log((p_v + p_g) / p_v)] if p_v != 0 else []
