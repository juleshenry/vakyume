from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class harmony_checks:
    def check_harmony(P, S_0, S_p, T_e, T_i, p_0, p_c, p_s, **kwargs):
        res = (S_0) - (S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6)
        # Ensure we return a magnitude (abs) to avoid truth-value ambiguity for complex residuals
        return abs(res)
