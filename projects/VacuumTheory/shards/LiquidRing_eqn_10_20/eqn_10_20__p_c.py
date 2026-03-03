from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

def eqn_10_20__p_c( self, P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_s: float, **kwargs, ):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    # Solve for p_c:
    R = (S_0 / (S_p)) ** (1.666667)
    # After clearing **0.6: R = (P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) 
    # (P - p_c) = R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(460 + T_i))
    # p_c = P - R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(460 + T_i))
    p_c = P - R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(460 + T_i))
    return [p_c]
