from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_19__P(
    self,
    S_Th: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_c: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    # Solve for P:
    # Step 1: (S_p / S_Th) ** (1.666666666666667) = (P - p_s)*(460 + T_i) / ( (P - p_c)*(460 + T_e) )
    R = (S_p / (S_Th)) ** (1.666666666666667)
    # Step 2: R * ((460 + T_e)) * (P - p_c) = ((460 + T_i)) * (P - p_s)
    # Step 3: P * (R * ((460 + T_e)) - ((460 + T_i))) = R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s
    # Step 4: P = (R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s) / (R * ((460 + T_e)) - ((460 + T_i)))
    P = (R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s) / (
        R * ((460 + T_e)) - ((460 + T_i))
    )
    return [P]
