from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_8_8__P_2(self, P_1: float, adiabatic_power_watts: float, f: float, **kwargs):
    # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
    # Solve for P_2:
    # Step 1: (P_2 / P_1) ** 0.286 - 1 = adiabatic_power_watts / (f / 12)
    # Step 2: (P_2 / P_1) ** 0.286 = adiabatic_power_watts / (f / 12) + 1
    # Step 3: P_2 / P_1 = (adiabatic_power_watts / (f / 12) + 1) ** (1.0 / 0.286)
    # Step 4: P_2 = P_1 * (adiabatic_power_watts / (f / 12) + 1) ** (1.0 / 0.286)
    P_2 = P_1 * (adiabatic_power_watts / (f / 12) + 1) ** (1.0 / 0.286)
    return [P_2]
