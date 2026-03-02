from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_10__Suc_Pres_cap import eqn_2_10__Suc_Pres
from .eqn_2_10__delta_P_cap import eqn_2_10__delta_P
from .eqn_2_10__oper_press import eqn_2_10__oper_press

class FluidFlowVacuumLines:
    eqn_2_10__Suc_Pres = eqn_2_10__Suc_Pres
    eqn_2_10__delta_P = eqn_2_10__delta_P
    eqn_2_10__oper_press = eqn_2_10__oper_press

    @kwasak_static
    def eqn_2_10(self, Suc_Pres=None, delta_P=None, oper_press=None):
        return
