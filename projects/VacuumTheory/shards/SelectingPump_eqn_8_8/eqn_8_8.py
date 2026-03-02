from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_8_8__P_1_cap import eqn_8_8__P_1
from .eqn_8_8__P_2_cap import eqn_8_8__P_2
from .eqn_8_8__adiabatic_power_watts import eqn_8_8__adiabatic_power_watts
from .eqn_8_8__f import eqn_8_8__f

class SelectingPump:
    eqn_8_8__P_1 = eqn_8_8__P_1
    eqn_8_8__P_2 = eqn_8_8__P_2
    eqn_8_8__adiabatic_power_watts = eqn_8_8__adiabatic_power_watts
    eqn_8_8__f = eqn_8_8__f

    @kwasak_static
    def eqn_8_8(self, P_1=None, P_2=None, adiabatic_power_watts=None, f=None):
        return
