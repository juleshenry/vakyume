from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_3__N_mfw_cap import eqn_10_3__N_mfw
from .eqn_10_3__Q_gas_cap import eqn_10_3__Q_gas
from .eqn_10_3__T_cap import eqn_10_3__T

class LiquidRing:
    eqn_10_3__N_mfw = staticmethod(eqn_10_3__N_mfw)
    eqn_10_3__Q_gas = staticmethod(eqn_10_3__Q_gas)
    eqn_10_3__T = staticmethod(eqn_10_3__T)

    @kwasak_static
    def eqn_10_3(N_mfw=None, Q_gas=None, T=None, **kwargs):
        return
