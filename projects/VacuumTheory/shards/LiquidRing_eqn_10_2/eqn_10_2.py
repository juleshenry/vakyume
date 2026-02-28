from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_2__PS_cap import eqn_10_2__PS
from .eqn_10_2__Q_gas_cap import eqn_10_2__Q_gas
from .eqn_10_2__V_cap import eqn_10_2__V
from .eqn_10_2__dP_cap import eqn_10_2__dP
from .eqn_10_2__dt import eqn_10_2__dt

class LiquidRing:
    eqn_10_2__PS = staticmethod(eqn_10_2__PS)
    eqn_10_2__Q_gas = staticmethod(eqn_10_2__Q_gas)
    eqn_10_2__V = staticmethod(eqn_10_2__V)
    eqn_10_2__dP = staticmethod(eqn_10_2__dP)
    eqn_10_2__dt = staticmethod(eqn_10_2__dt)

    @kwasak_static
    def eqn_10_2(PS=None, Q_gas=None, V=None, dP=None, dt=None, **kwargs):
        return
