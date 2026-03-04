from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_10_2__P_capS_cap import eqn_10_2__PS
from .eqn_10_2__Q_cap_gas import eqn_10_2__Q_gas
from .eqn_10_2__V_cap import eqn_10_2__V
from .eqn_10_2__dP_cap import eqn_10_2__dP
from .eqn_10_2__dt import eqn_10_2__dt

class LiquidRing:
    eqn_10_2__PS = eqn_10_2__PS
    eqn_10_2__Q_gas = eqn_10_2__Q_gas
    eqn_10_2__V = eqn_10_2__V
    eqn_10_2__dP = eqn_10_2__dP
    eqn_10_2__dt = eqn_10_2__dt

    @kwasak
    def eqn_10_2(self, PS=None, Q_gas=None, V=None, dP=None, dt=None):
        return
