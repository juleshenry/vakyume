from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_3__N_cap_mfw import eqn_10_3__N_mfw
from .eqn_10_3__Q_cap_gas import eqn_10_3__Q_gas
from .eqn_10_3__T_cap import eqn_10_3__T


class LiquidRing:
    eqn_10_3__N_mfw = eqn_10_3__N_mfw
    eqn_10_3__Q_gas = eqn_10_3__Q_gas
    eqn_10_3__T = eqn_10_3__T

    @kwasak
    def eqn_10_3(self, N_mfw=None, Q_gas=None, T=None):
        return
