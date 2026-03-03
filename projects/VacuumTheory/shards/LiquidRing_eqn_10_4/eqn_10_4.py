from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_4__Q_cap_gas import eqn_10_4__Q_gas
from .eqn_10_4__S_capP_cap_1 import eqn_10_4__SP_1
from .eqn_10_4__S_capP_cap_2 import eqn_10_4__SP_2
from .eqn_10_4__S_cap_p import eqn_10_4__S_p
from .eqn_10_4__V_cap import eqn_10_4__V
from .eqn_10_4__t import eqn_10_4__t

class LiquidRing:
    eqn_10_4__Q_gas = eqn_10_4__Q_gas
    eqn_10_4__SP_1 = eqn_10_4__SP_1
    eqn_10_4__SP_2 = eqn_10_4__SP_2
    eqn_10_4__S_p = eqn_10_4__S_p
    eqn_10_4__V = eqn_10_4__V
    eqn_10_4__t = eqn_10_4__t

    @kwasak
    def eqn_10_4(self, Q_gas=None, SP_1=None, SP_2=None, S_p=None, V=None, t=None):
        return
