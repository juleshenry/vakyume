from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_4_10__T import eqn_4_10__T
from .eqn_4_10__V import eqn_4_10__V
from .eqn_4_10__del_P import eqn_4_10__del_P
from .eqn_4_10__leakage import eqn_4_10__leakage
from .eqn_4_10__t import eqn_4_10__t

class AirLeak:
    eqn_4_10__T = eqn_4_10__T
    eqn_4_10__V = eqn_4_10__V
    eqn_4_10__del_P = eqn_4_10__del_P
    eqn_4_10__leakage = eqn_4_10__leakage
    eqn_4_10__t = eqn_4_10__t

    @kwasak
    def eqn_4_10(self, T=None, V=None, del_P=None, leakage=None, t=None):
        return
