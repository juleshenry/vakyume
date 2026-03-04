from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_10__bhp import eqn_10_10__bhp
from .eqn_10_10__bhp_0 import eqn_10_10__bhp_0
from .eqn_10_10__mu import eqn_10_10__mu
from .eqn_10_10__rho import eqn_10_10__rho

class LiquidRing:
    eqn_10_10__bhp = eqn_10_10__bhp
    eqn_10_10__bhp_0 = eqn_10_10__bhp_0
    eqn_10_10__mu = eqn_10_10__mu
    eqn_10_10__rho = eqn_10_10__rho

    @kwasak
    def eqn_10_10(self, bhp=None, bhp_0=None, mu=None, rho=None):
        return
