from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_3_4__KAPPA_cap import eqn_3_4__KAPPA
from .eqn_3_4__P_cap import eqn_3_4__P
from .eqn_3_4__V_cap import eqn_3_4__V

class PressMgmt:
    eqn_3_4__KAPPA = eqn_3_4__KAPPA
    eqn_3_4__P = eqn_3_4__P
    eqn_3_4__V = eqn_3_4__V

    @kwasak_static
    def eqn_3_4(self, KAPPA=None, P=None, V=None):
        return
