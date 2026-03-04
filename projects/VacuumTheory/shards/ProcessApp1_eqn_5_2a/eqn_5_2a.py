from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_5_2a__K_cap_1 import eqn_5_2a__K_1
from .eqn_5_2a__K_cap_2 import eqn_5_2a__K_2
from .eqn_5_2a__alpha_1_2 import eqn_5_2a__alpha_1_2


class ProcessApp1:
    eqn_5_2a__K_1 = eqn_5_2a__K_1
    eqn_5_2a__K_2 = eqn_5_2a__K_2
    eqn_5_2a__alpha_1_2 = eqn_5_2a__alpha_1_2

    @kwasak
    def eqn_5_2a(self, K_1=None, K_2=None, alpha_1_2=None):
        return
