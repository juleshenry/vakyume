from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_2a__K_1_cap import eqn_5_2a__K_1
from .eqn_5_2a__K_2_cap import eqn_5_2a__K_2
from .eqn_5_2a__alpha_1_2 import eqn_5_2a__alpha_1_2

class ProcessApp1:
    eqn_5_2a__K_1 = staticmethod(eqn_5_2a__K_1)
    eqn_5_2a__K_2 = staticmethod(eqn_5_2a__K_2)
    eqn_5_2a__alpha_1_2 = staticmethod(eqn_5_2a__alpha_1_2)

    @kwasak_static
    def eqn_5_2a(K_1=None, K_2=None, alpha_1_2=None, **kwargs):
        return
