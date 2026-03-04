from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_5_1__K_cap_i import eqn_5_1__K_i
from .eqn_5_1__x_i import eqn_5_1__x_i
from .eqn_5_1__y_i import eqn_5_1__y_i

class ProcessApp1:
    eqn_5_1__K_i = eqn_5_1__K_i
    eqn_5_1__x_i = eqn_5_1__x_i
    eqn_5_1__y_i = eqn_5_1__y_i

    @kwasak
    def eqn_5_1(self, K_i=None, x_i=None, y_i=None):
        """
        K_i := volatility
        y_i := component concentration, vapor
        x_i := component concentration, liquid
        """
        return
