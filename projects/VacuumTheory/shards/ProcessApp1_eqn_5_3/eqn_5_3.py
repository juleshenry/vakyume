from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_5_3__P_0_i import eqn_5_3__P_0_i
from .eqn_5_3__p_i import eqn_5_3__p_i
from .eqn_5_3__x_i import eqn_5_3__x_i


class ProcessApp1:
    eqn_5_3__P_0_i = eqn_5_3__P_0_i
    eqn_5_3__p_i = eqn_5_3__p_i
    eqn_5_3__x_i = eqn_5_3__x_i

    @kwasak
    def eqn_5_3(self, P_0_i=None, p_i=None, x_i=None):
        """
        p_i := component partial pressure
        x_i := liquid component mole fraction
        P_0_i := pure component vapor pressure at equilibrium temperature
        """
        return
