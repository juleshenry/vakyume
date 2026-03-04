from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_7__M import eqn_7_7__M
from .eqn_7_7__P import eqn_7_7__P
from .eqn_7_7__P_i_0 import eqn_7_7__P_i_0
from .eqn_7_7__W_air import eqn_7_7__W_air
from .eqn_7_7__W_i import eqn_7_7__W_i
from .eqn_7_7__epsilon_i import eqn_7_7__epsilon_i
from .eqn_7_7__p_c import eqn_7_7__p_c
from .eqn_7_7__x_i import eqn_7_7__x_i


class Precondensors:
    eqn_7_7__M = eqn_7_7__M
    eqn_7_7__P = eqn_7_7__P
    eqn_7_7__P_i_0 = eqn_7_7__P_i_0
    eqn_7_7__W_air = eqn_7_7__W_air
    eqn_7_7__W_i = eqn_7_7__W_i
    eqn_7_7__epsilon_i = eqn_7_7__epsilon_i
    eqn_7_7__p_c = eqn_7_7__p_c
    eqn_7_7__x_i = eqn_7_7__x_i

    @kwasak
    def eqn_7_7(
        self,
        M=None,
        P=None,
        P_i_0=None,
        W_air=None,
        W_i=None,
        epsilon_i=None,
        p_c=None,
        x_i=None,
    ):
        return
