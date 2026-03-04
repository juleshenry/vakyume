from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_6_4__Q_v import eqn_6_4__Q_v
from .eqn_6_4__delta_h_v import eqn_6_4__delta_h_v
from .eqn_6_4__w_v import eqn_6_4__w_v


class ProcessApp2:
    eqn_6_4__Q_v = eqn_6_4__Q_v
    eqn_6_4__delta_h_v = eqn_6_4__delta_h_v
    eqn_6_4__w_v = eqn_6_4__w_v

    @kwasak
    def eqn_6_4(self, Q_v=None, delta_h_v=None, w_v=None):
        return
