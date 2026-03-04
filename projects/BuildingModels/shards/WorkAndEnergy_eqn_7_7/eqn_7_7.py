from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_7_7__W_cap_F_cap import eqn_7_7__W_F
from .eqn_7_7__W_cap_f import eqn_7_7__W_f
from .eqn_7_7__W_cap_g import eqn_7_7__W_g
from .eqn_7_7__W_cap_net import eqn_7_7__W_net


class WorkAndEnergy:
    eqn_7_7__W_F = eqn_7_7__W_F
    eqn_7_7__W_f = eqn_7_7__W_f
    eqn_7_7__W_g = eqn_7_7__W_g
    eqn_7_7__W_net = eqn_7_7__W_net

    @kwasak
    def eqn_7_7(self, W_F=None, W_f=None, W_g=None, W_net=None):
        """
        W_net := net work done
        W_F := work done by applied force
        W_g := work done by force of gravity
        W_f := work done by force of friction
        """
        return
