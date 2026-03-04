from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_8_6__K_cap_A_cap import eqn_8_6__K_A
from .eqn_8_6__K_cap_B_cap import eqn_8_6__K_B
from .eqn_8_6__W_cap_net import eqn_8_6__W_net


class PotentialEnergyAndConservationOfEnergy:
    eqn_8_6__K_A = eqn_8_6__K_A
    eqn_8_6__K_B = eqn_8_6__K_B
    eqn_8_6__W_net = eqn_8_6__W_net

    @kwasak
    def eqn_8_6(self, K_A=None, K_B=None, W_net=None):
        """
        W_net := net work done
        K_A := initial kinetic energy
        K_B := final kinetic energy
        """
        return
