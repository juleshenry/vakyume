from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_8_2__hp import eqn_8_2__hp
from .eqn_8_2__installed_costs import eqn_8_2__installed_costs

class SelectingPump:
    eqn_8_2__hp = eqn_8_2__hp
    eqn_8_2__installed_costs = eqn_8_2__installed_costs

    @kwasak
    def eqn_8_2(self, hp=None, installed_costs=None):
        """
        hp:= horse power of pump
        """
        return
