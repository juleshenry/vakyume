from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_8_4__hp import eqn_8_4__hp
from .eqn_8_4__installed_costs import eqn_8_4__installed_costs

class SelectingPump:
    eqn_8_4__hp = eqn_8_4__hp
    eqn_8_4__installed_costs = eqn_8_4__installed_costs

    @kwasak_static
    def eqn_8_4(self, hp=None, installed_costs=None):
        return
