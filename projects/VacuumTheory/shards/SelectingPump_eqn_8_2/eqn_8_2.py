from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_8_2__hp import eqn_8_2__hp
from .eqn_8_2__installed_costs import eqn_8_2__installed_costs

class SelectingPump:
    eqn_8_2__hp = staticmethod(eqn_8_2__hp)
    eqn_8_2__installed_costs = staticmethod(eqn_8_2__installed_costs)

    @kwasak_static
    def eqn_8_2(hp=None, installed_costs=None, **kwargs):
        return
