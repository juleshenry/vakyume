from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class Basic:
    @kwasak
    def eqn_1_1(self, x=None, y=None, z=None):
        return

    def eqn_1_1__x(self, y: float, z: float, **kwargs):
        # z = x + y
        result = []
        x = -y + z
        result.append(x)
        return result

    def eqn_1_1__y(self, x: float, z: float, **kwargs):
        # z = x + y
        result = []
        y = -x + z
        result.append(y)
        return result

    def eqn_1_1__z(self, x: float, y: float, **kwargs):
        # z = x + y
        result = []
        z = x + y
        result.append(z)
        return result
