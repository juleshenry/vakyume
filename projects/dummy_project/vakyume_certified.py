from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak_static
from vakyume.config import UnsolvedException
import numpy as np

class Basic:
    @kwasak_static
    def eqn_1_1(x=None, y=None, z=None, **kwargs):
        return

    def eqn_1_1__x(y: float, z: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        x = -y + z
        result.append(x)
        return result
    def eqn_1_1__y(x: float, z: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        y = -x + z
        result.append(y)
        return result
    def eqn_1_1__z(x: float, y: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        z = x + y
        result.append(z)
        return result
