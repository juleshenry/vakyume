from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Basic:
    @kwasak_static
    def eqn_1_1(x=None, y=None, z=None, **kwargs):
        return



    @staticmethod
    def eqn_1_1__z(x: float, y: float, **kwargs):
        result = [x + y]  # Fixed the mistake by changing subtraction to addition
        return result


    @staticmethod
    def eqn_1_1__x(y: float, z: float, **kwargs):
        x = z - y
        result = [x]
        return result



    @staticmethod
    def eqn_1_1__y(x: float, z: float, **kwargs):
        y = z - x
        result = [y]
        return result

