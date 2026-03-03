from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class subshards:
    def check_harmony(F_ext, M, a_CM, **kwargs):
        return (F_ext) - (M * a_CM)
