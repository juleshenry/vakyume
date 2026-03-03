from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_5_2b__K_cap_1 import eqn_5_2b__K_1
from .eqn_5_2b__K_cap_2 import eqn_5_2b__K_2
from .eqn_5_2b__x_1 import eqn_5_2b__x_1
from .eqn_5_2b__x_2 import eqn_5_2b__x_2
from .eqn_5_2b__y_1 import eqn_5_2b__y_1
from .eqn_5_2b__y_2 import eqn_5_2b__y_2

class ProcessApp1:
    eqn_5_2b__K_1 = eqn_5_2b__K_1
    eqn_5_2b__K_2 = eqn_5_2b__K_2
    eqn_5_2b__x_1 = eqn_5_2b__x_1
    eqn_5_2b__x_2 = eqn_5_2b__x_2
    eqn_5_2b__y_1 = eqn_5_2b__y_1
    eqn_5_2b__y_2 = eqn_5_2b__y_2

    @kwasak
    def eqn_5_2b(self, K_1=None, K_2=None, x_1=None, x_2=None, y_1=None, y_2=None):
        return
