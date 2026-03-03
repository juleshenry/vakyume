from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_17__H_2_1_cap import eqn_5_17__H_2_1
from .eqn_5_17__H_2_3_cap import eqn_5_17__H_2_3
from .eqn_5_17__H_2_mi_cap import eqn_5_17__H_2_mi
from .eqn_5_17__x_1 import eqn_5_17__x_1
from .eqn_5_17__x_3 import eqn_5_17__x_3

class ProcessApp1:
    eqn_5_17__H_2_1 = eqn_5_17__H_2_1
    eqn_5_17__H_2_3 = eqn_5_17__H_2_3
    eqn_5_17__H_2_mi = eqn_5_17__H_2_mi
    eqn_5_17__x_1 = eqn_5_17__x_1
    eqn_5_17__x_3 = eqn_5_17__x_3

    @kwasak_static
    def eqn_5_17(self, H_2_1=None, H_2_3=None, H_2_mi=None, x_1=None, x_3=None):
        return
