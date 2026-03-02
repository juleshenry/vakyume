from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_6_9__A_cap import eqn_6_9__A
from .eqn_6_9__dV_dt_cap import eqn_6_9__dV_dt
from .eqn_6_9__delta_P_cap import eqn_6_9__delta_P
from .eqn_6_9__m import eqn_6_9__m
from .eqn_6_9__mu import eqn_6_9__mu
from .eqn_6_9__r import eqn_6_9__r
from .eqn_6_9__r_M_cap import eqn_6_9__r_M

class ProcessApp2:
    eqn_6_9__A = eqn_6_9__A
    eqn_6_9__dV_dt = eqn_6_9__dV_dt
    eqn_6_9__delta_P = eqn_6_9__delta_P
    eqn_6_9__m = eqn_6_9__m
    eqn_6_9__mu = eqn_6_9__mu
    eqn_6_9__r = eqn_6_9__r
    eqn_6_9__r_M = eqn_6_9__r_M

    @kwasak_static
    def eqn_6_9(self, A=None, dV_dt=None, delta_P=None, m=None, mu=None, r=None, r_M=None):
        return
