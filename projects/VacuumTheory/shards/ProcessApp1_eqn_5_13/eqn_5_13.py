from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_13__HETP_cap import eqn_5_13__HETP
from .eqn_5_13__H_p_cap import eqn_5_13__H_p
from .eqn_5_13__N_ES_cap import eqn_5_13__N_ES

class ProcessApp1:
    eqn_5_13__HETP = eqn_5_13__HETP
    eqn_5_13__H_p = eqn_5_13__H_p
    eqn_5_13__N_ES = eqn_5_13__N_ES

    @kwasak_static
    def eqn_5_13(self, HETP=None, H_p=None, N_ES=None):
        return
