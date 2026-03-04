from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_4_7__W_cap import eqn_4_7__W
from .eqn_4_7__W_cap_T_cap import eqn_4_7__W_T
from .eqn_4_7__sum_individual_leak_rates import eqn_4_7__sum_individual_leak_rates

class AirLeak:
    eqn_4_7__W = eqn_4_7__W
    eqn_4_7__W_T = eqn_4_7__W_T
    eqn_4_7__sum_individual_leak_rates = eqn_4_7__sum_individual_leak_rates

    @kwasak
    def eqn_4_7(self, W=None, W_T=None, sum_individual_leak_rates=None):
        return
