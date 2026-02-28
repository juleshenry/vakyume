from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_4_7__W_cap import eqn_4_7__W
from .eqn_4_7__W_T_cap import eqn_4_7__W_T
from .eqn_4_7__sum_individual_leak_rates import eqn_4_7__sum_individual_leak_rates

class AirLeak:
    eqn_4_7__W = staticmethod(eqn_4_7__W)
    eqn_4_7__W_T = staticmethod(eqn_4_7__W_T)
    eqn_4_7__sum_individual_leak_rates = staticmethod(eqn_4_7__sum_individual_leak_rates)

    @kwasak_static
    def eqn_4_7(W=None, W_T=None, sum_individual_leak_rates=None, **kwargs):
        return
