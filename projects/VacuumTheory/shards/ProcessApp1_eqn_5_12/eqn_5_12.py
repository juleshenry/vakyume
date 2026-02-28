from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_12__Eff_cap import eqn_5_12__Eff
from .eqn_5_12__N_ES_cap import eqn_5_12__N_ES
from .eqn_5_12__N_t_cap import eqn_5_12__N_t
from .eqn_5_12__T_cap import eqn_5_12__T

class ProcessApp1:
    eqn_5_12__Eff = staticmethod(eqn_5_12__Eff)
    eqn_5_12__N_ES = staticmethod(eqn_5_12__N_ES)
    eqn_5_12__N_t = staticmethod(eqn_5_12__N_t)
    eqn_5_12__T = staticmethod(eqn_5_12__T)

    @kwasak_static
    def eqn_5_12(Eff=None, N_ES=None, N_t=None, T=None, **kwargs):
        return
