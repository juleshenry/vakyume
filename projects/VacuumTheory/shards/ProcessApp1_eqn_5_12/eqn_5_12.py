from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_5_12__Eff import eqn_5_12__Eff
from .eqn_5_12__N_ES import eqn_5_12__N_ES
from .eqn_5_12__N_t import eqn_5_12__N_t
from .eqn_5_12__T import eqn_5_12__T

class ProcessApp1:
    eqn_5_12__Eff = eqn_5_12__Eff
    eqn_5_12__N_ES = eqn_5_12__N_ES
    eqn_5_12__N_t = eqn_5_12__N_t
    eqn_5_12__T = eqn_5_12__T

    @kwasak
    def eqn_5_12(self, Eff=None, N_ES=None, N_t=None, T=None):
        return
