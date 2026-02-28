from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_11_3__F_s_cap import eqn_11_3__F_s
from .eqn_11_3__t import eqn_11_3__t
from .eqn_11_3__t_c import eqn_11_3__t_c

class RotaryPistonVane:
    eqn_11_3__F_s = staticmethod(eqn_11_3__F_s)
    eqn_11_3__t = staticmethod(eqn_11_3__t)
    eqn_11_3__t_c = staticmethod(eqn_11_3__t_c)

    @kwasak_static
    def eqn_11_3(F_s=None, t=None, t_c=None, **kwargs):
        return
