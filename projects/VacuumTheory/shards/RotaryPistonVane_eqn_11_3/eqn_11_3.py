from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_11_3__F_s import eqn_11_3__F_s
from .eqn_11_3__t import eqn_11_3__t
from .eqn_11_3__t_c import eqn_11_3__t_c


class RotaryPistonVane:
    eqn_11_3__F_s = eqn_11_3__F_s
    eqn_11_3__t = eqn_11_3__t
    eqn_11_3__t_c = eqn_11_3__t_c

    @kwasak
    def eqn_11_3(self, F_s=None, t=None, t_c=None):
        """
        t:= actual evacuation time
        t_c:= calculated evacuation time using Eq 10.4
        F_s:= system factor, based on operating experience
        """
        return
