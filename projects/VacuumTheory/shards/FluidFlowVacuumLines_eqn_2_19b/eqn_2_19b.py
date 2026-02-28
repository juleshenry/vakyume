from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_19b__Re_cap import eqn_2_19b__Re
from .eqn_2_19b__h import eqn_2_19b__h
from .eqn_2_19b__mu import eqn_2_19b__mu
from .eqn_2_19b__rho import eqn_2_19b__rho
from .eqn_2_19b__v import eqn_2_19b__v
from .eqn_2_19b__w import eqn_2_19b__w

class FluidFlowVacuumLines:
    eqn_2_19b__Re = staticmethod(eqn_2_19b__Re)
    eqn_2_19b__h = staticmethod(eqn_2_19b__h)
    eqn_2_19b__mu = staticmethod(eqn_2_19b__mu)
    eqn_2_19b__rho = staticmethod(eqn_2_19b__rho)
    eqn_2_19b__v = staticmethod(eqn_2_19b__v)
    eqn_2_19b__w = staticmethod(eqn_2_19b__w)

    @kwasak_static
    def eqn_2_19b(Re=None, h=None, mu=None, rho=None, v=None, w=None, **kwargs):
        return
