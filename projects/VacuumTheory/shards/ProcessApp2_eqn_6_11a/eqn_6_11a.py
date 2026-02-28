from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_6_11a__A_d_cap import eqn_6_11a__A_d
from .eqn_6_11a__delta_T_cap import eqn_6_11a__delta_T
from .eqn_6_11a__delta_h_i import eqn_6_11a__delta_h_i
from .eqn_6_11a__delta_m import eqn_6_11a__delta_m
from .eqn_6_11a__h_d import eqn_6_11a__h_d
from .eqn_6_11a__m_b import eqn_6_11a__m_b
from .eqn_6_11a__t_R_cap import eqn_6_11a__t_R

class ProcessApp2:
    eqn_6_11a__A_d = staticmethod(eqn_6_11a__A_d)
    eqn_6_11a__delta_T = staticmethod(eqn_6_11a__delta_T)
    eqn_6_11a__delta_h_i = staticmethod(eqn_6_11a__delta_h_i)
    eqn_6_11a__delta_m = staticmethod(eqn_6_11a__delta_m)
    eqn_6_11a__h_d = staticmethod(eqn_6_11a__h_d)
    eqn_6_11a__m_b = staticmethod(eqn_6_11a__m_b)
    eqn_6_11a__t_R = staticmethod(eqn_6_11a__t_R)

    @kwasak_static
    def eqn_6_11a(A_d=None, delta_T=None, delta_h_i=None, delta_m=None, h_d=None, m_b=None, t_R=None, **kwargs):
        return
