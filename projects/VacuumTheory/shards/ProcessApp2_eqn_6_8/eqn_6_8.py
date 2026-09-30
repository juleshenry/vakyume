from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_6_8__C_cap_1 import eqn_6_8__C_1
from .eqn_6_8__C_cap_2 import eqn_6_8__C_2
from .eqn_6_8__T_cap_1 import eqn_6_8__T_1
from .eqn_6_8__T_cap_2 import eqn_6_8__T_2
from .eqn_6_8__c_p import eqn_6_8__c_p
from .eqn_6_8__delta_h_c import eqn_6_8__delta_h_c
from .eqn_6_8__delta_h_v import eqn_6_8__delta_h_v
from .eqn_6_8__delta_t import eqn_6_8__delta_t
from .eqn_6_8__m_b import eqn_6_8__m_b
from .eqn_6_8__w_v import eqn_6_8__w_v


class ProcessApp2:
    eqn_6_8__C_1 = eqn_6_8__C_1
    eqn_6_8__C_2 = eqn_6_8__C_2
    eqn_6_8__T_1 = eqn_6_8__T_1
    eqn_6_8__T_2 = eqn_6_8__T_2
    eqn_6_8__c_p = eqn_6_8__c_p
    eqn_6_8__delta_h_c = eqn_6_8__delta_h_c
    eqn_6_8__delta_h_v = eqn_6_8__delta_h_v
    eqn_6_8__delta_t = eqn_6_8__delta_t
    eqn_6_8__m_b = eqn_6_8__m_b
    eqn_6_8__w_v = eqn_6_8__w_v

    @kwasak
    def eqn_6_8(
        self,
        C_1=None,
        C_2=None,
        T_1=None,
        T_2=None,
        c_p=None,
        delta_h_c=None,
        delta_h_v=None,
        delta_t=None,
        m_b=None,
        w_v=None,
    ):
        return
