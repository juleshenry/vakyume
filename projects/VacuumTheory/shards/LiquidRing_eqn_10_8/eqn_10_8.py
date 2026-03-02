from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_8__bhp import eqn_10_8__bhp
from .eqn_10_8__c_p import eqn_10_8__c_p
from .eqn_10_8__delta_T_cap import eqn_10_8__delta_T
from .eqn_10_8__delta_h_i import eqn_10_8__delta_h_i
from .eqn_10_8__f_a import eqn_10_8__f_a
from .eqn_10_8__rho import eqn_10_8__rho
from .eqn_10_8__w_i import eqn_10_8__w_i

class LiquidRing:
    eqn_10_8__bhp = eqn_10_8__bhp
    eqn_10_8__c_p = eqn_10_8__c_p
    eqn_10_8__delta_T = eqn_10_8__delta_T
    eqn_10_8__delta_h_i = eqn_10_8__delta_h_i
    eqn_10_8__f_a = eqn_10_8__f_a
    eqn_10_8__rho = eqn_10_8__rho
    eqn_10_8__w_i = eqn_10_8__w_i

    @kwasak_static
    def eqn_10_8(self, bhp=None, c_p=None, delta_T=None, delta_h_i=None, f_a=None, rho=None, w_i=None):
        return
