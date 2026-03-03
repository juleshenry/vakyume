from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_6_1__T_cap_1 import eqn_6_1__T_1
from .eqn_6_1__T_cap_2 import eqn_6_1__T_2
from .eqn_6_1__T_cap_R_cap import eqn_6_1__T_R
from .eqn_6_1__c_p import eqn_6_1__c_p
from .eqn_6_1__del_h_v import eqn_6_1__del_h_v
from .eqn_6_1__w_1 import eqn_6_1__w_1
from .eqn_6_1__w_2 import eqn_6_1__w_2
from .eqn_6_1__w_v import eqn_6_1__w_v

class ProcessApp2:
    eqn_6_1__T_1 = eqn_6_1__T_1
    eqn_6_1__T_2 = eqn_6_1__T_2
    eqn_6_1__T_R = eqn_6_1__T_R
    eqn_6_1__c_p = eqn_6_1__c_p
    eqn_6_1__del_h_v = eqn_6_1__del_h_v
    eqn_6_1__w_1 = eqn_6_1__w_1
    eqn_6_1__w_2 = eqn_6_1__w_2
    eqn_6_1__w_v = eqn_6_1__w_v

    @kwasak
    def eqn_6_1(self, T_1=None, T_2=None, T_R=None, c_p=None, del_h_v=None, w_1=None, w_2=None, w_v=None):
        return
