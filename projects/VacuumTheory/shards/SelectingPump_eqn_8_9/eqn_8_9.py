from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_8_9__E_cap_j import eqn_8_9__E_j
from .eqn_8_9__E_cap_m import eqn_8_9__E_m
from .eqn_8_9__e import eqn_8_9__e
from .eqn_8_9__r import eqn_8_9__r
from .eqn_8_9__s import eqn_8_9__s


class SelectingPump:
    eqn_8_9__E_j = eqn_8_9__E_j
    eqn_8_9__E_m = eqn_8_9__E_m
    eqn_8_9__e = eqn_8_9__e
    eqn_8_9__r = eqn_8_9__r
    eqn_8_9__s = eqn_8_9__s

    @kwasak
    def eqn_8_9(self, E_j=None, E_m=None, e=None, r=None, s=None):
        """
        E_j:=ejector thermal efficiency
        e:=electrical cost, cents per kWh
        s:=steam cost, dollar per 1000 lb
        E_m:=mechanical pump thermal efficiency
        """
        return
