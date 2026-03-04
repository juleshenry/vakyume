from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_3_2__G import eqn_3_2__G
from .eqn_3_2__G_C import eqn_3_2__G_C
from .eqn_3_2__H import eqn_3_2__H
from .eqn_3_2__P import eqn_3_2__P
from .eqn_3_2__rho import eqn_3_2__rho


class PressMgmt:
    eqn_3_2__G = eqn_3_2__G
    eqn_3_2__G_C = eqn_3_2__G_C
    eqn_3_2__H = eqn_3_2__H
    eqn_3_2__P = eqn_3_2__P
    eqn_3_2__rho = eqn_3_2__rho

    @kwasak
    def eqn_3_2(self, G=None, G_C=None, H=None, P=None, rho=None):
        return
