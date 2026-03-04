from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_24_2__c import eqn_24_2__c
from .eqn_24_2__u_0 import eqn_24_2__u_0
from .eqn_24_2__u_x import eqn_24_2__u_x
from .eqn_24_2__v import eqn_24_2__v
from .eqn_24_2__x import eqn_24_2__x


class TheTheoryOfSpecialRelativity:
    eqn_24_2__c = eqn_24_2__c
    eqn_24_2__u_0 = eqn_24_2__u_0
    eqn_24_2__u_x = eqn_24_2__u_x
    eqn_24_2__v = eqn_24_2__v
    eqn_24_2__x = eqn_24_2__x

    @kwasak
    def eqn_24_2(self, c=None, u_0=None, u_x=None, v=None, x=None):
        """
        FE := force per unit length
        l := length
        x_0 := initial position
        v := velocity
        t := time
        x := position
        x_0 := initial position
        v_0 := initial velocity
        c := speed of light
        x_0 := initial position
        v_0 := initial velocity
        t := time
        x := position
        x_0 := initial position
        v_0 := initial velocity
        c := speed of light
        u_0 := initial velocity
        v := velocity
        x := position
        t := time
        c := speed of light
        u_0 := initial velocity
        v := velocity
        x := position
        t := time
        c := speed of light
        """
        return
