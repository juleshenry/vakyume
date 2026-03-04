from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class WorkAndEnergy:
    @kwasak
    def eqn_7_1(self, F=None, W=None, d=None):
        """
        W := work
        F := force
        d := displacement
        x_0 := initial position
        x_1 := final position
        """
        return

    def eqn_7_1__F(self, W: float, d: float, **kwargs):
        # W = F * d
        result = []
        F = W / d
        result.append(F)
        return result

    def eqn_7_1__W(self, F: float, d: float, **kwargs):
        # W = F * d
        result = []
        W = F * d
        result.append(W)
        return result

    def eqn_7_1__d(self, F: float, W: float, **kwargs):
        # W = F * d
        result = []
        d = W / F
        result.append(d)
        return result

    @kwasak
    def eqn_7_2(self, F1=None, F2=None, F3=None, W_tot=None, x=None):
        """
        W_tot := total work
        """
        return

    def eqn_7_2__F1(self, F2: float, F3: float, W_tot: float, x: float, **kwargs):
        # W_tot = F1 * x + F2 * x + F3 * x
        result = []
        F1 = -F2 - F3 + W_tot / x
        result.append(F1)
        return result

    def eqn_7_2__F2(self, F1: float, F3: float, W_tot: float, x: float, **kwargs):
        # W_tot = F1 * x + F2 * x + F3 * x
        result = []
        F2 = -F1 - F3 + W_tot / x
        result.append(F2)
        return result

    def eqn_7_2__F3(self, F1: float, F2: float, W_tot: float, x: float, **kwargs):
        # W_tot = F1 * x + F2 * x + F3 * x
        result = []
        F3 = -F1 - F2 + W_tot / x
        result.append(F3)
        return result

    def eqn_7_2__W_tot(self, F1: float, F2: float, F3: float, x: float, **kwargs):
        # W_tot = F1 * x + F2 * x + F3 * x
        result = []
        W_tot = x * (F1 + F2 + F3)
        result.append(W_tot)
        return result

    def eqn_7_2__x(self, F1: float, F2: float, F3: float, W_tot: float, **kwargs):
        # W_tot = F1 * x + F2 * x + F3 * x
        result = []
        x = W_tot / (F1 + F2 + F3)
        result.append(x)
        return result

    @kwasak
    def eqn_7_6(self, F_=None, W_f=None, d=None):
        """
        W_f := work done by force of friction
        F_k := force of kinetic friction
        d := displacement
        """
        return

    def eqn_7_6__F_(self, W_f: float, d: float, **kwargs):
        # W_f = F_ * d
        result = []
        F_ = W_f / d
        result.append(F_)
        return result

    def eqn_7_6__W_f(self, F_: float, d: float, **kwargs):
        # W_f = F_ * d
        result = []
        W_f = F_ * d
        result.append(W_f)
        return result

    def eqn_7_6__d(self, F_: float, W_f: float, **kwargs):
        # W_f = F_ * d
        result = []
        d = W_f / F_
        result.append(d)
        return result

    @kwasak
    def eqn_7_7(self, W_F=None, W_f=None, W_g=None, W_net=None):
        """
        W_net := net work done
        W_F := work done by applied force
        W_g := work done by force of gravity
        W_f := work done by force of friction
        """
        return

    def eqn_7_7__W_F(self, W_f: float, W_g: float, W_net: float, **kwargs):
        # W_net = W_F + W_g + W_f
        result = []
        W_F = -W_f - W_g + W_net
        result.append(W_F)
        return result

    def eqn_7_7__W_f(self, W_F: float, W_g: float, W_net: float, **kwargs):
        # W_net = W_F + W_g + W_f
        result = []
        W_f = -W_F - W_g + W_net
        result.append(W_f)
        return result

    def eqn_7_7__W_g(self, W_F: float, W_f: float, W_net: float, **kwargs):
        # W_net = W_F + W_g + W_f
        result = []
        W_g = -W_F - W_f + W_net
        result.append(W_g)
        return result

    def eqn_7_7__W_net(self, W_F: float, W_f: float, W_g: float, **kwargs):
        # W_net = W_F + W_g + W_f
        result = []
        W_net = W_F + W_f + W_g
        result.append(W_net)
        return result
