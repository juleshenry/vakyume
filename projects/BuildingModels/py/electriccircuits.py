from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class ElectricCircuits:
    @kwasak
    def eqn_20_3(self, I1=None, I2=None, I3=None):
        """
        I1 := current 1
        I2 := current 2
        I3 := current 3
        """
        return

    def eqn_20_3__I1(self, I2: float, I3: float, **kwargs):
        # I1 = I2 + I3
        result = []
        I1 = I2 + I3
        result.append(I1)
        return result

    def eqn_20_3__I2(self, I1: float, I3: float, **kwargs):
        # I1 = I2 + I3
        result = []
        I2 = I1 - I3
        result.append(I2)
        return result

    def eqn_20_3__I3(self, I1: float, I2: float, **kwargs):
        # I1 = I2 + I3
        result = []
        I3 = I1 - I2
        result.append(I3)
        return result

    @kwasak
    def eqn_20_4(self, I=None, Reff=None, Vvoltmeter=None):
        """
        Vvoltmeter := voltage measured by voltmeter
        Reff := effective resistance
        I := current
        RV := voltmeter resistance
        """
        return

    def eqn_20_4__I(self, Reff: float, Vvoltmeter: float, **kwargs):
        # Vvoltmeter = I * Reff
        result = []
        I = Vvoltmeter / Reff
        result.append(I)
        return result

    def eqn_20_4__Reff(self, I: float, Vvoltmeter: float, **kwargs):
        # Vvoltmeter = I * Reff
        result = []
        Reff = Vvoltmeter / I
        result.append(Reff)
        return result

    def eqn_20_4__Vvoltmeter(self, I: float, Reff: float, **kwargs):
        # Vvoltmeter = I * Reff
        result = []
        Vvoltmeter = I * Reff
        result.append(Vvoltmeter)
        return result

    @kwasak
    def eqn_20_5(self, C=None, IR=None, Q=None, V=None):
        """
        V := voltage
        I := current
        R := resistance
        Q := charge
        C := capacitance
        """
        return

    def eqn_20_5__C(self, IR: float, Q: float, V: float, **kwargs):
        # V = IR + Q / C
        result = []
        C = -Q / (IR - V)
        result.append(C)
        return result

    def eqn_20_5__IR(self, C: float, Q: float, V: float, **kwargs):
        # V = IR + Q / C
        result = []
        IR = V - Q / C
        result.append(IR)
        return result

    def eqn_20_5__Q(self, C: float, IR: float, V: float, **kwargs):
        # V = IR + Q / C
        result = []
        Q = C * (-IR + V)
        result.append(Q)
        return result

    def eqn_20_5__V(self, C: float, IR: float, Q: float, **kwargs):
        # V = IR + Q / C
        result = []
        V = IR + Q / C
        result.append(V)
        return result
