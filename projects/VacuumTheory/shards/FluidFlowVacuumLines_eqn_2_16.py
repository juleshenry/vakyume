from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_16(Re=None, f=None, **kwargs):
        return

    @staticmethod
    def eqn_2_16__Re(f: float, **kwargs):
        # [.pyeqn] f = 64 / Re
        result = []
        Re = 64/f
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_16__f(Re: float, **kwargs):
        # [.pyeqn] f = 64 / Re
        result = []
        f = 64/Re
        result.append(f)
        return result

