from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_15(Re: float = None, f: float = None, **kwargs):
        return

    @staticmethod
    def eqn_2_15__Re(f: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        Re = 0.009971220736 / f**4
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_15__f(Re: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        f = 0.316 / Re ** (1 / 4)
        result.append(f)
        return result


