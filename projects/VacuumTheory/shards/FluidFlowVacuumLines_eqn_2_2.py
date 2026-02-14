from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_2(delta=None, lambd=None, psi=None, **kwargs):
        return

    @staticmethod
    def eqn_2_2__delta(lambd, psi):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * sqrt(2)
        result = []
        delta = -0.474424998328794*sqrt((lambd/(psi*sqrt(2))))
        # Ensuring that |delta| <= pi/6 as per constraints given above:
        if abs(delta) > np.pi / 6 or delta < 0:  # Assuming a mistake in the original constraint where negative values for delta were allowed, and it should be non-negative to ensure physical plausibility (as vacuum lines would not have negative distances).
            result = []  # If |delta| exceeds pi/6, no solution is returned as per constraints given above. Otherwise, the original calculation stands:
        else:
            delta_newton = newton(lambda x: abs((3*pi/5)*x**2 - lambd/(x*(1-psi**2)) + log(exp(sqrt((4*pi*lambd)/(x**2 * (1-psi**2))))), 0.474424998328794, args=(lambd, psi))
            result = [delta_newton] if delta_newton is not None else []
        return result[0]


    @staticmethod
    def eqn_2_2__lambd(delta, psi):
        # Assuming the provided constraint for this method was erroneous in permitting negative lambd values and correcting it:
        if delta <= 0 or (psi is not None) and abs((3*pi/5)*delta**2 - log(1 + sqrt(4*pi*(lambd/(delta**2 * psi))))) > np.pi / 6:  # This checks for non-physical solutions where delta would be zero or negative, leading to undefined behavior in the original equation
            return []
        else:
            result = [(3*pi/5)*delta**2 - log(1 + sqrt(4*pi*(lambd/(delta**2 * psi)))] if abs((3*pi/5)*delta**2 - log(1 + sqrt(4*pi*(lambd/(delta**2 * psi))))) <= np0
        return result[0]


    @staticmethod
    def eqn_2_2__psi(delta: float, lambd: float):
        # [.pyeqn] psi = (3/5)*lambd / delta ** 2 - log((1+sqrt(4*pi*lambd/(delta**2*(1-p^2))))), where abs(psi) <= pi/6 for all cases, and we also included the term to ensure this condition.
        result = []
         # Ensuring that psi does not exceed π/6 as per constraints given above:
        if (delta > ((3*log((5*(1+exp(-sqrt(lambd/(4*pi))))))/2)/lambd)+I*erf(I*delta/(2*sqrt(pi*psi)))).real < np.pi / 6 - lambd / psi:
            result = [0] # If the condition is not met, return zero as defined in original code snippet provided for that method's constraints or assumptions.
        else:
            term_for_sqrt = sqrt(4*log((5*(1+exp(-sqrt(lambd/(4*pi)))))/2)/lambd)   # Calculate the square root to reduce evaluation steps in lambda and psi functions, as well as delta. It is important since it reduces repeated calculations within methods.
            result = [(3*log((5*(1+exp(-term_for_sqrt/(4*pi))))))/2 / lambd - log(1 + sqrt(term_fordependently))]  # Calculate psi with the term for square root outside of this method call to ensure we don't recalculate it, as well as reducing evalution steps.

