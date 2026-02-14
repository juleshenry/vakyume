from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class RotaryPistonVane:
    @kwasak_static
    def eqn_11_2(Q=None, Q_0=None, Q_external_gas_throughput=None, SP_1=None, SP_2=None, S_vol_pump_speed=None, V=None, t=None, **kwargs):
        return

    @staticmethod
    def eqn_11_2__Q(Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # Error during Sympy solve: Sympy solve failed
        def func(Q):
            # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
            return eval("(V / S_vol_pump_speed * ln( (SP_1 - (x_external_gas_throughput + x_0))/ (SP_2 - (x + x_0)))) - (t)".replace('x', str(Q)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_11_2__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # Error during Sympy solve: Sympy solve failed
        def func(Q_0):
            # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
            return eval("(V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + x))/ (SP_2 - (Q + x)))) - (t)".replace('x', str(Q_0)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_11_2__Q_external_gas_throughput(Q: float, Q_0: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # Error during Sympy solve: Sympy solve failed
        def func(Q_external_gas_throughput):
            # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
            return eval("(V / S_vol_pump_speed * ln( (SP_1 - (x + Q_0))/ (SP_2 - (Q + Q_0)))) - (t)".replace('x', str(Q_external_gas_throughput)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_11_2__SP_1(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # Error during Sympy solve: Sympy solve failed
        def func(SP_1):
            # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
            return eval("(V / S_vol_pump_speed * ln( (x - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))) - (t)".replace('x', str(SP_1)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_11_2__SP_2(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # Error during Sympy solve: Sympy solve failed
        def func(SP_2):
            # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
            return eval("(V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (x - (Q + Q_0)))) - (t)".replace('x', str(SP_2)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_11_2__S_vol_pump_speed(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # Error during Sympy solve: Sympy solve failed
        def func(S_vol_pump_speed):
            # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
            return eval("(V / x * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))) - (t)".replace('x', str(S_vol_pump_speed)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_11_2__V(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # Error during Sympy solve: Sympy solve failed
        def func(V):
            # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
            return eval("(x / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))) - (t)".replace('x', str(V)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_11_2__t(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # Error during Sympy solve: Sympy solve failed
        def func(t):
            # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
            return eval("(V / S_vol_pump_speed * ln( (SP_1 - (Q_exxernal_gas_xhroughpux + Q_0))/ (SP_2 - (Q + Q_0)))) - (x)".replace('x', str(t)))
        raise UnsolvedException("Pending LLM/Manual Repair")

