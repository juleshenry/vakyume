import os
import sys
from vakyume.verifier import Verify
from vakyume.master import PipelineContext


# Mocking enough for Verify to work
class RotaryPistonVane:
    @staticmethod
    def eqn_11_2(
        Q=None,
        Q_0=None,
        Q_external_gas_throughput=None,
        SP_1=None,
        SP_2=None,
        S_vol_pump_speed=None,
        V=None,
        t=None,
        **kwargs,
    ):
        return

    @staticmethod
    def eqn_11_2__t(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        **kwargs,
    ):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        from math import log

        # This is the "messed up" version from test_dummy.py
        t = (
            V
            / S_vol_pump_speed
            * log((SP_1 + (Q_external_gas_throughput + Q_0)) / (SP_2 + (Q + Q_0)))
        )
        return [t]

    @staticmethod
    def eqn_11_2__V(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        t: float,
        **kwargs,
    ):
        raise Exception("Pending Repair")


# Mocking the context for subshards
ctx = PipelineContext("dummy_project")
pyeqn = "t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))"

v = Verify(RotaryPistonVane, pyeqn=pyeqn, subshards_dir=ctx.subshards_dir)
results = v.verify()
import json

print(json.dumps(results, indent=2))
