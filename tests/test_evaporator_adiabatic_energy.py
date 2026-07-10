import numpy as np
import pytest


pytest.importorskip("assimulo")

from PharmaPy.Evaporators import ContinuousEvaporator


pytestmark = [pytest.mark.assimulo, pytest.mark.unit]


class _EnthalpySource:
    def __init__(self, enthalpy):
        self.enthalpy = enthalpy

    def getEnthalpy(self, *args, **kwargs):
        return self.enthalpy


def test_adiabatic_energy_residual_includes_vapor_enthalpy():
    evaporator = ContinuousEvaporator.__new__(ContinuousEvaporator)
    evaporator._Inlet = _EnthalpySource(10.0)
    evaporator.Liquid_1 = _EnthalpySource(20.0)
    evaporator.Vapor_1 = _EnthalpySource(30.0)
    evaporator.reflux_ratio = 0
    evaporator.adiabatic = True
    evaporator.vol_tot = 2.0

    result = evaporator.energy_balances(
        time=0.0,
        flow_liq=1.0,
        flow_vap=2.0,
        vol_liq=1.0,
        u_int=40.0,
        temp=350.0,
        x_liq=np.array([1.0]),
        y_vap=np.array([1.0]),
        mol_i=np.array([3.0]),
        mol_liq=3.0,
        mol_vap=4.0,
        pres=5.0,
        u_inputs={
            "mole_flow": 4.0,
            "mole_frac": np.array([1.0]),
            "temp": 300.0,
        },
    )

    np.testing.assert_allclose(result, [-40.0, 130.0])
