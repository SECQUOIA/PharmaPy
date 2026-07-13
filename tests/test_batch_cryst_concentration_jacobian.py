import sys
from types import ModuleType

import numpy as np
import pytest


pytestmark = pytest.mark.unit


def _stub_assimulo_modules(monkeypatch):
    assimulo = ModuleType("assimulo")

    solvers = ModuleType("assimulo.solvers")
    solvers.CVode = object

    problem = ModuleType("assimulo.problem")
    problem.Explicit_Problem = object

    monkeypatch.setitem(sys.modules, "assimulo", assimulo)
    monkeypatch.setitem(sys.modules, "assimulo.solvers", solvers)
    monkeypatch.setitem(sys.modules, "assimulo.problem", problem)


def _import_batch_cryst(monkeypatch):
    try:
        from PharmaPy.Crystallizers import BatchCryst
    except ModuleNotFoundError as exc:
        if exc.name != "assimulo":
            raise
        _stub_assimulo_modules(monkeypatch)

        from PharmaPy.Crystallizers import BatchCryst

    return BatchCryst


class _Liquid:
    def getDensity(self, temp=None):
        return 1000.0


class _Solid:
    kv = 2.0

    def getDensity(self, temp=None):
        return 1500.0


class _Kinetics:
    prim_nucl = 0.0
    sec_nucl = 0.0
    growth = 2.0e12
    params = {
        "growth": [0.0, 0.0, 1.0],
        "nucl_prim": [0.0, 0.0, 1.0],
        "nucl_sec": [0.0, 0.0, 0.0, 0.0],
    }

    def get_solubility(self, temp, conc):
        return 0.25


def test_batch_cryst_concentration_jacobian_diagonal_matches_material_balance(monkeypatch):
    BatchCryst = _import_batch_cryst(monkeypatch)

    crystallizer = BatchCryst.__new__(BatchCryst)
    crystallizer.num_distr = 4
    crystallizer.num_species = 2
    crystallizer.target_ind = 0
    crystallizer.kron_jtg = np.array([1.0, 0.0])
    crystallizer.controls = {
        "temp": {
            "fun": lambda time: 298.15,
            "args": (),
            "kwargs": {},
        }
    }
    crystallizer.Liquid_1 = _Liquid()
    crystallizer.Solid_1 = _Solid()
    crystallizer._Kinetics = _Kinetics()

    moments = np.array([1.0, 2.0, 4.0, 8.0])
    mass_conc = np.array([0.55, 0.20])
    vol_liq = 2.0
    states = np.concatenate((moments, mass_conc, [vol_liq]))

    jacobian = crystallizer.jac_states(
        time=0.0,
        states=states,
        params=None,
        return_only=False,
    )

    rho_l = crystallizer.Liquid_1.getDensity(temp=298.15)
    rho_c = crystallizer.Solid_1.getDensity(temp=298.15)
    tr = (
        3
        * crystallizer.Solid_1.kv
        * crystallizer.Kinetics.growth
        * moments[2]
        * rho_c
        * (1e-6) ** 3
    )
    dtr_dconc_tg = crystallizer.Kinetics.params["growth"][-1] * tr / (
        mass_conc[crystallizer.target_ind] - crystallizer.Kinetics.get_solubility(298.15, mass_conc)
    )
    first_conc = np.outer(
        crystallizer.kron_jtg - mass_conc / rho_l,
        crystallizer.kron_jtg,
    )
    expected_conc_block = -1 / vol_liq * (
        dtr_dconc_tg * first_conc - tr / rho_l * np.eye(len(mass_conc))
    )

    np.testing.assert_allclose(
        jacobian[crystallizer.num_distr:-1, crystallizer.num_distr:-1],
        expected_conc_block,
        rtol=1e-12,
        atol=1e-12,
    )
