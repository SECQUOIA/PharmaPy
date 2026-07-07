from types import SimpleNamespace

import numpy as np

from PharmaPy.Drying_Model import Drying


def test_material_balance_converts_molar_drying_rate_to_mass_for_saturation():
    dryer = Drying(number_nodes=2, supercrit_names=["nitrogen"])
    dryer.idx_volatiles = np.array([0, 2])
    dryer.porosity = 0.4
    dryer.rho_liq = np.array([800.0, 900.0])
    dryer.dz = np.ones(2)
    dryer.Liquid_1 = SimpleNamespace(mw=np.array([18.0, 28.0, 46.0]))

    satur = np.array([0.6, 0.8])
    temp_gas = np.array([300.0, 305.0])
    temp_sol = np.array([299.0, 304.0])
    y_gas = np.array([
        [0.02, 0.96, 0.02],
        [0.03, 0.94, 0.03],
    ])
    x_liq = np.array([
        [0.25, 0.75],
        [0.40, 0.60],
    ])
    u_gas = np.zeros(2)
    dens_gas = np.ones(2)
    dry_rate = np.array([
        [2.0, 0.0, 4.0],
        [3.0, 0.0, 5.0],
    ])
    inputs = {"mass_frac": np.array([0.01, 0.98, 0.01])}

    dsat_dt, _, _ = dryer.material_balance(
        time=0.0,
        satur=satur,
        temp_gas=temp_gas,
        temp_sol=temp_sol,
        y_gas=y_gas,
        x_liq=x_liq,
        u_gas=u_gas,
        dens_gas=dens_gas,
        dry_rate=dry_rate,
        inputs=inputs,
    )

    mass_rate = (dry_rate[:, dryer.idx_volatiles] *
                 dryer.Liquid_1.mw[dryer.idx_volatiles] / 1000).sum(axis=1)
    expected_dsat_dt = -mass_rate / dryer.rho_liq / dryer.porosity

    np.testing.assert_allclose(dsat_dt, expected_dsat_dt)
