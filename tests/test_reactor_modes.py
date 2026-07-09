# -*- coding: utf-8 -*-

import copy
import importlib.util
import json

import numpy as np
import pytest


HAS_ASSIMULO = importlib.util.find_spec("assimulo") is not None

pytestmark = [
    pytest.mark.assimulo,
    pytest.mark.integration,
    pytest.mark.skipif(
        not HAS_ASSIMULO,
        reason="assimulo is not installed; solver-backed integration tests skipped",
    ),
]

if HAS_ASSIMULO:
    from PharmaPy.Kinetics import RxnKinetics
    from PharmaPy.Phases import LiquidPhase
    from PharmaPy.Reactors import CSTR, BatchReactor, PlugFlowReactor
    from PharmaPy.Reactors import SemibatchReactor
    from PharmaPy.Streams import LiquidStream
    from PharmaPy.Utilities import CoolingWater

    SENSITIVITY_REACTORS = [CSTR, SemibatchReactor]
else:
    SENSITIVITY_REACTORS = []


def _load_pfr_config(data_path):
    with open(data_path["integration"] / "pfr_test_constructor_kwargs.json") as f:
        config = json.load(f)

    config = copy.deepcopy(config)
    tau = config["inlet"].pop("tau")
    config["inlet"]["vol_flow"] = config["phase"]["vol"] / tau

    datapath = str(data_path["integration"] / "pfr_test_pure_comp.json")
    config["kinetics"].update({
        "stoich_matrix": [[-1, -1, 1], [0, -1, 1]],
        "k_params": [40 / 60, 10 / 60],
        "ea_params": [2e3, 1e3],
        "delta_hrxn": [-5e3, -2.5e3],
    })
    config["kinetics"]["path"] = datapath

    return config, datapath


def _reactor_objects(data_path, reactor):
    config, datapath = _load_pfr_config(data_path)

    inlet = LiquidStream(datapath, **config["inlet"])
    phase = LiquidPhase(datapath, **config["phase"])
    kinetics = RxnKinetics(**config["kinetics"])
    utility = CoolingWater(**config["utility"])

    reactor.Inlet = inlet
    reactor.Phases = phase
    reactor.Kinetics = kinetics
    reactor.Utility = utility

    return reactor


def test_pfr_solve_steady_reads_inlet_mole_conc(data_path):
    config, _ = _load_pfr_config(data_path)
    reactor = _reactor_objects(
        data_path, PlugFlowReactor(**config["reactor"])
    )

    vol_position, states = reactor.solve_steady(reactor.Liquid_1.vol)

    vol_position = np.asarray(vol_position)

    assert vol_position.size > 1
    assert states.shape[0] == vol_position.size
    assert reactor.concProfSteady.shape[0] == vol_position.size
    assert reactor.Kinetics.num_rxns == 2
    assert reactor.tempProfSteady[-1] > reactor.Inlet.temp
    assert reactor.tempProfSteady[-1] < reactor.temp_ht_steady


@pytest.mark.parametrize("reactor_cls", SENSITIVITY_REACTORS)
def test_sensitivity_mode_refuses_unsupported_reactors(data_path, reactor_cls):
    if reactor_cls is CSTR:
        reactor = reactor_cls()
    else:
        reactor = reactor_cls(vol_tank=0.002)

    reactor = _reactor_objects(data_path, reactor)

    with pytest.raises(NotImplementedError, match="sensitivity.*not supported"):
        reactor.solve_unit(runtime=1, eval_sens=True, verbose=False)


def test_coil_ht_mode_refuses_unsupported_heat_transfer():
    reactor = BatchReactor(isothermal=False, ht_mode="coil")

    with pytest.raises(NotImplementedError, match="coil.*not supported"):
        reactor.heat_transfer(np.array([300.0]), np.array([290.0]), 0.002)


@pytest.mark.parametrize("reactor_cls", SENSITIVITY_REACTORS)
def test_coil_ht_mode_refuses_through_solve_unit(data_path, reactor_cls):
    if reactor_cls is CSTR:
        reactor = reactor_cls(isothermal=False, ht_mode="coil")
    else:
        reactor = reactor_cls(
            vol_tank=0.002, isothermal=False, ht_mode="coil")

    reactor = _reactor_objects(data_path, reactor)

    with pytest.raises(NotImplementedError, match="coil.*not supported"):
        reactor.solve_unit(runtime=1, verbose=False)
