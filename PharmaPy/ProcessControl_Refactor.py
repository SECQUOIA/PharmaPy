from PharmaPy.DataClasses import (PhaseConnection,PhaseMapping,PhaseRef,PhaseStateCollection,PhaseStateVariable,
                                  IntraPhaseProcess,StateCollection,StateKey,StateVariable,StreamConnection,
                                  OperatingKey)
from typing import Any
import numpy as np


class Controller:

    def __init__(self):
        self.states = {}
        self.operating_conditions = {}

    def reset(self):
        self.states.clear()
        self.operating_conditions.clear()

    def compute_states(
        self,
        time,
        completed_state,
        unit,
    ):

        self.states = {}

        self.update_state(
            time,
            completed_state,
            unit,
        )

        return self.states.copy()


    def observe(
        self,
        time,
        completed_state,
        unit,
        resolved_inlets=None,
        resolved_outlets=None,
    ):
        """
        Update controller measurements/internal state.

        This should not modify operating conditions.
        """

        pass


    def compute_operating_conditions(
        self,
        time,
        completed_state,
        unit,
        resolved_inlets=None,
        resolved_outlets=None,
    )->dict[OperatingKey,Any]:

        self.operating_conditions = {}

        self.actuate(
            time,
            completed_state,
            unit,
            resolved_inlets,
            resolved_outlets,
        )

        return self.operating_conditions.copy()


    def actuate(
        self,
        time,
        completed_state,
        unit,
        resolved_inlets=None,
        resolved_outlets=None,
    ):
        """
        Set manipulated variables through operating_conditions.
        """

        pass


    def update_state(
        self,
        time,
        completed_state,
        unit,
    ):
        pass
class DefaultContinuousVesselVolume(Controller):

    def __init__(
        self,
        target_volume=None,
        K=1e4,
    ):

        super().__init__()

        self.target_volume = target_volume
        self.K = K


    def observe(
        self,
        time,
        completed_state,
        unit,
        resolved_inlets=None,
        resolved_outlets=None,
    ):

        if self.target_volume is None:
            self.target_volume = unit.Phases.vol


    def actuate(
        self,
        time,
        completed_state,
        unit,
        resolved_inlets=None,
        resolved_outlets=None,
    ):

        if resolved_inlets is None:
            return

        inlet_flow = sum(
            stream.stream.vol_flow
            for stream in resolved_inlets.streams
        )

        volume_error = unit.Phases.vol - self.target_volume

        outlet_flow = inlet_flow + self.K * volume_error

        self.operating_conditions[
            OperatingKey(
                "vol_flow",
                connection=0,
                port='outlet'
            )
        ] = max(outlet_flow,0.0)

class TankLevelController(Controller):

    def update_operating_conditions(
            self,
            time,
            completed_state,
            unit,
        ):

        vessel_vol = unit.Phases.vol

        inlet_flow = sum(
            connection.stream.vol_flow
            for connection in unit.inlet_connections
        )

        outlet_flow = (
            0.0
            if vessel_vol < 2
            else inlet_flow
        )

        self.operating_conditions[
            OperatingKey(
                "vol_flow",
                connection=0,
            )
        ] = outlet_flow

class ComplexController(Controller):

    def update_operating_conditions(
            self,
            time,
            completed_state,
            unit,
        ):

        self.update_inlet_conditions(
            time,
            completed_state,
            unit,
        )

        self.update_outlet_conditions(
            time,
            completed_state,
            unit,
        )

    def update_inlet_conditions(self,time,completed_state,unit):
        ...

    def update_outlet_conditions(self,time,completed_state,unit):
        ...
        