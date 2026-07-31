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
            unit
        )->dict[StateKey,Any]:

        self.states ={}

        self.update_state(
            time,
            completed_state,
            unit,
        )

        return self.states.copy()
    def compute_operating_conditions(self,time,completed_state,unit)->dict[OperatingKey,Any]:

         self.operating_conditions = {}
         self.update_operating_conditions(time,completed_state,unit)
         return self.operating_conditions.copy()
    
    def update_operating_conditions(self,time,completed_state,unit):
         pass
    
    def update_state(
            self,
            time,
            completed_state,
            unit,
        ):
            pass
    
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
        