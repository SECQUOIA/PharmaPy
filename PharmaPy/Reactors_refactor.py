from PharmaPy.MultiPhaseVessel import MultiPhaseVessel
from Mechanisms import ReactionMechanism
from PharmaPy.DataClasses import *

class _BaseReactor(MultiPhaseVessel):
    def __init__(self, integrator=None, temp_ref=273.15, isothermal=False, reset_states=False, controls={}, h_conv=0, state_events={}, adiabatic=False, jac_type="AD", Phases=None, basis='mass_j', ht_mode="jacket", diam=0, area_base=0):
        super().__init__(integrator, temp_ref, isothermal, reset_states, controls, h_conv, state_events, adiabatic, jac_type, Phases, basis, ht_mode, diam, area_base)
        if isothermal:
            assert adiabatic != 1, "Cannot be isothermal and adiabatic with a reaction present"
    @property
    def RxnKinetics(self):
        raise AttributeError(
            "RxnKinetics is a convenience initializer only. Modify self.reaction_regions directly instead"
            "Use reaction_regions[ind].kinetics if kinetics access are desired."
        )

    @RxnKinetics.setter
    def RxnKinetics(self, instance):
        """Conveniene initializer to instantiate reaction kinetics in their proper regimes"""
        if not isinstance(instance,list):
            instance = [instance]
        self._create_intraphaseProcesses_from_RxnKinetics(instance)
        self._post_RxnKinetics_setter(instance)

    def _create_intraphaseProcesses_from_RxnKinetics(self,RKs):
        '''Creates intraphase_processes from RxnKinetics, assuming everything is in the liquid phase corresponding to the index of the RxnKinetic
        This is because RxnKinetics handles multiple reactions, but assumes a single phase
        If reaction_regions is already defined, this step is ignored
        For reactions in other phases, the user must define reaction_regions directly
        '''
        if len(self.intraphase_processes) >0:
            raise RuntimeError(
                "reaction_regions already defined. "
                "Cannot use RxnKinetics convenience API."
            )
        processes = []
        for i,rk in enumerate(RKs):
            assert len(self.Phases.Liquids)>=i-1, "The number of reaction kinetics must match or be less than the number of liquid phases or you must specify reaction_regions manually"
            region = IntraPhaseProcess(phase=PhaseRef("liquid",i),
                                    mechanism=ReactionMechanism(rk))
            processes.append(region)
        self.intraphase_processes = processes

    def _post_RxnKinetics_setter(self,RKs):
        "Place holder in case future children need special behavior"
        self.nomenclature(overwrite=True)
class BatchReactor(_BaseReactor):

    oper_mode = "batch"

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

class SemiBatchReactor(_BaseReactor):

    oper_mode = "semibatch"

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

class ContinuousReactor(_BaseReactor):

    oper_mode = "continuous"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def build_outlet_stream(self, time=None):

        # ------------------------------------------------------------------
        # Create default outlet if none exists (old PharmaPy behavior)
        # ------------------------------------------------------------------

        if len(self.outlet_connections) == 0:

            if len(self.inlet_connections) == 0:
                raise RuntimeError(
                    "Continuous reactor requires an inlet before "
                    "an outlet stream can be generated."
                )

            connections = []

            for inlet_connection in self.inlet_connections:

                stream = copy.deepcopy(
                    inlet_connection.stream
                )

                connections.append(
                    StreamConnection(
                        stream=stream,
                        phase_mappings=copy.deepcopy(
                            inlet_connection.phase_mappings
                        ),
                        split_fraction=inlet_connection.split_fraction
                    )
                )

            self.outlet_connections = connections

        # ------------------------------------------------------------------
        # Total inlet flow
        # ------------------------------------------------------------------

        total_inlet_flow = 0.0

        for connection in self.inlet_connections:

            for mapping in connection.phase_mappings:

                phase = connection.stream.get_phase_from_ref(
                    mapping.source_phase
                )

                total_inlet_flow += phase.mass_flow

        total_split = sum(
            connection.split_fraction
            for connection in self.outlet_connections
        )

        # ------------------------------------------------------------------
        # Populate outlet streams
        # ------------------------------------------------------------------

        for connection in self.outlet_connections:

            outlet_flow = (
                total_inlet_flow
                * connection.split_fraction
                / total_split
            )

            for mapping in connection.phase_mappings:

                vessel_phase = self.Phases.get_phase_from_ref(
                    mapping.source_phase
                )

                stream_phase = connection.stream.get_phase_from_ref(
                    mapping.sink_phase
                )

                updates = vessel_phase.get_state_dict(
                    self.phase_states[mapping.source_phase]
                )

                updates["mass_flow"] = outlet_flow

                stream_phase.updatePhase(**updates)