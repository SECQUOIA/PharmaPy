
"""
Created on Fri July 10 2026

@author: zhillma
Refactored the code by dcasasor
"""
from PharmaPy.Phases import classify_phases, SolidPhase, LiquidPhase, VaporPhase
from PharmaPy.Streams import LiquidStream, SolidStream, VaporStream
from PharmaPy.MixedPhases import Slurry, SlurryStream, MixedPhase, MixedStream

from PharmaPy.ProcessControl_Refactor import Controller
from PharmaPy.Results import DynamicResult



import copy
import numpy as np
import os
from PharmaPy.DataClasses import *



eps = np.finfo(float).eps

class MultiPhaseVessel():
    def __init__(self,integrator=None,temp_ref=273.15,
     isothermal=False, reset_states=False, controller=Controller(), h_conv=0, 
      state_events={},
      adiabatic=False,jac_type="AD",Phases=None,
      basis='mass_j',ht_mode="jacket",diam=0,area_base=0):
      

        if isothermal and controller is not None:
            assert 'global_temp' not in controller.states and 'temp' not in controller.states, "Cannot change the temperature of an isothermal unit"

        self.basis = basis
        self.adiabatic = adiabatic
        self.isothermal = isothermal

        self.jac_type = jac_type
        

        self.controller = controller #TODO ZZ refactor analyze_controls to give a Controls dataclass, an empty one if controls None
        self.oper_mode = None #This is not called within the class, but is used by pharmapy to handle connections (either 'batch' or 'continuous', etc.)
        
        
        # Phase init
        self._phase_connections = []
        if Phases is not None:
            self.Phases = Phases
        self._intraphase_processes = []

        #heat transfer
        self.area_ht = None
        self._Utility = None
        self.ht_mode = ht_mode
        self.h_conv = h_conv
        self.diam = diam
        self.area_base = area_base
        self.temp_ref = temp_ref #enthalpy ref

        #State events
        if state_events is None:
            state_events = []
        self.state_event_list = state_events

        #state initialization, all types
        self._initialize_states(reset_states)
        
        #port initialization
        self._inlet_connections = []
        self._outlet_connections = []

        #Integrator
        self.integrator = integrator

        self.debug=False


    @property
    def Phases(self):
        return self._Phases
    
    @Phases.setter
    def Phases(self, phases):
        self._Phases = MixedPhase(phases)
        self.__initial_phases = copy.deepcopy(self._Phases)
        self._post_set_phases()
        
    

    def _basis_units(self):

        units = {
            "mass_j": "kg",
            "mass_conc": "kg/m3",
            "mole_j": "kmol",
            "mole_conc": "kmol/m3"
        }

        return units[self.basis]
    
    def _material_state_definition(self)->StateVariable:

        return StateVariable(
            name=self.basis,
            dim=self.Phases.num_species,
            units=self._basis_units(),
            index=self.Phases.name_species,
            
        )

    def _post_set_phases(self):
        self.define_material_states()
        self.initialize_defualt_states()
        self.configure_default_connections()
        self.nomenclature() 
    def configure_default_connections(self):
        """Hook for subclasses to create default outlet connections."""
        pass
    def initialize_defualt_states(self):
        self.default_diff_states_from_phases()
    def _initialize_state_collections(self):

        self.phase_states = PhaseStateCollection()

        # States exposed to solver
        self.solver_state_collection = StateCollection()

        # States exposed as outputs
        self.output_state_collection = StateCollection()

    @property
    def Inlet(self):
        raise AttributeError("Inlet is a convenience API for setting inlet_connections")


    @Inlet.setter
    def Inlet(self, inlet):

        inlet= inlet if isinstance(inlet,(list,tuple)) else [inlet]

        self._create_default_connections('inlet_connections',inlet)
    
    @property
    def inlet_connections(self)->list[StreamConnection]:
        return self._inlet_connections


    @inlet_connections.setter
    def inlet_connections(self, connections):

        if not isinstance(connections, (list,tuple)):
            raise TypeError(
                "inlet_connections should be a list or tuple"
            )

        if not all(
            isinstance(c, StreamConnection)
            for c in connections
        ):
            raise TypeError(
                "inlet_connections must contain StreamConnection objects")

        self._inlet_connections = connections
    def _create_default_connections(self,connection_attr, inlet_streams):

        connections = []

        for stream in inlet_streams:

            if not isinstance(stream, MixedStream):
                stream = MixedStream(stream)

            mappings = []

            counts = {}

            for phase in stream:

                phase_type = (phase.phase_family.lower())

                idx = counts.get(phase_type, 0)
                counts[phase_type] = idx + 1

                ref = PhaseRef(phase_type, idx)

                mappings.append(
                    PhaseMapping(
                        source_phase=ref,
                        sink_phase=ref
                    )
                )

            connections.append(
                StreamConnection(
                    stream=stream,
                    phase_mappings=mappings
                )
            )
        setattr(self,connection_attr,connections)

    @property
    def outlet_connections(self)->list[StreamConnection]:
        return self._outlet_connections
    
    @outlet_connections.setter
    def outlet_connections(self,connections):

        if not isinstance(connections, (list,tuple)):
            raise TypeError("outlet_connections should be a list or tuple")
        
        if not all(isinstance(c, StreamConnection) for c in connections):
            raise TypeError(
                "outlet_connections must contain StreamConnection objects")

        self._outlet_connections = connections

    @property
    def Outlet(self):

        if self.outlet_conditions is not None:

            if len(self.outlet_conditions.streams) == 1:
                return self.outlet_conditions.streams[0].stream

            raise AttributeError(
                "Multiple outlet streams exist."
            )

        # Backward compatibility before solve
        if len(self.outlet_connections) == 1:
            return self.outlet_connections[0].stream

        raise AttributeError(
            "Multiple outlet connections exist. "
            "Use outlet_connections instead.")
    
    @Outlet.setter
    def Outlet(self,outlet):
        outlet = outlet if isinstance(outlet,(list,tuple)) else [outlet]
        self._create_default_connections('outlet_connections',outlet)
    

    def define_material_states(self):


        material_state = self._material_state_definition()
        counts = {}
        for i, phase in enumerate(self.Phases):
            phase_type = phase.phase_family.lower()
            idx=counts.get(phase_type,0)
            counts[phase_type] = idx+1
            phase_ref = PhaseRef(
                phase_type=phase_type,
                index=idx
            )

            self.phase_states.add(
                phase_ref,
                copy.deepcopy(material_state)
            )

    def default_diff_states_from_phases(self):
        """
        Mark the default material state for every phase as differential.

        Unit operations that require different behavior should override this
        method or modify `phase_states` after phase initialization.
        """
        # do nothing if the phases are already marked diff
        if any(
            state.state_type == "diff"
            for phasestate in self.phase_states
            for state in self.phase_states[phasestate.phase].states.values()
        ):
            return
        for phase in self.phase_states:
            if phase.state.name==self.basis:
                phase.state.update_variable('state_type','diff')

    @property
    def phase_connections(self):
        return self._phase_connections
    
    @phase_connections.setter
    def phase_connections(self,connections:list):
        if not all([isinstance(c,PhaseConnection) for c in connections]):
            raise TypeError("phase_connections should all be PhaseConnection objects")
        if not isinstance(connections,list):
            raise TypeError("phase_connections is expected to be a list")
        self._phase_connections = connections
        self.nomenclature(overwrite=True)

    @property
    def intraphase_processes(self):
        return self._intraphase_processes
    
    @intraphase_processes.setter
    def intraphase_processes(self, regions):

        if not isinstance(regions, list):
            raise TypeError(
                "reaction_regions is expected to be a list"
            )

        if not all(
            isinstance(r, IntraPhaseProcess)
            for r in regions
        ):
            raise TypeError(
                "reaction_regions should all be "
                "ReactionRegion objects"
            )

        self._intraphase_processes = regions

        self.nomenclature(overwrite=True)
    
    @property
    def Utility(self):
        return self._Utility

    @Utility.setter
    def Utility(self, utility):
        self.u_ht = 1 / (1 / self.h_conv + 1 / utility.h_conv)
        self._Utility = utility
        self.output_state_collection.add(
            StateVariable(
                name="q_ht",
                dim=1,
                units="W",
                state_type="post",
                compute_history=self.compute_qht_history
            ),overwrite=True
        )
        

    
    def complete_state(self, state, time)->dict[StateKey]:

        completed = state.copy()

        controlled = self.controller.compute_states(
            time=time,
            completed_state=completed,
            unit=self,
        )

        completed.update(controlled)

        required_states = (
            self.solver_state_collection.states
            | self.output_state_collection.states
        )

        for key, variable in required_states.items():

            if key in completed:
                continue

            
            # Pull from phases/state sources
            if key == StateKey("global_temp"):
                value = self.get_default_temperature()
            else:
                value = self.get_state_value(key)

            if value is not None:
                completed[key] = value

        return completed
    
    def get_default_temperature(self):
        return self.Phases[0].temp
    
    def get_state_value(self, key):

        # phase-associated state
        if key.phase is not None:

            phase = self.Phases.get_phase_from_ref(key.phase)

            return getattr(phase, key.name, None)

        # unit operation state
        return getattr(self, key.name, None)
            
    def __getattr__(self, name):
        # For Backward compatability 
        # You should not use phase_# explicitly, 
        # everywhere should always iterate over all phases
        #Exception is setting the default phase_connections 
        # since those use the same assumptions as PharmaPy 1.0
        if name.startswith("Liquid_"):
            idx = int(name.split("_")[1]) - 1
            return self.Phases.Liquids[idx]

        if name.startswith("Solid_"):
            idx = int(name.split("_")[1]) - 1
            return self.Phases.Solids[idx]

        if name.startswith("Vapor_"):
            idx = int(name.split("_")[1]) - 1
            return self.Phases.Vapors[idx]

        raise AttributeError(name) #TODO Check if this raises unwanted errors
       
    def _initialize_states(self,reset=False):

        self.reset_states = reset

        self.state_variables = StateCollection()

        self.input_states = StateCollection()

        self.output_states = StateCollection()

        self.elapsed_time = 0

        self.result = None
        self._initialize_state_collections()

        
        
    def reset(self):

        self.Phases = self.__initial_phases

        self.elapsed_time = 0
        self.result = None

        for process in self.intraphase_processes:
            process.mechanism.reset()

        for connection in self.phase_connections:
            if connection.mechanism is not None:
                connection.mechanism.reset()
    
    @property
    def has_energy_balance(self):
        return (
            not self.isothermal
            and "global_temp" not in self.controller.states and "temp" not in self.controller.states
        )
    @property
    def has_utility_balance(self):
        return (
            self.has_energy_balance
            and self.Utility is not None
            and not self.adiabatic
        )
    
    def define_solver_states(self,overwrite=False):

        
        for phase_state in self.phase_states:
            if phase_state.state.state_type!='diff': 
                continue
            state_copy = copy.deepcopy(phase_state.state)
            state_copy.phase = (phase_state.phase)
            self.solver_state_collection.add(state_copy,overwrite)
        
        if self.has_energy_balance:

            self.solver_state_collection.add(
                StateVariable(
                    name="global_temp",
                    dim=1,
                    units="K",
                    state_type='diff'
                ),overwrite
            )
            
        for connection in self.phase_connections:

            if connection.mechanism is not None:
                connection.mechanism.add_solver_state_variables(
                    self.solver_state_collection,overwrite
                )

        for process in self.intraphase_processes:
            process.mechanism.add_solver_state_variables(
                self.solver_state_collection,overwrite
            )

    def define_output_states(self,overwrite=False):

        self.output_state_collection.add(
            StateVariable(
                name="Total_m_in_vessel",
                dim=1,
                units="kg/s",
                state_type="post",
                compute_history=self.compute_outlet_massflow_history
            ),overwrite
        )

        
        for process in self.intraphase_processes:

            phase = self.Phases.get_phase_from_ref(process.phase)
            process.mechanism.add_output_state_variables(
                self.output_state_collection,
                overwrite=overwrite,
                process=process
            )
            self.output_state_collection.add(
                StateVariable(
                    name="mole_conc",
                    dim=phase.num_species,
                    units="kmol/m3",
                    state_type="post",
                    index=phase.name_species,
                    phase=process.phase,
                    compute_history=self.compute_mole_conc_history
                ),
                overwrite
            )

        

        for conn in self.phase_connections:

            if conn.mechanism is not None:
                conn.mechanism.add_output_state_variables(
                    self.output_state_collection,overwrite
                )

    def nomenclature(self,overwrite=False):

        self.define_solver_states(overwrite)
        self.define_output_states(overwrite)
        self.name_states = self.solver_state_collection.names()
        self.dim_states = self.solver_state_collection.dims()

    @staticmethod
    def compute_outlet_massflow_history(state_var,time, solver_history,context):
        """computes mass flow across all outlets for all time
        """
        m_flow = np.zeros_like(time)
        for key,value in solver_history.items():
            if not isinstance(key,StateKey):continue
            if context.basis in key.name:
                m_flow += value.sum(axis=1)

        return m_flow
    
    @staticmethod
    def compute_qht_history(state_var,time,solver_history,context):
        q_ht = context.get_heat_transfer_rate(solver_history[StateKey('global_temp')],context.Utility.temp_in)
        return q_ht
    
    @staticmethod
    def compute_mole_conc_history(
            state_var,
            time,
            solver_history,
            context
        ):

        phase = state_var.phase

        mass_key = StateKey(
            context.basis,
            phase
        )

        mass_history = solver_history[mass_key]

        # shape:
        # time x species

        phase_obj = context.Phases.get_phase_from_ref(
            phase
        )

        mw = phase_obj.mw

        vol = phase_obj.vol


        # kmol/m3
        mole_conc = (
            mass_history
            / mw
            / vol
        )

        return mole_conc
    def update_phases_from_state(self, completed_state):
        global_temp = completed_state.get(StateKey("global_temp"))

        for phase_ref, collection in self.phase_states.phasestates.items():
            phase = self.Phases.get_phase_from_ref(phase_ref)
            updates = {}

            for variable in collection.states.values():
                key = StateKey(variable.name, phase_ref)
                if key in completed_state:
                    updates[variable.name] = completed_state[key]
            
            if StateKey("temp", phase_ref) not in completed_state:
                if global_temp is not None:
                    updates["temp"] = global_temp
            if updates:
                phase.updatePhase(**updates)


        for connection in self.phase_connections:
            mechanism = connection.mechanism
            if mechanism is None:
                continue
            mechanism.update_state(completed_state)

        for process in self.intraphase_processes:
            process.mechanism.update_state(completed_state)

    def pack_state_rates(
            self,
            material_rates=None,
            global_rates=None
        ):

        if material_rates is None:
            material_rates = {}

        if global_rates is None:
            global_rates = {}

        packed = []

        for key,state in self.solver_state_collection.states.items():

            if state.state_type != "diff":
                continue

            if state.phase is not None:

                if key not in material_rates:
                    continue

                packed.extend(np.asarray(material_rates[key]).flatten())

            else:

                if key not in global_rates:
                    continue

                packed.extend(np.asarray(global_rates[key]).flatten())

        return np.asarray(packed)
    def debug_print(self,material_contributions):
        rates = self.sum_material_contributions(material_contributions)
        
        print("packed dm =", rates[self.material_key(PhaseRef("liquid",0))].sum())
        print("Mass balance terms")

        print("inlet species")
        print(material_contributions['inlet'][self.material_key(PhaseRef("liquid",0))])

        print("outlet species")
        print(material_contributions['outlet'][self.material_key(PhaseRef("liquid",0))])

        print("reaction species")
        print(material_contributions['intraphase'][self.material_key(PhaseRef("liquid",0))])

        print("sum inlet   ", material_contributions['inlet'][self.material_key(PhaseRef("liquid",0))].sum())
        print("sum outlet  ", material_contributions['outlet'][self.material_key(PhaseRef("liquid",0))].sum())
        print("sum reaction", material_contributions['intraphase'][self.material_key(PhaseRef("liquid",0))].sum())
    
    def unit_model(self, time, states, params=None, sw=None,
                    mat_bce=False, enrgy_bce=False,alg_bce=False):
        unpacked_state = self.solver_state_collection.unpack(states)
        completed_state = self.complete_state(unpacked_state,time)
        self.update_phases_from_state(completed_state)
        # Balances
        material_rates, material_contributions, aux = self.material_balances(
            time,completed_state)
        

        if mat_bce:
            return self.pack_state_rates(material_rates)
        global_rates = {}
        energy_rates = self.energy_balances(
                time,completed_state, aux)
        global_rates.update(energy_rates)

        # utility_rates = self.utility_energy_balance(
        #     time,completed_state)
        # global_rates.update(utility_rates)

        if enrgy_bce:
            return self.pack_state_rates(global_rates=global_rates)

        balances = self.pack_state_rates(material_rates=material_rates,
                                        global_rates=global_rates)
        assert len(balances) == len(states), (
            f"Returned {len(balances)} derivatives "
            f"for {len(states)} solver states."
        )
        if self.debug:
            self.debug_print(material_contributions)
        self.derivatives = balances
        return balances
        # algebraic_residuals = self.algebraic_balances(
        #     time,
        #     completed_state,
        #     operating_conditions,
        # )

        # if alg_bce:
        #     return algebraic_residuals
        # if not algebraic_residuals:
        #     # legacy behavior
        #     return balances
        # else:
        #     return balances, algebraic_residuals
    def compile_integrator(self, **kwargs):
        return self.integrator.compile_integrator(
            self,
            **kwargs
        )

    def configure_solver(self):
        pass

    def create_solver_init_states(self):
        return self.solver_state_collection.pack(
            self.complete_state({},0))

    def solve_unit(
        self,
        runtime=None,
        time_grid=None,
        **kwargs,
    ):
        return self.integrator.solve(
            self,
            runtime=runtime,
            time_grid=time_grid,
            **kwargs,
        )

    
    def initialize_material_rate_dictionary(self)->dict[StateKey,Any]:

        rates = {}

        for key, state in self.solver_state_collection.states.items():

            if key.phase is None:
                continue

            if key.name != self.basis:
                continue

            rates[key] = np.zeros(state.dim)

        return rates
    def get_operating_conditions(self,time:float,completed_state:dict[StateKey])->tuple[StreamConditions,StreamConditions,dict[OperatingKey]]:
        # --------------------------------------------
        # Controller sees current vessel state
        # --------------------------------------------

        self.controller.observe(time,completed_state,self)

        operating_conditions=self.controller.compute_operating_conditions(
                time,
                completed_state,
                self
            )
        # --------------------------------------------
        # Resolve inlet after possible inlet control
        # --------------------------------------------
        resolved_inlets = self.resolve_inlets(completed_state,operating_conditions)
        # --------------------------------------------
        # Controller observes actual inlet
        # --------------------------------------------
        self.controller.observe(time,completed_state,self,resolved_inlets=resolved_inlets)


        operating_conditions.update(
            self.controller.compute_operating_conditions(
                time,completed_state,self,resolved_inlets=resolved_inlets))


        # Re-resolve inlet in case controller changed it
        resolved_inlets = self.resolve_inlets(completed_state,operating_conditions)

        return resolved_inlets,operating_conditions
        
    def material_balances(
        self,
        time:float,
        completed_state:dict[StateKey],
    ):

        resolved_inlets,operating_conditions = self.get_operating_conditions(time,completed_state)
        # --------------------------------------------
        # Material balances
        # --------------------------------------------
        contributions, aux = self.limit_material_rates(time,completed_state,resolved_inlets,operating_conditions)


        rates = self.sum_material_contributions(contributions)

        return rates, contributions, aux

    def check_negative_inventory(
            self,
            rates,
            completed_state,
            limiter_dt=1.0):

        violations = {}

        for state_key, rate in rates.items():

            inventory = completed_state[
                state_key
            ]

            mask = (
                inventory + rate * limiter_dt
                < -eps*10
            )

            if np.any(mask):
                violations[state_key] = mask

        return violations
    
    def calculate_scale(
            self,
            violations,
            completed_state,
            contributions,
            dt=1.0
        ):

        scalable_terms = [
            "transfer",
            "outlet",
        ]

        fixed_terms = [
            "inlet",
            "intraphase",
        ]

        scales = {}

        # Only phases that violated need correction
        for state_key in violations:

            phase_ref = state_key.phase

            phase = self.Phases.get_phase_from_ref(
                phase_ref
            )

            inventory = completed_state[state_key]

            fixed = np.zeros_like(inventory)
            scalable = np.zeros_like(inventory)

            for term in fixed_terms:
                fixed += contributions[term].get(
                    state_key,
                    0,
                )

            for term in scalable_terms:
                scalable += contributions[term].get(
                    state_key,
                    0,
                )


            # How much inventory remains after unavoidable mechanisms
            allowable = inventory + fixed*dt


            species_scales = np.ones_like(
                inventory,
                dtype=float
            )


            consuming = scalable < 0

            species_scales[consuming] = (
                allowable[consuming]
                /
                (-scalable[consuming]*dt + eps)
            )


            # Any species that is already impossible due to
            # fixed terms should force complete removal of scalable sinks.
            species_scales = np.clip(
                species_scales,
                0.0,
                1.0,
            )


            # A phase outlet/transfer cannot selectively remove
            # only the offending species. It must scale the entire
            # phase movement.
            phase_scale = np.min(
                species_scales
            )


            scales[state_key] = np.full_like(
                inventory,
                phase_scale,
                dtype=float
            )


        return scales
    
    def scale_phase_inventory(
            self,
            scales):

        for state_key, scale_vector in scales.items():

            phase = self.Phases.get_phase_from_ref(
                state_key.phase
            )

            new_mass = (
                phase.mass_j
                * scale_vector
            )

            phase.updatePhase(
                mass_j=new_mass
            )
                
    def restore_effective_inventory(
            self,
            original_inventory:dict[PhaseRef],
            scales,
        ):

        for phase_ref, mass in original_inventory.items():

            phase = self.Phases.get_phase_from_ref(
                phase_ref
            )

            phase.updatePhase(
                **{
                    self.basis:
                        mass * scales[phase_ref]
                }
            )
    
    def limit_material_rates(
                self,
                time,
                completed_state,
                resolved_inlets,
                operating_conditions,
                limiter_dt=1.0,
            ):

            resolved_outlets = self._resolve_outlets(
                completed_state,
                operating_conditions
            )

            contributions, aux = self.calculate_material_contributions(
                time,
                completed_state,
                resolved_inlets,
                resolved_outlets,
            )

            rates = self.sum_material_contributions(
                contributions
            )

            contributions,aux = self.apply_linear_rate_scaling(
                contributions,
                aux,
                rates,
                completed_state,
                limiter_dt
            )

            return contributions, aux
    def apply_linear_rate_scaling(
            self,
            contributions,
            aux,
            rates,
            completed_state,
            dt=1.0,
        ):
        """
        Linearly scale consuming contributions to prevent negative inventory.

        Assumes dt is the expected integration step. This does not modify
        phases and does not regenerate nonlinear mechanisms.
        """

        violations = self.check_negative_inventory(
            rates,
            completed_state,
            limiter_dt=dt
        )

        if not violations:
            return contributions, aux

        scales = self.calculate_scale(
            violations,
            completed_state,
            contributions,
        )

        scalable_terms = [
            "transfer",
            "outlet",
        ]

        for term in scalable_terms:

            for state_key, rate in contributions[term].items():

                if state_key not in scales:
                    continue

                contributions[term][state_key] *= scales[state_key]
        for item in aux["outlet"]:
            phase_ref = item["mapping"].source_phase
            key = self.material_key(phase_ref)

            if key in scales:
                item["species_flow"] *= scales[key]
                item["outlet_phase"].updatePhase(
                    mass_j= getattr(item['outlet_phase'],self.basis) * scales[key]
                )
        return contributions, aux

    def sum_material_contributions(
        self,
        contributions:dict[str,dict[StateKey]]):

        rates = self.initialize_material_rate_dictionary()

        for contribution_name,contribution in contributions.items():

            for state_key, rate in contribution.items():

                rates[state_key] += rate

        return rates
    def material_key(self, phase_ref:PhaseRef):
        return StateKey(self.basis, phase_ref)
    
    def resolve_outlet_flows(
        self,
        operating_conditions:dict[OperatingKey,Any],
    ):

        flows = {}
        for i, connection in enumerate(self.outlet_connections):
            key = OperatingKey("vol_flow",connection=i,port='outlet')
            if key in operating_conditions:
                flows[i] = operating_conditions[key]

        return flows
    
    def get_total_inlet_vol_flow(self,resolved_inlets:StreamConditions)->float:

        total = 0.0
        for inlet in resolved_inlets.streams:
            total += inlet.stream.vol_flow

        return total
    
    def get_phase_operating_conditions(
        self,
        operating_conditions:dict[OperatingKey,Any],
        connection:int,
        phase_ref:PhaseRef,
        port:str
    )->dict[OperatingKey,Any]:

        updates = {}
        port = port.lower()

        for key, value in operating_conditions.items():

            if key.connection != connection:
                continue

            if (key.phase is not None
                and key.phase != phase_ref):
                continue

            if (key.port is not None
                and key.port!=port):
                continue

            updates[key.name] = value

        return updates

    def compute_requested_phase_outlet_flow(
            self,
            vessel_phase,
            total_outlet_flow,
            connection,
        ):
        """Compute desired outlet flow for a vessel phase.

        The default implementation distributes the specified outlet
        connection flow among the mapped phases in proportion to their
        current volume within the vessel.

        If no total outlet flow is specified, no outlet flow is requested.
        """

        if total_outlet_flow is None:
            return 0.0

        total_vessel_flow = sum(
            phase.vol
            for phase in self.Phases
        )

        if total_vessel_flow <= 0:
            return 0.0

        fraction = vessel_phase.vol / total_vessel_flow

        return fraction * total_outlet_flow
    
    def compute_actual_phase_outlet_flow(
        self,
        vessel_phase,
        requested_flow,
    ):
        """
        Compute the physically achievable outlet flow for a vessel phase.

        The default implementation assumes the requested flow is
        achievable. Subclasses may override this to enforce additional
        constraints (e.g. settling, phase disengagement, hydraulics).
        """

        return min(max(requested_flow, 0.0),vessel_phase.vol)
    
    def _resolve_outlets(
        self,
        completed_state:dict,
        operating_conditions:dict[OperatingKey],
    )->StreamConditions:

        resolved = []

        outlet_flows = self.resolve_outlet_flows(operating_conditions)

        for connection_num, connection in enumerate(self.outlet_connections):# iterate over outlet streams

            outlet_stream = copy.deepcopy(connection.stream)

            total_outlet_flow = outlet_flows.get(connection_num)

            for mapping in connection.phase_mappings: #iterate over each phase in that stream

                vessel_phase = self.Phases.get_phase_from_ref(mapping.source_phase)
                outlet_phase = outlet_stream.get_phase_from_ref(mapping.source_phase)
                
                # Default outlet request
                requested_flow = self.compute_requested_phase_outlet_flow(vessel_phase,total_outlet_flow,connection)

                # Controller (or other operating conditions) may override the request
                ops =self.get_phase_operating_conditions(
                    operating_conditions,
                    connection_num,
                    mapping.source_phase,
                    "outlet",
                )
                requested_flow = ops.pop("vol_flow", requested_flow)

                # Apply physical limits once
                actual_flow = self.compute_actual_phase_outlet_flow(vessel_phase,requested_flow)
                updates = vessel_phase.state_dict

                # Amounts are determined from the resolved outlet flow
                for name in outlet_phase.amount_names:
                    updates.pop(name, None)

                # Keep only the preferred composition representation
                for name in outlet_phase.composition_names:
                    if name != outlet_phase.default_composition_name:
                        updates.pop(name, None)

                # Add any remaining operating-condition overrides
                updates.update(
                    self.get_phase_operating_conditions(
                        operating_conditions,
                        connection_num,
                        mapping.source_phase,
                        "outlet",
                    )
                )

                # Physical limit always wins
                updates["vol_flow"] = actual_flow

                outlet_phase.updatePhase(**updates)

            resolved.append(ResolvedStreamConnection(connection,outlet_stream))

        return StreamConditions(resolved)
        
    def resolve_inlets(self,completed_state,operating_conditions)->StreamConditions:

        resolved = []

        for connection_num, connection in enumerate(self.inlet_connections):

            inlet_stream = copy.deepcopy(connection.stream)

            for mapping in connection.phase_mappings:

                inlet_phase = inlet_stream.get_phase_from_ref(mapping.source_phase)

                updates = self.get_phase_operating_conditions(
                    operating_conditions,
                    connection_num,
                    mapping.source_phase,
                    port='inlet'
                )

                inlet_phase.updatePhase(**updates)


            resolved.append(
                ResolvedStreamConnection(
                    connection=connection,
                    stream=inlet_stream,
                )
            )

        return StreamConditions(resolved)
    
    def calculate_material_contributions(self,time,completed_state,resolved_inlets,resolved_outlets)->tuple[dict[str,dict[StateKey,Any]],dict[str:list]]:

        contributions = {
            "inlet": self.initialize_material_rate_dictionary(),
            "intraphase": self.initialize_material_rate_dictionary(),
            "transfer": self.initialize_material_rate_dictionary(),
            "outlet": self.initialize_material_rate_dictionary()
        }
        aux = {
            "inlet": [],
            "intraphase": [],
            "transfer": [],
            "outlet": []
        }

        self.add_inlet_terms(
            contributions["inlet"],
            aux['inlet'],
            time,
            completed_state,
            resolved_inlets
        )
        self.add_intraphase_terms(
            contributions["intraphase"],
            aux['intraphase'],
            time,
            completed_state
        )
        self.add_transfer_terms(
            contributions["transfer"],
            aux['transfer'],
            time,
            completed_state
        )

        self.add_outlet_terms(
            contributions["outlet"],
            aux['outlet'],
            time,
            completed_state,
            resolved_outlets
        )
        if self.debug:
            phase2 = copy.deepcopy(self.Phases[0])
            print(f'Time:{round(time,2)}, rho:{phase2.getDensity()}, vol:{phase2.vol},mass:{phase2.mass},mf:{phase2.mass_frac}')
        return contributions, aux
    
    def add_inlet_terms(
            self,
            rates:dict,
            aux:list,
            time,
            completed_state,
            resolved_inlets
        ):

        for inlet in resolved_inlets.streams:
            stream = inlet.stream
            for mapping in inlet.connection.phase_mappings:
                stream_phase = stream.get_phase_from_ref(mapping.source_phase)

                species_flow = getattr(stream_phase,self.basis+'_flow')
                state_key = self.material_key(mapping.source_phase)
                if state_key not in rates:
                    continue
                rates[state_key] += species_flow

                aux.append({
                    "connection": inlet.connection,
                    "mapping": mapping,
                    "stream_phase": stream_phase,
                    "vessel_phase": self.Phases.get_phase_from_ref(
                        mapping.sink_phase
                    ),
                    "species_flow": species_flow,
                })

        

    def add_outlet_terms(
            self,
            rates,
            aux:list,
            time,
            completed_state,
            resolved_outlets:StreamConditions
        ):

        

        for outlet in resolved_outlets.streams:

            for mapping in outlet.connection.phase_mappings:

                outlet_phase = outlet.stream.get_phase_from_ref(mapping.source_phase)
                species_flow = getattr(outlet_phase,self.basis + "_flow")
                state_key = self.material_key(mapping.source_phase)
                if state_key not in rates:
                    continue
                rates[state_key] -= species_flow
                if self.debug:
                    print(f"Time:{time}, in_vessel:{self.Phases[0].mass_j}, outlet:{outlet_phase.mass_j}")
                    print(f"time:{time}, Outlet props:",
                        outlet_phase.vol_flow,
                        outlet_phase.mass_flow,
                        outlet_phase.getDensity(),
                        outlet_phase.vol_flow*outlet_phase.getDensity(),
                        f"Sum species_flow:{ species_flow.sum()}",
                        )

                aux.append({
                    "outlet": outlet,
                    "mapping": mapping,
                    "outlet_phase": outlet_phase,
                    "species_flow": species_flow,
                })


    
    def add_intraphase_terms(self,
            rates,
            aux,
            time,
            completed_state
        ):
        for process in self.intraphase_processes:
            phase = self.Phases.get_phase_from_ref(process.phase)

            species_rate, mech_aux = process.mechanism.get_material_rate(
                process=process,
                phase=phase,
                time=time,
                completed_state=completed_state
            )

            
            rates[self.material_key(process.phase)] += species_rate

            aux.append(mech_aux)
            
    def add_transfer_terms(
            self,
            rates,
            aux:list,
            time,
            completed_state
        ):

        # temp = completed_state["temp"]

        for connection in self.phase_connections:

            source = self.Phases.get_phase_from_ref(connection.source_phase)

            sink = self.Phases.get_phase_from_ref(connection.sink_phase)
            if connection.mechanism is None:
                connection.mechanism = DirectTransfer()

            if not connection.active_condition(source,sink):
                continue

            species_rate, mech_aux = connection.mechanism.get_material_rate(
                source_phase=source,
                sink_phase=sink,
                connection=connection,
                completed_state=completed_state,
                time=time
            )


            rates[self.material_key(connection.source_phase)] -= species_rate
            rates[self.material_key(connection.sink_phase)] += species_rate

            
            aux.append(mech_aux)
    
    def energy_balances(
            self,
            time,
            completed_state,
            aux
        ):
        """
        Energy contributions are accumulated in SI units (joules).

        Positive contributions add energy to the vessel.
        Negative contributions remove energy from the vessel.
        """
        contributions = {
            "inlet": 0,
            "intraphase": 0,
            "transfer": 0,
            "outlet": 0,
            "utility": 0,
            "mixing":0,
            "shaftwork":0
        }

        self.add_inlet_energy_terms(
            contributions,
            aux["inlet"],
            time,
            completed_state
        )

        self.add_intraphase_energy_terms(
            contributions,
            aux["intraphase"],
            time,
            completed_state
        )

        self.add_transfer_energy_terms(
            contributions,
            aux["transfer"],
            time,
            completed_state
        )

        self.add_outlet_energy_terms(
            contributions,
            aux["outlet"],
            time,
            completed_state
        )

        self.add_utility_energy_terms(
            contributions,
            time,
            completed_state
        )

        self.add_mixing_energy_terms(
            contributions,
            time,
            completed_state
        )

        self.add_shaftwork_energy_terms(
            contributions,
            time,
            completed_state
        )

        qdot = sum(contributions.values())  # [W]
        basis = 'mass' if self.basis=='mass_j' else self.basis

        # getCp resolves through MixedPhase.__getattr__, which mass-weights
        # the per-phase values and so returns a specific heat [J/kg/K]. The
        # balance is m*cp*dT/dt = qdot, so the holdup mass is required to
        # form the total heat capacity; without it dT/dt carries a spurious
        # mass factor and is only correct for a 1 kg holdup.
        specific_heat = self.Phases.getCp(basis=basis)  # [J/kg/K]
        total_mass = self.Phases.mass  # [kg], summed over every phase

        if total_mass <= 0:
            # An empty vessel has no thermal inertia, so the temperature is
            # undefined and we should not attempt to solve for it.
            raise ValueError(
                "Cannot evaluate the energy balance with a total holdup of "
                f"{total_mass} kg. Charge the vessel before solving, or set "
                "isothermal=True to drop the temperature state."
            )

        heat_capacity = total_mass * specific_heat  # [J/K]

        dtemp_dt = qdot / heat_capacity  # [K/s]
        return {StateKey("global_temp"): dtemp_dt}

    def add_inlet_energy_terms(
            self,
            contributions,
            aux,
            time,
            completed_state
        ):
        """Computes enthalpy effects"""
        for inlet in aux:

            phase = inlet["stream_phase"]

            h_in = phase.getEnthalpy(
                phase.temp,
                temp_ref=self.temp_ref,
                total_h=True,
                basis='mass'
            )  # J/kg mixture
            vessel_phase = inlet["vessel_phase"]

            h_vessel = vessel_phase.getEnthalpy(
                vessel_phase.temp,
                temp_ref=self.temp_ref,
                total_h=True,
                basis="mass",
            )

            contributions["inlet"] += (
                inlet["species_flow"] * (h_in - h_vessel)
            ).sum()

    def add_transfer_energy_terms(
            self,
            contributions,
            aux,
            time,
            completed_state
        ):
        
        for transfer in aux:
            connection = transfer['connection']
            contributions['transfer'] += connection.mechanism.get_heat_generation(aux=transfer,
                                                             completed_state=completed_state,
                                                             time=time)
    def add_intraphase_energy_terms(self,
                                  contributions,
                                  aux,
                                  time,
                                  completed_state
                                ):
        for process_aux in aux:
            process = process_aux['process']
            q = process.mechanism.get_heat_generation(aux=process_aux,
                                                    completed_state=completed_state,
                                                    time=time)
            contributions['intraphase'] += q

    def add_outlet_energy_terms(
            self,
            contributions,
            aux,
            time,
            completed_state
        ):
        """Computes enthalpy effects"""
        for outlet in aux:
            phase = outlet["outlet_phase"]
            temp= phase.temp
            h_out = phase.getEnthalpy(
                temp,
                temp_ref=self.temp_ref,
                total_h=True,
                basis='mass'
            )

            contributions["outlet"] -= (outlet["species_flow"] * h_out).sum()

    def get_heat_transfer_temperature(self):
        """The temperature to use to determine heat transfer from the vessel. The default assumption is to use the temperature of the first liquid phase"""
        return self.Phases[0].temp
    
    def add_utility_energy_terms(
            self,
            contributions,
            time,
            completed_state
        ):

        if not self.has_utility_balance:
            return

        temp = self.get_heat_transfer_temperature()
        temp_ht = self.Utility.temp_in#completed_state[StateKey("temp_ht")]

        qdot = self.get_heat_transfer_rate(temp, temp_ht)
                                           
        self.Utility.temp_out = temp_ht +qdot/(self.Utility.mass_flow*self.Utility.cp)

        contributions["utility"] -= qdot

    def add_mixing_energy_terms(self,contributions,time,completed_state):
        "Used to add heat of mixing terms"
        #should be += self.Phases.getHeatOfMixing()
        contributions['mixing']+=0

    def add_shaftwork_energy_terms(self,contributions,time,completed_state):
        contributions['shaftwork']+=0

    def get_heat_transfer_rate(
            self,
            temp,
            temp_ht
        ):

        if self.ht_mode == "coil":
            raise NotImplementedError

        area = self.get_heat_transfer_area()

        return self.u_ht * area * (temp - temp_ht)
    
        
    def get_heat_transfer_area(self):

        liquid = self.Phases.get_phase_from_ref(
            PhaseRef("liquid",0)
        )
        if self.diam <= 0:
            raise ValueError(
                "Heat transfer requires reactor diameter > 0."
            )

        return 4 * liquid.vol / self.diam + self.area_base
    
    
    def build_output_history(
            self,
            time,
            solver_history
        ):

        history = {}

        for key,state in self.output_state_collection.states.items():

            history[key] = state.compute_history(
                time,
                solver_history,
                self)

        return history
    def build_solver_history(
            self,
            time,
            solver_states
        ):

        history = self.solver_state_collection.unpack_history(
            solver_states
        )

        history["time"] = np.asarray(time)

        return history
    
    def update_final_state(self, solver_history):

        final_state = {
            key: value[-1]
            for key, value in solver_history.items()
            if key != "time"
        }

        self.update_phases_from_state(final_state)
        return final_state
    
    def update_final_conditions(self,completed_state,time,solver_history,output_history):
        completed_state = self.complete_state(completed_state,time[-1])
        resolved_inlets,operating_conditions = self.get_operating_conditions(time,completed_state)
        resolved_outlets = self._resolve_outlets(completed_state,operating_conditions)
        self.outlet_conditions =resolved_outlets
        
        self.elapsed_time = time[-1]

        
    def retrieve_results(self, time, solver_states):

        solver_history = self.build_solver_history(time,solver_states)
        output_history = self.build_output_history(time,solver_history)

        completed_state = self.update_final_state(solver_history)
        self.update_final_conditions(completed_state,time,solver_history,output_history)
        self.result = self.build_dynamic_result(time,
                    solver_history,
                    output_history)
        return self.result
    

    def build_dynamic_result(
        self,
        time,
        solver_history,
        output_history,
    ):

        data = {"time": np.asarray(time)}

        data.update(
            self.solver_state_collection.flatten(solver_history)
        )

        data.update(
            self.output_state_collection.flatten(output_history)
        )

        states_di = {}

        for key, state in self.solver_state_collection.states.items():

            states_di[
                self.solver_state_collection.format_key(key)
            ] = state.as_dict()

        fstates_di = {}

        for key, state in self.output_state_collection.states.items():

            fstates_di[
                self.output_state_collection.format_key(key)
            ] = state.as_dict()

        return DynamicResult(
            states_di,
            fstates_di,
            **data
        )
    @property
    def name_species(self):
        #Backward Compatibility
        return self.Phases.name_species

    @property
    def num_species(self):
        #backward compatibility
        return self.Phases.num_species
        
