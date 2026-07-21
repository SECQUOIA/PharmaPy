
"""
Created on Fri July 10 2026

@author: zhillma
Refactored the code by dcasasor
"""
from PharmaPy.Phases import classify_phases, SolidPhase, LiquidPhase, VaporPhase
from PharmaPy.Streams import LiquidStream, SolidStream
from PharmaPy.MixedPhases import Slurry, SlurryStream
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem
from assimulo.solvers.sundials import CVodeError
from PharmaPy.Commons import (reorder_sens, plot_sens, trapezoidal_rule,
                              upwind_fvm, high_resolution_fvm,
                              eval_state_events, handle_events,
                              unpack_states, complete_dict_states,
                              flatten_states)

from PharmaPy.ProcessControl import analyze_controls

from PharmaPy.jac_module import numerical_jac, numerical_jac_central, dx_jac_x
from PharmaPy.Connections import get_inputs, get_inputs_new

from PharmaPy.Results import DynamicResult
import PharmaPy.Kinetics as pk



import copy
import string
import numpy as np
import os
from dataclasses import dataclass, field
from typing import Optional, Sequence, Any
from types import MethodType
from collections import OrderedDict
from collections.abc import Callable


eps = np.finfo(float).eps
# gas_ct = 8.314  # J/mol/K
class TransferMechanism:
    def add_state_variables(self,collection,overwrite=False):
        pass
    def add_output_state_variables(self, outputs,overwrite=False):
        pass
    def transfer(self,mass_j):
        return mass_j
class DirectTransfer(TransferMechanism):
    pass

## helpers
def build_transfer_vectors(paths):

    vectors = []

    for path in paths:
        vec = np.asarray(path.species_weights)

        vec = vec / vec.sum()

        vectors.append(vec)

    return vectors

def get_phase_from_ref(
        phases,
        phase_ref
    ):

    candidates = []

    for phase in phases:

        phase_type = (
            phase.__class__.__name__
            .replace("Phase", "")
            .lower()
        )

        if phase_type == phase_ref.phase_type:
            candidates.append(phase)

    return candidates[phase_ref.index]

## Dataclasses
@dataclass(frozen=True)
class PhaseRef:
    phase_type: str
    index: int
    def __post_init__(self):
        object.__setattr__(self, "phase_type", str(self.phase_type).lower())
    def __eq__(self,otherPhaseRef):
        return self.phase_type==otherPhaseRef.phase_type and self.index==otherPhaseRef.index

@dataclass
class PhaseConnection:
    #TODO move this to connections when done
    #active_condition checks the source sink and temp and must return a boolean
    source_phase: PhaseRef
    sink_phase: PhaseRef
    kinetics:pk.CrystKinetics|pk.RxnKinetics
    species_weights: np.ndarray | None = None
    active_condition: callable=lambda source_phase,sink_phase:True
    mechanism:TransferMechanism|None=DirectTransfer()


@dataclass
class InletConnection:
    source_phase: PhaseRef
    sink_phase: PhaseRef

@dataclass
class ReactionRegion:
    phase:PhaseRef
    kinetics:pk.RxnKinetics

@dataclass
class StateVariable:
    name: str
    dim: int
    units: str
    state_type: str = "post_calc"
    index: Optional[Sequence] = None
    depends_on: tuple = ("time",)
    stream:Optional[str]=None
    phase: Optional[PhaseRef] = None
    compute_history:Callable[[Any,np.ndarray, dict, Any], np.ndarray] | None = None 
    """
    Parameters
    ----------
    name
        Name of the state variable.
    dim
        Number of dimensions.
    units
        Physical units.
    compute_history : callable, optional

        Function with signature

            compute_history(state_var,
                            time,
                            solver_history,
                            context=None)

        returning the complete history of the state.
    """
    def __post_init__(self):
        if self.compute_history is None:
            self.compute_history = self.default_history
        else:
            self.compute_history = MethodType(self.compute_history, self)
    def as_dict(self):
        """Backward compatibility."""
        out = {
            "dim": self.dim,
            "units": self.units,
            "type": self.state_type,
            "depends_on": list(self.depends_on)
        }

        if self.index is not None:
            out["index"] = self.index

        return out
    def update_variable(self,variable_name,new_value):
        setattr(self,variable_name,new_value)

    def default_history(self,time,solver_history,context=None):
        try:
            return solver_history[self.name]
        except KeyError:
            raise KeyError(f"'{self.name}' is not present in the solver history.")
    

@dataclass
class StateCollection:
    states: dict[str, StateVariable] = field(default_factory=dict)

    def _same_phase(self, existing, new):

        if existing.phase is None and new.phase is None:
            return True

        if existing.phase is None:
            return True

        if new.phase is None:
            return True

        return existing.phase == new.phase
    def add(self, state: StateVariable,overwrite=False,error_on_conflict=False):

        existing = self.states.get(state.name)

        if existing is None:
            self.states[state.name] = state
            return

        same_phase = self._same_phase(existing,state)

        if not same_phase:
            self.states[state.name] = state
            return

        same = state == existing

        if same:
            return

        if overwrite:
            self.states[state.name] = state
            return

        if error_on_conflict:
            raise ValueError(
                f"State {state.name} already exists "
                f"for phase {state.phase} and overwrite was False"
            )

    def names(self):
        return list(self.states.keys())

    def dims(self):
        return [state.dim for state in self.states.values()]

    def __contains__(self, name):
        return name in self.states
    
    def unpack(self, y):

        states = {}

        start = 0

        for state in self.states.values():

            end = start + state.dim

            value = y[start:end]

            if state.dim == 1:
                value = value[0]

            states[state.name] = value

            start = end

        return states
    def pack(self, state_dict):

        values = []

        for state in self.states.values():

            value = np.asarray(
                state_dict[state.name]
            ).flatten()

            values.extend(value)

        return np.asarray(values)
    
        

@dataclass
class PhaseStateVariable:
    phase: PhaseRef
    state: StateVariable

@dataclass
class PhaseStateCollection:
    phasestates: dict[PhaseRef, StateCollection] = field(default_factory=dict)

    def add(self, phase: PhaseRef, state: StateVariable):
        if phase not in self.phasestates:
            self.phasestates[phase] = StateCollection()

        self.phasestates[phase].add(state)

    def __getitem__(self, phase):
        return self.phasestates[phase]
    def __iter__(self):
        for phase, collection in self.phasestates.items():
            for state in collection.states.values():
                yield PhaseStateVariable(phase, state)

class MultiPhaseVessel():
    def __init__(self,target_comp, temp_ref,
     isothermal, reset_states, controls, h_conv, 
      state_events,population_balance_method,scale,
      adiabatic,jac_type,
      basis='mass_j',ht_mode="jacket",diam=0,area_base=0):
        """ Construct a Reactive Crystallizer Object

    Parameters
    ----------
    mask_params_rxn : list of bool (optional, default = None)
        Binary list of which parameters to exclude from the kinetics
        computations for the reactions
    mask_params_cryst : list of bool (optional, default = None)
        Binary list of which parameters to exclude from the kinetics
        computations for the crystallizations
    method : str
        Choice of the numerical method. Options are: 'moments', '1D-FVM'
    target_comp : str, list of strings
        Name of the crystallizing compound(s) from .json file. 
    scale : float
        Scaling factor by which crystal size distribution will be
        multiplied.
    controls : dict of dicts(funcs) (optional, default = None)
        Dictionary with keys representing the state(e.g.'Temp') which is
        controlled and the value indicating the function to use
        while computing the variable. Functions are of the form
        f(time) = state_value
    adiabatic : bool (optional, default=True)
        Boolean value indicating whether the heat transfer of
        the crystallization is considered.
    rad_zero : float (optional, default=TODO)
        TODO Size of the first bin of the CSD discretization [m]
    reset_states : bool (optional, default = False)
        Boolean value indicating whether the states should be reset
        before simulation
    h_conv : float (optional, default = 0)
        Convective h value for heat transfer to determine U
    
    basis : str (optional, default = mass_j)
        TODO Options :'mass_j'
    jac_type : str
        TODO Options: 'AD'
    state_events : lsit of dict(s)
        list of dictionaries, each one containing the specification of a
        state event
    param_wrapper : callable (optional, default = None)
        function with the signature

            param_wrapper(states, sens)

        Useful when the parameter estimation problem is a function of the
        states y -h(y)- rather than y itself.

        'states' is a DynamicResult object and 'sens' is a dictionary
        that contains N_y number of sensitivity arrays, representing
        time-depending sensitivities. Each array in sens has dimensions
        num_times x num_params. 'param_wrapper' has to return two outputs,
        one array containing h(y) and list of arrays containing
        sens(h(y))
    """

        if isothermal:
            assert adiabatic != 1, "Cannot be isothermal and adiabatic with a reaction present"
            if controls is not None:
                assert 'temp' not in controls.keys(), "Cannot change the temperature of an isothermal unit"

        self.basis = basis
        self.adiabatic = adiabatic
        self.isothermal = isothermal

        self.jac_type = jac_type
        if isinstance(target_comp, str):
            target_comp = [target_comp]
        self.target_comp = target_comp

        self.scale = scale
        self.controls = analyze_controls(controls) #TODO ZZ refactor analyze_controls to give a Controls dataclass, an empty one if controls None
        self.population_balance_method = population_balance_method # TODO ZZ should be either MOM or FVM objects, but children should handle string conversion to object
        self.sensit = None # TODO ZZ may not need a default since plotting removed
        self.jac_states_vals = None # TODO ZZ may not need a default since plotting removed
        self.oper_mode = None
        
        
        # Phase init
        self._Phases = None
        self.Slurry = None
        self._phase_connections = []


        # Parameters for optimization
        self.params_iter = None

        # Crystallization
        self._CrystKinetics = None
        self.material_from_upstream = False
        
        #heat transfer
        self.area_ht = None
        self._Utility = None
        self.ht_mode = ht_mode
        self.h_conv = h_conv
        self.diam = diam
        self.area_base = area_base

        #Reaction
        self.temp_ref = temp_ref
        self._RxnKinetics = None
        self._reaction_regions = []
        
        #State events
        if state_events is None:
            state_events = []
        self.state_event_list = state_events

        #state initialization, all types
        self._initialize_states(reset_states)
        self._Inlet = None
        # self._balances = []


    @property
    def Phases(self):
        return self._Phases
    
    @Phases.setter
    def Phases(self, phases):#TODO turn into mixedphase
        self._normalize_phases(phases)
        self.__original_phase_dict__ = [copy.deepcopy(phase.__dict__) for phase in self._Phases]
        self._post_set_phases()
        
    def _normalize_phases(self,phases):
        if isinstance(phases, (list, tuple)): 
            self._Phases = phases
        elif isinstance(phases, Slurry):
            self._Phases = phases.Phases
        elif phases.__module__ == 'PharmaPy.Phases':
            if self._Phases is None:
                self._Phases = [phases]
            else:
                self._Phases.append(phases)
        else:
            raise RuntimeError('Please provide a list or tuple of phases '
                               'objects')

    def _initialize_slurry(self):
        #deprecated
        if isinstance(self._Phases, Slurry):
            self.Slurry = self._Phases
        elif isinstance(self._Phases, (list, tuple)):
            if len(self._Phases) > 1:
                # Mixed phase
                self.Slurry = Slurry()
                self.Slurry.Phases = self._Phases
        
        self.vol_slurry = copy.copy(self.Slurry.vol)
        if isinstance(self.vol_slurry, np.ndarray):
            self.vol_phase = self.vol_slurry[0]
        else:
            self.vol_phase = self.vol_slurry

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
            dim=self.num_species,
            units=self._basis_units(),
            index=self.name_species,
            
        )

    def _post_set_phases(self):
        
        # Names and target compounds for crystallization
        self.name_species = self.mother_liquor.name_species
        self.num_species = len(self.name_species)
        if self.target_comp is not None:
            self.target_ind = []
            for tc in self.target_comp:
                name_bool = [name == tc for name in self.name_species] #TODO check that it selects correctly
                self.target_ind.append(np.where(name_bool)[0][0])
        # # Input defaults TODO delete this if not necessary
        # self.input_defaults = {
        #     'distrib': np.zeros_like(self.Solid_1.distrib)}
        # Species
        self.define_material_states()
        self.default_diff_states_from_phases()
        self.nomenclature() 

    def _initialize_state_collections(self):

        self.phase_states = PhaseStateCollection()

        # Stream variables are not phases. They describe inlet/outlet connections.
        self.stream_states = StateCollection()

        # States exposed to solver
        self.solver_state_collection = StateCollection()

        # States exposed as outputs
        self.output_state_collection = StateCollection()

    @property
    def Inlet(self):
        return self._Inlet


    @Inlet.setter
    def Inlet(self, inlet):

        self._Inlet = inlet if isinstance(inlet,list) else [inlet]

        self._create_default_inlet_connections()
    
    @property
    def inlet_connections(self):
        return self._inlet_connections


    @inlet_connections.setter
    def inlet_connections(self, connections):

        if not isinstance(connections, list):
            raise TypeError(
                "inlet_connections should be a list"
            )

        if not all(
            isinstance(c, InletConnection)
            for c in connections
        ):
            raise TypeError(
                "inlet_connections must contain "
                "InletConnection objects"
            )

        self._inlet_connections = connections
    def _create_default_inlet_connections(self):

        self._inlet_connections = []

        inlet = getattr(self, "_Inlet", None)

        if inlet is None:
            return

        inlet_phases = (
            inlet.Phases
            if hasattr(inlet, "Phases")
            else [inlet]
        )

        type_counter = {}

        for phase in inlet_phases:

            phase_type = (
                phase.__class__.__name__
                .replace("Phase", "")
                .lower()
            )

            idx = type_counter.get(
                phase_type,
                0
            )

            type_counter[phase_type] = idx + 1

            source_ref = PhaseRef(
                phase_type,
                idx
            )

            sink_ref = PhaseRef(
                phase_type,
                idx
            )

            self._inlet_connections.append(
                InletConnection(
                    source_phase=source_ref,
                    sink_phase=sink_ref
                )
            )
    def define_material_states(self):


        material_state = self._material_state_definition()

        for i, phase in enumerate(self.Phases):

            phase_ref = PhaseRef(
                phase_type=phase.__class__.__name__.replace("Phase","").lower(),
                index=i
            )

            self.phase_states.add(
                phase_ref,
                copy.deepcopy(material_state)
            )
    def default_diff_states_from_phases(self):
        """These are which states from the base phases should be in the solver
        Either override this function for different defaults for children, or modify phase_states directly after phase defintion for complex behavior"""
        for phase in self.phase_states:
            if phase.phase==PhaseRef('liquid',0) and phase.state.name==self.basis:
                phase.state.update_variable('state_type','diff')


    def define_stream_states(self):

        self.stream_states.add(
            StateVariable(
                name="vol_flow",
                dim=1,
                units="m3/s",
                stream="inlet"
            )
        )

        self.stream_states.add(
            StateVariable(
                name="temp",
                dim=1,
                units="K",
                stream="inlet"
            )
        )

        # self.stream_states.add(
        #     StateVariable(
        #         name=self.basis,
        #         dim=self.num_species,
        #         units=self._material_state_definition().units,
        #         stream="inlet"
        #     )
        # )
        self.stream_states.add(
            StateVariable(
                name='mole_conc',
                dim=self.num_species,
                units='kmol/m3',
                stream='inlet'
            )
        )
    @property
    def phases_by_type(self):
        out = {
            "liquid": [],
            "solid": [],
            "vapor": []
        }

        for phase in self._Phases:
            if isinstance(phase, LiquidPhase):
                out["liquid"].append(phase)

            elif isinstance(phase, SolidPhase):
                out["solid"].append(phase)

            elif isinstance(phase, VaporPhase):
                out["vapor"].append(phase)

        return out
    @property
    def Liquids(self):
        return self.phases_by_type["liquid"]

    @property
    def Solids(self):
        return self.phases_by_type["solid"]

    @property
    def Vapors(self):
        return self.phases_by_type["vapor"]
    @property
    def mother_liquor(self):
        return self.Liquids[0]
    
    @property
    def CrystKinetics(self):
        raise AttributeError(
            "CrystKinetics is a convenience initializer only. Modify self.phase_connections directly instead."
            "Use phase_connections[ind].kinetics if the kinetics are desired."
        )

    @CrystKinetics.setter
    def CrystKinetics(self, instance: pk.CrystKinetics|list):
        ''' CrystKinetics is only a convenience initializer ONLY that assumes transfer between liquid1 and solid1 for each Crystkinetic in the list
            and that those connections are always active. It also assumes that target_comp is the only species moving for that index
            THIS CANNOT BE SET BEFORE PHASES
        If any more complex behavior is needed, or if phases need to be set after this, the user should set phase_connections directly using a list of PhaseConnection objects'''
        assert self._Phases is not None, 'Phases must be set before using the crystkinetics convenienve API. Otherwise, you must set phase_connections directly'
        if not isinstance(instance,list) and not isinstance(instance,pk.CrystKinetics):
            raise TypeError("CrystKinetics must be set by either a CrystKinetics object or a list of CrystKinetics objects")
        if not isinstance(instance,list):
            instance = [instance]
        self._CrystKinetics = instance
        self._create_default_phase_connections()
        self._post_CrystKinetics_setter()

    def _post_CrystKinetics_setter(self):
        "Place holder in case future children need special behavior"
        pass
        

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

    def _create_default_phase_connections(self):
        ''' Assumes everything occurs between the first liquid and the first solid phases
        Assumes that whatever the target index is for that Crystkinetic, the massfrac is 100% that compound and 0 everything else
        Only runs if phase_connections are not already set'''

        if len(self.phase_connections)>0:
            raise RuntimeError(
                "phase_connections already defined. "
                "Cannot use CrystKinetics convenience API."
            )
        connections = []
        for i,ck in enumerate(self._CrystKinetics):
            weights = np.zeros(self.num_species)
            weights[self.target_ind[i]] = 1
            if ck.supports(['growth','nucl_prim','nucl_sec']):
                # liquid to solid because crystallization is valid
                connection = PhaseConnection(source_phase=PhaseRef("liquid",0),
                                             sink_phase=PhaseRef('solid',0),
                                             kinetics=ck,
                                             species_weights=weights,
                                             active_condition=lambda source,sink:True,
                                             mechanism=self.Solid_1.transfer_mechanism
                                             )
                connections.append(connection)
            if ck.supports('dissolution'):
                #solid to liquid because dissolution
                connection = PhaseConnection(source_phase=PhaseRef("solid",0),
                                             sink_phase=PhaseRef('liquid',0),
                                             kinetics=ck,
                                             species_weights=weights,
                                             active_condition=lambda source,sink:True,
                                             mechanism=DirectTransfer()

                                             )
                connections.append(connection)
        

    @property
    def reaction_regions(self):
        return self._reaction_regions
    @reaction_regions.setter
    def reaction_regions(self, regions):

        if not isinstance(regions, list):
            raise TypeError(
                "reaction_regions is expected to be a list"
            )

        if not all(
            isinstance(r, ReactionRegion)
            for r in regions
        ):
            raise TypeError(
                "reaction_regions should all be "
                "ReactionRegion objects"
            )

        self._reaction_regions = regions

        self.nomenclature(overwrite=True)
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
        self._create_ReactionRegions_from_RxnKinetics(instance)
        self._post_RxnKinetics_setter(instance)

    def _create_ReactionRegions_from_RxnKinetics(self,RKs):
        '''Creates reaction_regions from RxnKinetics, assuming everything is in the liquid phase corresponding to the index of the RxnKinetic
        This is because RxnKinetics handles multiple reactions, but assumes a single phase
        If reaction_regions is already defined, this step is ignored
        For reactions in other phases, the user must define reaction_regions directly
        '''
        if len(self.reaction_regions) >0:
            raise RuntimeError(
                "reaction_regions already defined. "
                "Cannot use RxnKinetics convenience API."
            )
        reaction_regions = []
        for i,rk in enumerate(RKs):
            assert len(self.Liquids)>=i-1, "The number of reaction kinetics must match or be less than the number of liquid phases or you must specify reaction_regions manually"
            region = ReactionRegion(phase=PhaseRef("liquid",i),
                                    kinetics=rk)
            reaction_regions.append(region)
        self.reaction_regions = reaction_regions

    def _post_RxnKinetics_setter(self,RKs):
        "Place holder in case future children need special behavior"
        pass

    @property
    def Utility(self):
        return self._Utility

    @Utility.setter
    def Utility(self, utility):
        self.u_ht = 1 / (1 / self.h_conv + 1 / utility.h_conv)
        self._Utility = utility
    
    
    def complete_state(self, state, time):

        completed = state.copy()

        required_states = (
            list(self.solver_state_collection.states.values())
            + list(self.output_state_collection.states.values())
        )

        for variable in required_states:

            name = variable.name

            if name in completed:
                continue

            if name in self.controls:
                control = self.controls[name]

                completed[name] = control["fun"](
                    time,
                    *control["args"],
                    **control["kwargs"]
                )

                continue

            default = getattr(self, name, None)

            if default is not None:
                completed[name] = default

        return completed
    def __getattr__(self, name):
        if name.startswith("Liquid_"):
            idx = int(name.split("_")[1]) - 1
            return self.Liquids[idx]

        if name.startswith("Solid_"):
            idx = int(name.split("_")[1]) - 1
            return self.Solids[idx]

        if name.startswith("Vapor_"):
            idx = int(name.split("_")[1]) - 1
            return self.Vapors[idx]

        raise AttributeError(name) #TODO Check if this raises unwanted errors
       
    def _initialize_states(self,reset=False):

        self.reset_states = reset

        self.state_variables = StateCollection()

        self.input_states = StateCollection()

        self.output_states = StateCollection()

        self.elapsed_time = 0

        self.outputs = None
        self._initialize_state_collections()

        
        #TODO add what's needed to solver_states (temp etc.) based on nomenclature logic here
        
    def reset(self):
        self._initialize_states(self.reset_states)
        self.nomenclature()

        for phase, di in zip(self.Phases, self.__original_phase_dict__):
            phase.__dict__.update(di)

        self.profiles_runs = []
        # TODO does not work if more than 1 liquid and 1 solid
    
    def _eval_state_events(self, time, states, sw):
        # TODO reactor version changes discretized_model to True if PFR (cobc in our case)
        events = eval_state_events(
            time, states, sw, self.len_states,
            self.solver_states, self.state_event_list, sdot=self.derivatives,
            discretized_model=False)

        return events
    
    def heat_transfer(self, temp, temp_ht, vol):
        # Heat transfer area ##r
        heat_transf = self.u_ht * self.area_ht * (temp - temp_ht)
        return heat_transf
    
    def define_solver_states(self,overwrite=False):

        
        for phase_state in self.phase_states:
            if phase_state.state.state_type!='diff': continue
            state_copy = copy.deepcopy(
                phase_state.state
            )
            state_copy.phase = (
                phase_state.phase
            )
            self.solver_state_collection.add(
               state_copy,overwrite
            )
        
        if self.adiabatic:

            self.solver_state_collection.add(
                StateVariable(
                    name="temp",
                    dim=1,
                    units="K",
                    state_type='diff'
                ),overwrite
            )
            
        elif "temp" not in self.controls:

            self.solver_state_collection.add(
                StateVariable(
                    name="temp",
                    dim=1,
                    units="K",
                    state_type='diff'
                ),overwrite
            )
            self.solver_state_collection.add(
                StateVariable(
                    name="temp_ht",
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

    def define_output_states(self,overwrite=False):

        self.output_state_collection.add(
            StateVariable(
                name="q_rxn",
                dim=1,
                units="W",
                state_type="alg",
                compute_history=self.compute_qrxn_history
            ),overwrite
        )

        self.output_state_collection.add(
            StateVariable(
                name="q_ht",
                dim=1,
                units="W",
                state_type="alg",
                compute_history=self.compute_qht_history
            ),overwrite
        )

        self.output_state_collection.add(
            StateVariable(
                name="m_flow",
                dim=1,
                units="kg/s",
                state_type="alg",
                compute_history=self.compute_outlet_massflow_history
            ),overwrite
        )

        for region in self.reaction_regions:

            region.kinetics.add_output_state_variables(
                self.output_state_collection,
                reaction_region=region,
                overwrite=overwrite
            )

        

        for conn in self.phase_connections:

            if conn.mechanism is not None:
                conn.mechanism.add_output_state_variables(
                    self.output_state_collection,overwrite
                )

    def nomenclature(self,overwrite=False):

        self.define_solver_states(overwrite)
        self.define_output_states(overwrite)
        self.define_stream_states()
        self.build_states_in_dict()
        self.name_states = self.solver_state_collection.names()
        self.dim_states = self.solver_state_collection.dims()

    def _evaluate_stream(self, stream_name, stream, time):

        inputs = {}

        stream_states = [
            state for state in self.stream_states.states.values()
            if state.stream == stream_name
        ]

        for state in stream_states:

            value = getattr(stream, state.name)

            if callable(value):
                value = value(time)

            inputs[state.name] = np.asarray(value)

        return inputs
    
    def get_current_inlet(self, time):

        inlet = copy.deepcopy(self.Inlet)

        for state in self.stream_states.states.values():

            if not hasattr(inlet, state.name):
                continue

            value = getattr(inlet, state.name)

            if callable(value):
                value = value(time)

            setattr(inlet, state.name, value)

        return inlet
    def get_phase(self,phase_ref:PhaseRef):
        return get_phase_from_ref(self.Phases,phase_ref)
    def build_states_in_dict(self):
        #exists for backwards compatiblity
        self.states_in_dict = {}

        for phase_state in self.phase_states:

            phase_name = self.phase_ref_to_legacy_name(
                phase_state.phase
            )

            self.states_in_dict.setdefault(
                phase_name,
                {}
            )

            self.states_in_dict[phase_name][
                phase_state.state.name
            ] = phase_state.state.dim

        for state in self.stream_states.states.values():

            stream_name = (
                state.stream.capitalize()
                if state.stream is not None
                else "Inlet"
            )

            self.states_in_dict.setdefault(
                stream_name,
                {}
            )

            self.states_in_dict[stream_name][
                state.name
            ] = state.dim
    def update_phases_from_state(self, completed_state):

        for phase_ref, collection in self.phase_states.phasestates.items():

            phase = self.get_phase(phase_ref)

            updates = {}

            for variable in collection.states.values():

                if variable.name in completed_state:
                    updates[variable.name] = completed_state[variable.name]

            if updates:
                phase.updatePhase(**updates)


        for connection in self.phase_connections:

            mechanism = connection.mechanism

            if mechanism is None:
                continue

            mechanism.update_state(
                completed_state
            )
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

        for state in self.solver_state_collection.states.values():

            if state.state_type != "diff":
                continue

            if state.phase is not None:

                if state.phase not in material_rates:
                    continue

                packed.extend(np.asarray(material_rates[state.phase]).flatten())

            else:

                if state.name not in global_rates:
                    continue

                packed.extend(np.asarray(global_rates[state.name]).flatten())

        return np.asarray(packed)
    
    def unit_model(self, time, states, params=None, sw=None,
                    mat_bce=False, enrgy_bce=False):
        unpacked_state = self.solver_state_collection.unpack(states)
        completed_state = self.complete_state(unpacked_state,time)
        self.update_phases_from_state(completed_state)
        
        u_input = self.get_current_inlet(time)

        # Balances
        material_rates, material_contributions = self.material_balances(
            time,completed_state, params, u_input)
        

        if mat_bce:
            return self.pack_state_rates(material_rates)
        global_rates = {}
        energy_rates = self.energy_balances(
                time,completed_state, params, u_input, material_contributions)
        global_rates.update(energy_rates)
        if "temp_ht" in self.solver_state_collection:
            utility_rates = self.utility_energy_balance(time,completed_state,params,u_input)

        global_rates.update(utility_rates)
        if enrgy_bce:
            return self.pack_state_rates(global_rates=global_rates)

        balances = self.pack_state_rates(material_rates=material_rates,
                                         global_rates=global_rates)

        self.derivatives = balances
        return balances

    

    def unit_jacobian(self, t, y):
        return self.jac_states_fun(t, y)

    def jac_states_numerical(self, time, states, params, return_only=True):
        #TODO check if necessary
        if return_only:
            return self.jac_states_vals
        else:
            def wrap_states(st): return self.unit_model(time, st, params)

            abstol = self.sundials_opt['atol']
            reltol = self.sundials_opt['rtol']
            jac_states = numerical_jac_central(wrap_states, states,
                                               dx=dx_jac_x,
                                               abs_tol=abstol, rel_tol=reltol)

            return jac_states

    def jac_params_numerical(self, time, states, params):
        #TODO check if necessary
        def wrap_params(theta): return self.unit_model(time, states, theta)

        abstol = self.sundials_opt['atol']
        reltol = self.sundials_opt['rtol']
        p_bar = self.sundials_opt['pbar']

        dp = np.abs(p_bar) * np.sqrt(max(reltol, eps))

        jac_params = numerical_jac_central(wrap_params, params,
                                           dx=dp,
                                           abs_tol=abstol, rel_tol=reltol)

        return jac_params
    
    def rhs_sensitivity(self, time, states, sens, params):

        jac_params_vals = self.jac_params_fn(time, states, params)

        jac_states_vals = self.jac_states_fn(time, states, params,
                                             return_only=False)

        rhs_sens = np.dot(jac_states_vals, sens) + jac_params_vals

        self.jac_states_vals = jac_states_vals

        return rhs_sens
    
    def compile_integrator(self,eval_sens=False,
                           jac_v_prod=False,
                           sundials_opts=None,
                           verbose=True,any_event=True, return_flag=False):
        '''Builds the unit and compiles the integrator without calling it'''
        # chatGPT helped with this refactor 20260226
        self._compiled=True
        self._compiled_val_sens= eval_sens
        self._compiled_jac_v_prod = jac_v_prod
        self.reset()
        states_init = self.solver_state_collection.pack(self.complete_state({}, 0))#TODO check that it is correct
        merged_params = None
        problem = self.set_ode_problem(
            eval_sens,
            states_init,
            merged_params,
            jac_v_prod
        )
        if len(self.state_event_list) > 0:
            def new_handle(solver, info):
                return handle_events(solver, info, self.state_event_list,
                                     any_event=any_event)

            problem.state_events = self._eval_state_events
            problem.handle_event = new_handle

        solver = CVode(problem)
        solver.iter = 'Newton'
        solver.discr = 'BDF'

        if sundials_opts:
            for name, val in sundials_opts.items():
                setattr(solver, name, val)
                if name == 'time_limit':
                    solver.report_continuously = True

        if not verbose:
            solver.verbosity = 50
        # solver.verbosity = 0 #todo see why the other one did not work

        
        if eval_sens:
            solver.sensmethod = 'SIMULTANEOUS'
            solver.suppress_sens = False
            solver.report_continuously = True

        if self.population_balance_method == '1D-FVM':
            solver.linear_solver = 'SPGMR'  # large, sparse systems

        if not verbose:
            solver.verbosity = 50
        

        self._problem = problem
        self._solver = solver
        if return_flag:
            return states_init, merged_params

    def set_ode_problem(self, eval_sens, states_init, params_mergd,
                        jacv_prod):
        #eval_sens is populated by sim_exec and so must be accepted in order, but is no longer useful to this model
        if self.state_event_list is None:
            def model(time, states, params=params_mergd):
                return self.unit_model(time, states, params)

            problem = Explicit_Problem(model, states_init,
                                        t0=self.elapsed_time)
        else:
            sw0 = [True] * len(self.state_event_list) #switches, currently unused in unit_model
            def model(time, states, sw=None):#equivalent to fobj in reactor
                return self.unit_model(time, states, params_mergd, sw)

            problem = Explicit_Problem(model, states_init,
                                        t0=self.elapsed_time, sw0=sw0)

        # ----- Jacobian callables
        if self.population_balance_method == 'moments':
            # w.r.t. states
            # problem.jac = lambda time, states: \
            #     self.unit_jacobians(time, states, None, params_mergd,
            #                         None, None)

            pass

        elif self.population_balance_method == 'fvm':
            # J*v product (AD, slower than the one used by SUNDIALS)
            if jacv_prod:
                problem.jacv = lambda time, states, fy, v: \
                    self.unit_jacobians(time, states, None, params_mergd,
                                        fy, v)

        return problem
    
    def _build_initial_state_and_params(self):
        self.reset()
        return self.solver_state_collection.pack(
            self.complete_state({}, 0)
        )
    
    def _fast_solve(self, runtime, time_grid,
                  verbose=True, test=False,
                   sundials_opts=None, any_event=True):

        
        self.reset()
        states_init= self.solver_state_collection.pack(self.complete_state({}, 0))

        # update time
        if runtime is not None:
            final_time = runtime + self.elapsed_time
        if time_grid is not None:
            final_time = time_grid[-1]

        # reinitialize
        self._solver.t = self.elapsed_time if time_grid is None else time_grid[0]
        self._solver.y = states_init
        self._solver.initialize()
        

        time, states = self._solver.simulate(final_time, ncp_list=time_grid)

        # DO NOT call retrieve_results()

        return time, states
    
    def solve_unit(self, runtime=None, time_grid=None,
                   eval_sens=False,
                   jac_v_prod=False, verbose=True, test=False,
                   sundials_opts=None, any_event=True):
        """
        runtime : float (default = None)
            Value for the total unit runtime
        time_grid : list of float (optional, dafault = None)
            Optional list of time values for the integrator to use
            during simulation
        eval_sens : bool (optional, default = False)
            Boolean value indicating whether the parametric
            sensitivity system will be included during simulation.
            Must be True to access sensitivity information.
        jac_v_prod :
            TODO
        verbose : bool (optional, default = True)
            Boolean value indicating whether the simulator will
            output run statistics after simulation is complete.
            Use True if you want to see the number of function
            evaluations and wall-clock runtime for the unit.
        test :
            TODO
        sundials_opts :
            TODO
        any_event :
            TODO
        """
        if hasattr(self, "_compiled") and self._compiled:
            if eval_sens != self._compiled_val_sens:
                raise RuntimeError("Must recompile integrator when eval_sens changes")
            if jac_v_prod != self._compiled_jac_v_prod:
                raise RuntimeError("Must recompile integrator when jac_v_prod changes")
            return self._fast_solve(runtime, time_grid,
                                    verbose, test,
                                    sundials_opts, any_event)
        # merged_params = np.append(merged_params,self.RxnKinetics.concat_params()[self.mask_params_rxn])
        # states_init = np.append(states_init,self.Liquid_1.mole_conc)
        # ---------- Create problem
        states_init,merged_params = self.compile_integrator(eval_sens=eval_sens,
                                jac_v_prod=jac_v_prod,
                                sundials_opts=sundials_opts,
                                verbose=verbose,any_event=any_event,return_flag=True)

        self.derivatives = self._problem.rhs(self.elapsed_time, states_init,
                                       merged_params)

        self.sundials_opt = self._solver.get_options()
        if runtime is not None:
            final_time = runtime + self.elapsed_time
        if time_grid is not None:
            final_time = time_grid[-1]
       
        time, states= self._solver.simulate(final_time, ncp_list=time_grid)
        
        self.retrieve_results(time, states)

        return time, states

    def flatten_states(self):
        out = flatten_states(self.profiles_runs)

        return out
    def initialize_phase_rate_dictionary(self):
        material_rates = {}

        for phase_ref, phase_states in self.phase_states.phasestates.items():

            material_state = next(
                s for s in phase_states.states.values()
                if s.name == self.basis
            )

            material_rates[phase_ref] = np.zeros(material_state.dim)
        return material_rates
    
    def material_balances(self,time,completed_state, params, u_input):

        
        contributions, aux = self.limit_material_rates(
            time,
            completed_state,
            params,
            u_input
        )

        rates = self.sum_material_contributions(
            contributions
        )

        return rates, aux

    def check_negative_inventory(
            self,
            rates,
            completed_state):

        violations = {}

        for phase_ref, rate in rates.items():

            phase = self.get_phase(
                phase_ref
            )

            inventory = getattr(
                phase,
                self.basis
            )

            mask = (inventory + rate < -eps)

            if np.any(mask):

                violations[phase_ref] = mask

        return violations
    def calculate_scale(
            self,
            violations,
            contributions):

        scales = {}

        scalable_terms = (
            "reaction",
            "transfer"
        )

        for phase_ref, mask in violations.items():

            phase = self.get_phase(
                phase_ref
            )

            inventory = getattr(
                phase,
                self.basis
            )

            total_consumption = np.zeros_like(
                inventory
            )

            for term in scalable_terms:

                total_consumption += np.minimum(
                    contributions[term][phase_ref],
                    0
                )

            scale = np.ones_like(
                inventory
            )

            violating_species = np.where(mask)[0]

            for ind in violating_species:

                consumed = -total_consumption[ind]

                if consumed <= 0:
                    continue

                scale[ind] = min(
                    1.0,
                    inventory[ind] / (consumed+eps)
                )

            scales[phase_ref] = scale

        return scales
    def scale_phase_inventory(
            self,
            scales):

        for phase_ref, scale_vector in scales.items():

            phase = self.get_phase(
                phase_ref
            )

            new_mass = (
                phase.mass_j
                * scale_vector
            )

            phase.updatePhase(
                mass_j=new_mass
            )
                
    
    
    def limit_material_rates(
            self,
            time,
            completed_state,
            params,
            u_input
        ):
        """Soft constraint to ensure mass conservation is not violated between the generation and phase transference. If violation is detected, the  
        mass movement is scaled back by scaling the effective mass available until no violations are detected. The original mass is then restored from completed_state"""
        
        for iteration in range(5):
            contributions,aux= self.calculate_material_contributions(time,completed_state, params, u_input)

            rates = self.sum_material_contributions(contributions)

            violations = self.check_negative_inventory(rates,completed_state)

            if len(violations) == 0:
                self.update_phases_from_state(completed_state)
                return contributions, aux
            scales = self.calculate_scale(violations,contributions)

            self.scale_phase_inventory(scales)
        self.update_phases_from_state(completed_state)
        raise RuntimeError(
            "Material balance limiter failed "
            "to converge after 5 iterations."
        )

    def sum_material_contributions(
        self,
        contributions):

        rates = self.initialize_phase_rate_dictionary()

        for contribution in contributions.values():

            for phase_ref, rate in contribution.items():

                rates[phase_ref] += rate

        return rates
    def calculate_material_contributions(self,time,completed_state, params, u_input):

        contributions = {
            "inlet": self.initialize_phase_rate_dictionary(),
            "reaction": self.initialize_phase_rate_dictionary(),
            "transfer": self.initialize_phase_rate_dictionary(),
            "outlet": self.initialize_phase_rate_dictionary()
        }
        aux = {
            "inlet": [],
            "reaction": [],
            "transfer": [],
            "outlet": []
        }

        self.add_inlet_terms(
            contributions["inlet"],
            aux['inlet'],
            time,
            completed_state,
            params,
            u_input
        )

        self.add_reaction_terms(
            contributions["reaction"],
            aux['reaction'],
            time,
            completed_state,
            params,
            u_input
        )

        self.add_transfer_terms(
            contributions["transfer"],
            aux['transfer'],
            time,
            completed_state,
            params,
            u_input
        )

        self.add_outlet_terms(
            contributions["outlet"],
            aux['outlet'],
            time,
            completed_state,
            params,
            u_input
        )

        return contributions, aux
    
    def add_inlet_terms(
            self,
            rates:dict,
            aux:list,
            time,
            completed_state,
            params,
            u_input
        ):

        if u_input is None:
            return

        inlet_phases = (
            u_input.Phases
            if hasattr(u_input, "Phases")
            else [u_input]
        )

        for connection in self.inlet_connections:
            inlet_phase =get_phase_from_ref(
                inlet_phases,
                connection.source_phase
            )
            rates[connection.sink_phase] += inlet_phase.material_flow(self.basis)
            aux.append({
                "connection": connection,
                "phase": inlet_phase,
                "flow": inlet_phase.material_flow(self.basis),
                'note':'flow is in basis units'
                })

    def add_outlet_terms(
            self,
            rates,
            aux:list,
            time,
            completed_state,
            params,
            u_input
        ):

        if not hasattr(self,"Outlet") or self.Outlet is None:
            return

        for phase_ref in rates:

            phase = self.get_phase(phase_ref)

            rates[phase_ref] -= (phase.material_flow(self.basis))
            aux.append({
                "phase_ref": phase_ref,
                "phase": phase,
                "flow": phase.material_flow(self.basis),
                'note':'flow is in basis units'
                })

    def add_reaction_terms(
            self,
            rates,
            aux,
            time,
            completed_state,
            params,
            u_input,
            molarity_molPerLiter=True
        ):

        # temp = completed_state["temp"]
        mole_adjust = 1000 if molarity_molPerLiter else 1

        for region in self.reaction_regions:

            phase = self.get_phase(
                region.phase
            )
            temp = phase.temp

            rk = region.kinetics

            mask = np.array([
                species in rk.partic_species
                for species in self.name_species
            ])

            conc = phase.mole_conc
            if rk.keq_params is not None:

                deltah_rxn = phase.getHeatOfRxn((
                        rk.stoich_matrix,
                        temp,
                        mask,
                        rk.delta_hrxn,
                        rk.tref_hrxn
                    ))

            else:

                deltah_rxn = None
            #if deltaH_rxn is NOne, then the get_rxn_rates does not evaluate Keq. 
            # Because that determines how much material changes, 
            # this is a material question rather than an energy balance question 
            # even though deltaH is normally an energy balance concern
            reaction_rates = rk.get_rxn_rates(
                conc,
                temp,
                overall_rates=False,
                delta_hrxn=deltah_rxn
            )
            

            species_rates = rk.get_rxn_rates(
                conc[mask],
                temp,
                overall_rates=True
            )

            species_massPerVol_rates = np.zeros(self.num_species)
            species_massPerVol_rates[mask] = species_rates
            species_massPerVol_rates *= phase.mw
            species_mass_rates = species_massPerVol_rates* phase.vol*mole_adjust

            aux.append({
                "region": region,
                "phase": phase,
                "rxn_rates": reaction_rates,
                'species_molar_rates':species_rates,
                "species_massPerVol_rates":species_massPerVol_rates,
                "species_mass_rates":species_mass_rates
                }) 

            rates[region.phase] += species_mass_rates

    def add_transfer_terms(
            self,
            rates,
            aux:list,
            time,
            completed_state,
            params,
            u_input
        ):

        # temp = completed_state["temp"]

        for connection in self.phase_connections:

            source = self.get_phase(connection.source_phase)

            sink = self.get_phase(connection.sink_phase)

            if not connection.active_condition(source,sink):
                continue

            transfer_rate = connection.mechanism.get_transfer_rate(
                source_phase=source,
                sink_phase=sink,
                connection=connection,
                completed_state=completed_state,
                time=time,
                params=params
            )

            if connection.species_weights is None:

                species_rate = (transfer_rate* source.mass_frac)

            else:

                species_rate = (transfer_rate* np.asarray(connection.species_weights))

            rates[connection.source_phase] -= species_rate

            rates[connection.sink_phase] += species_rate
            aux.append({
                "connection": connection,
                "source": source,
                "sink": sink,
                "transfer_rate": transfer_rate,
                "species_rate": species_rate,
                'note':'rates are in mass'
            })
    
    def energy_balances(
            self,
            time,
            completed_state,
            params,
            u_input,
            aux
        ):

        contributions = {
            "inlet": 0,
            "reaction": 0,
            "transfer": 0,
            "outlet": 0,
            "utility": 0
        }

        self.add_inlet_energy_terms(
            contributions,
            aux["inlet"],
            time,
            completed_state,
            params,
            u_input
        )

        self.add_reaction_energy_terms(
            contributions,
            aux["reaction"],
            time,
            completed_state,
            params,
            u_input
        )

        self.add_transfer_energy_terms(
            contributions,
            aux["transfer"],
            time,
            completed_state,
            params,
            u_input
        )

        self.add_outlet_energy_terms(
            contributions,
            aux["outlet"],
            time,
            completed_state,
            params,
            u_input
        )

        self.add_utility_energy_terms(
            contributions,
            time,
            completed_state,
            params,
            u_input
        )

        qdot = sum(
            contributions.values()
        )

        heat_capacity = (
            self.Phases.getCP()
        )
        dtemp_dt = qdot / heat_capacity
        return {"temp": dtemp_dt}

    def add_inlet_energy_terms(
            self,
            contributions,
            aux,
            time,
            completed_state,
            params,
            u_input
        ):

        for inlet in aux:

            phase = inlet["phase"]

            h_in = phase.getEnthalpy(
                phase.temp,
                temp_ref=self.temp_ref,
                total_h=True,
                basis='mass'
            )  # J/kg mixture

            contributions["inlet"] += (
                inlet["flow"] * h_in
            )
    def add_transfer_energy_terms(
            contributions,
            aux,
            time,
            completed_state,
            params,
            u_input
        ):
        
        for transfer in aux:
            connection = transfer['connection']
            contributions['transfer'] += connection.mechanism.get_heat_generation(transfer_rate=transfer['transfer_rate'],
                                                             source=transfer['source'],
                                                             sink=transfer['sink'],
                                                             connection=connection,
                                                             completed_state=completed_state,
                                                             time=time)
    def add_reaction_energy_terms(self,
                                  contributions,
                                  aux,
                                  time,
                                  completed_state,
                                  params,
                                  u_input,
                                  molarity_molPerLiter=True
                                ):
        mole_adjust = 1000 if molarity_molPerLiter else 1
        for rxn in aux:
            phase = rxn['phase']
            temp = phase.temp
            rk = rxn['region'].kinetics
            mask = np.array([
                species in rk.partic_species
                for species in self.name_species
            ])

            deltah_rxn = (
                phase.getHeatOfRxn(
                    rk.stoich_matrix,
                    temp,
                    mask,
                    rk.delta_hrxn,
                    rk.tref_hrxn
                ))
            contributions["reaction"] += -(deltah_rxn * rxn["rxn_rates"]).sum()* phase.vol *mole_adjust # molarity here is mol/L, but J/kg is expected for internal consistency

    def add_outlet_energy_terms(
            self,
            contributions,
            aux,
            time,
            completed_state,
            params,
            u_input
        ):

        for outlet in aux:

            phase = outlet["phase"]
            temp= phase.temp
            h_out = phase.getEnthalpy(
                temp,
                temp_ref=self.temp_ref,
                mass_frac=phase.mass_frac,
                total_h=True,
                basis='mass'
            )

            contributions["outlet"] -= (
                outlet["flow"] * h_out
            )
    def get_heat_transfer_temperature(self):
        """The temperature to use to determine heat transfer from the vessel. The default assumption is to use the temperature of the first liquid phase"""
        return self.get_phase(PhaseRef("liquid", 0)).temp
    
    def add_utility_energy_terms(
            self,
            contributions,
            time,
            completed_state,
            params,
            u_input
        ):

        if self.Utility is None or self.isothermal:
            return

        temp = self.get_heat_transfer_temperature()
        temp_ht = completed_state["temp_ht"]

        qdot = self.get_heat_transfer_rate(
            temp,
            temp_ht
        )

        contributions["utility"] -= qdot
    def get_heat_transfer_rate(
            self,
            temp,
            temp_ht
        ):

        if self.ht_mode == "coil":
            raise NotImplementedError

        area = self.get_heat_transfer_area()

        return (
            self.u_ht
            * area
            * (temp - temp_ht)
        )
    def get_heat_transfer_area(self):

        liquid = self.get_phase(
            PhaseRef("liquid",0)
        )

        return (
            4 * liquid.vol / self.diam
            + self.area_base
        )
    def utility_energy_balance(
            self,
            time,
            completed_state,
            params,
            u_input
        ):

        if (self.Utility is None
            or "temp_ht" not in completed_state):
            return {}

        temp = self.get_heat_transfer_temperature()
        temp_ht = completed_state["temp_ht"]

        ht_controls = self.Utility.get_inputs(time)

        flow_ht = ht_controls["vol_flow"]
        temp_in = ht_controls["temp_in"]

        cp = self.Utility.cp
        rho = self.Utility.rho

        if "vol" in completed_state:
            vol_ht = self.vol_ht
        else:
            liquid = self.get_phase(PhaseRef("liquid",0))
            vol_ht = liquid.vol * 0.15

        qdot = self.get_heat_transfer_rate(temp,temp_ht)
        dtemp_ht_dt = (flow_ht / vol_ht
            * (temp_in - temp_ht)+
            qdot/(rho * cp * vol_ht))
        return {"temp_ht":dtemp_ht_dt}
    def build_outlet_stream(self, time):
        #TODO when multiple outlets needed, do for outlet in self.Outlets
        if self.Outlet is None:
            return

        outlet = copy.deepcopy(self.Outlet)

        outlet.Phases = tuple(
            copy.deepcopy(p)
            for p in self.Phases
        )

        self.Outlet = outlet
    def build_output_history(
            self,
            time,
            solver_history
        ):

        history = {}

        for state in self.output_state_collection.states.values():

            history[state.name] = state.compute_history(
                time,
                solver_history,
                self
            )

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
    def retrieve_results(self, time, solver_states):

        solver_history = self.build_solver_history(time,solver_states)

        output_history = self.build_output_history(time,solver_history)

        self.outputs = DynamicResult(**solver_history,**output_history)

        self.update_final_state(solver_history)

        self.elapsed_time = time[-1]

        return self.outputs

class ReactiveMSMPR(MultiPhaseVessel):
    """
    Assumes:
        constant volume
        constant solid density
        metric units
        
    """
    def __init__(self, target_comp, mask_params_rxn=None,mask_params_cryst=None, temp_ref=298.15, isothermal=True,
                  reset_states=False, controls=None, h_conv=1000, ht_mode='jacket',
                  return_sens=True, state_events=None, method='1D-FVM',
                  scale=1, vol_tank=None, adiabatic=False, rad_zero=0, vol_ht=None,
                  basis='mass_conc', jac_type=None, param_wrapper=None, num_interp_points=3, grid_size=500, param_estimation_run=False):
        super().__init__(target_comp,mask_params_rxn,mask_params_cryst, temp_ref, isothermal, reset_states, controls, h_conv, ht_mode, return_sens, state_events, method, scale, vol_tank, adiabatic, rad_zero, vol_ht, basis, jac_type, param_wrapper)
        self.is_continuous = True
        self.oper_mode = 'Continuous'
        self._Inlet = None
        self.vol_offset = 0.75
        self.num_interp_points = num_interp_points
        self.mydistrib = np.zeros(grid_size)
        # self.checker = progress_checker(flag='Solver')
        self.kin_array = {}
        self.param_estimation_run=param_estimation_run
        self._zero_flow=False
    @property
    def Inlet(self):
        return self._Inlet

    @Inlet.setter
    def Inlet(self, inlet_object):
        self._Inlet = inlet_object
        self._Inlet.num_interpolation_points = self.num_interp_points

    def _get_tau(self):
        time_upstream = getattr(self.Inlet, 'time_upstream', None)
        if time_upstream is None:
            time_upstream = [0]

        inputs = self.get_inputs(time_upstream[-1])

        volflow_in = inputs['Inlet']['vol_flow']
        tau = self.Liquid_1.vol / volflow_in

        self.tau = tau
        return tau
    # TODO solve_steady_state

    def material_balances(self, time, params, u_inputs, rhos, mu_n,
                          distrib, mass_conc, mole_conc, temp, temp_ht, vol, phi_in):

        rho_sol = rhos[0][1]
        # self.checker.check(time)

        input_flow = u_inputs['Inlet']['vol_flow']

        input_conc = u_inputs['Liquid_1']['mass_conc']
        # input_mole = u_inputs['Inlet']['mole_conc']
        ## Reactive terms:
        if self.RxnKinetics.keq_params is None:
            rate = self.RxnKinetics.get_rxn_rates(mole_conc[self.mask_species],temp)
        else:
            deltah_rxn = self.Liquid_1.getHeatOfRxn(temp,
                                                    self.Kinetics.tref_hrxn)

            rate = self.RxnKinetics.get_rxn_rates(mole_conc[self.mask_species],
                                               temp,
                                               deltah_rxn)
        rates = np.zeros_like(mole_conc)
        rates[self.mask_species] = rate
        rxn_term = rates*self.Liquid_1.mw #calc rates as moles convert to mass (mol/t to kg/t) zzz
        new_mass_conc = np.clip(mass_conc + rxn_term, 0, None)
        
        if self.population_balance_method == 'moments':
            input_distrib = u_inputs['Inlet']['mu_n'] * (1e6)**np.arange(self.num_distr)#* self.scale
            ddistr_dt, transf = self.method_of_moments(distrib, new_mass_conc, temp,
                                                       params, rho_sol)
        elif self.population_balance_method == '1D-FVM':
            input_distrib = u_inputs['Inlet']['distrib'] * self.scale
            if True:#time < 5200:
                ddistr_dt, transf = self.fvm_method(distrib, mu_n, new_mass_conc, temp,
                                                    params, rho_sol)
                nuclp,sec, growth, dissol = self.CrystKinetics.get_kinetics(new_mass_conc, temp,
                                                            self.Solid_1.kv, mu_n,nucl_sec_out=True)
                # self.oldparams = nuclp,sec,growth,dissol,ddistr_dt,transf
            else:
                nuclp,sec,growth,dissol,ddistr_dt,transf = self.oldparams
            self.kin_array[time] = [nuclp,sec,growth,float(transf)]
            


            self.Solid_1.moments[[2, 3]] = mu_n[[2, 3]]
        
        # ---------- Add flow terms
        # Distribution
        tau_inv = input_flow / vol #theta in many nomenclatures
        flow_distrib = tau_inv * (input_distrib - distrib)

        ddistr_dt = ddistr_dt + flow_distrib
        # Liquid phase
        phi = 1 - self.Solid_1.kv * mu_n[3] #epsilon in documentation

        c_tank = new_mass_conc
        # Re derive MSMPR to match basis and add reaction here
        #check how handle multiple species (if not array follows)
        # check how incorporate stoichs
        # check if need stoichs in mom
        # check if handles multiple cryst species # needed
        
        flow_term = tau_inv * (input_conc*phi_in[0] - c_tank*phi) #check phi_in[0] or just phi_in
        transf_term = transf * (self.kron_jtg - c_tank / rho_sol)
        # check if units right
        dcomp_dt = 1 / phi * (flow_term - transf_term + rxn_term)

        if self.basis == 'mass_frac':
            rho_liq = self.Liquid_1.getDensity()
            dcomp_dt *= 1 / rho_liq

        dmaterial_dt = np.concatenate((ddistr_dt, dcomp_dt))
        self.mydistrib = np.append(self.mydistrib,self.mydistrib[-1]+ddistr_dt).reshape(-1,ddistr_dt.shape[0])
        # self.mass_flow_out = self.Solid_1.mass_flow
        return dmaterial_dt, transf
    def energy_balances(self, time,params, cryst_rate, u_inputs, rhos, mu_n,
                        distrib, mass_conc,mole_conc, temp, temp_ht, vol,
                        h_in, heat_prof=False):

        temp = np.atleast_1d(temp)
        rho_susp, rho_in = rhos

        input_flow = u_inputs['Inlet']['vol_flow']
        input_conc = u_inputs['Inlet']['mole_conc']
        input_temp = u_inputs['Inlet']['temp']
        # Thermodynamic properties (basis: slurry volume)
        phi_liq = 1 - self.Solid_1.kv * mu_n[3]

        phis = [phi_liq, 1 - phi_liq]
        h_sp = self.Slurry.getEnthalpy(temp, phis, rho_susp)
        capacitance = self.Slurry.getCp(temp, phis, rho_susp)  # J/m**3/K
        # Heat of rxn
        deltah_ref = self.RxnKinetics.delta_hrxn
        tref_dh = self.RxnKinetics.tref_hrxn

        deltah_rxn = self.Liquid_1.getHeatOfRxn(
            self.RxnKinetics.stoich_matrix, temp, self.mask_species,
            deltah_ref, tref_dh) 
        rates = self.RxnKinetics.get_rxn_rates(mole_conc.T[self.mask_species].T,
                                            temp, overall_rates=False,
                                            delta_hrxn=deltah_rxn)
        # Renaming
        dh_cryst = -1.46e4  # J/kg  # TODO: read this from json file
        # dh_cryst = -self.Liquid_1.delta_fus[self.target_ind] / \
        #     self.Liquid_1.mw[self.target_ind] * 1000  # J/kg

        height_liq = vol / (np.pi/4 * self.diam_tank**2)
        area_ht = np.pi * self.diam_tank * height_liq + self.area_base  # m**2

        # Energy terms (W)
        flow_term = input_flow * (h_in - h_sp)
        cryst_term = dh_cryst*cryst_rate * vol
        rxn_term = -(deltah_rxn * rates).sum(axis=1) * vol * 1000  # mol/Ls * J/mol * vol(m3) * 1000L/m3 -> J/s

        if 'temp' in self.controls.keys():
            ht_term = capacitance * vol  # return capacitance TODO check if works with RC
        elif 'temp' in self.solver_states:
            ht_term = self.u_ht*area_ht*(temp - temp_ht) if not self.isothermal else (flow_term + rxn_term - cryst_term)
        if heat_prof:
            heat_components = np.hstack([cryst_term, ht_term, flow_term, rxn_term])
            return heat_components
        else:
            # Balance inside the tank
            dtemp_dt = (flow_term + rxn_term - cryst_term - ht_term) / vol / capacitance

            # Balance in the jacket
            ht_media = self.Utility.get_inputs(time)
            flow_ht = ht_media['vol_flow']
            tht_in = ht_media['temp_in']

            cp_ht = self.Utility.cp
            rho_ht = self.Utility.rho

            vol_ht = self.vol_tank*0.14  # m**3

            dtht_dt = flow_ht / vol_ht * (tht_in - temp_ht) - \
                self.u_ht*area_ht*(temp_ht - temp) / rho_ht/vol_ht/cp_ht

            return dtemp_dt, dtht_dt
        
    def retrieve_results(self, time, states):
        time = np.array(time)
        # self.checker = progress_checker(max_count=30000, flag='retrieve')
        # ---------- Create result object
        inputs = self.get_inputs(time)
        volflow = inputs['Inlet']['vol_flow']

        dp = unpack_states(states, self.dim_states, self.name_states)
        if self.basis == 'mass_j':
            dp['mole_j'] = dp['mass_j']/self.Liquid_1.mw*1000
            dp['liq_mass_tot'] = np.sum(dp['mass_j'],axis=1)
            dp['liq_moles_tot'] = np.sum(dp['mole_j'],axis=1)
            dp['mass_frac'] = dp['mass_j']/dp['liq_mass_tot'][:,None]
            dp['mole_frac'] = dp['mole_j']/dp['liq_moles_tot'][:,None]
            dp['vol'] = np.sum(dp['mass_j'] / self.Liquid_1.getDensityPure()[0],axis=1)
            dp['rho'] = dp['liq_mass_tot']/dp['vol']
            dp['mass_conc'] = dp['mass_j']/dp['vol'][:,None]

        dp['time'] = time
        dp['vol_flow'] = volflow
        dp['x_cryst'] = self.x_grid
        dp['mole_conc']= dp['mass_conc']/self.Liquid_1.mw
        
        if 'temp' in self.controls:
            control = self.controls['temp']
            dp['temp'] = control['fun'](time, *control['args'], **control['kwargs'])
        mass_conc_sat = dp['mass_conc'] if self.basis!='mass_j' else dp['mass_j']/dp['mass_j'][:,self.Liquid_1.ind_solv][:, np.newaxis]*self.Liquid_1.getDensityPure()[0][self.Liquid_1.ind_solv]
        sat_conc = self.CrystKinetics.get_solubility(dp['temp'], mass_conc_sat)

        supersat = mass_conc_sat[:, self.target_ind] - sat_conc 
        if self.CrystKinetics.sup_sat_type == 'relative':
            supersat = supersat / sat_conc

        if self.CrystKinetics.sup_sat_type == 'ratio':
            supersat = supersat / sat_conc + 1

        dp['solubility'] = sat_conc
        dp['supersat'] = supersat

        if self.population_balance_method == '1D-FVM':
            dp['distrib'] *= 1 / self.scale
            moms = self.Solid_1.getMoments(distrib=dp['distrib'])
            dp['mu_n'] = moms

            dp['vol_distrib'] = self.Solid_1.convert_distribution(
                num_distr=dp['distrib'])

            if type(self) == ReactiveMSMPR:
                vol_slurry = self.Slurry.vol
                self.Solid_1.updatePhase(distrib=dp['distrib'][-1] * vol_slurry)

        if self.population_balance_method == 'moments':
            dp['mu_n'] = dp['mu_n'] * (1e-6)**np.arange(self.num_distr)

        if 'semi' in self.__class__.__name__.lower() :
            dp['total_distrib'] = dp['distrib']

        self.profiles_runs.append(dp)
        dp = self.flatten_states()
        self.get_heat_duty(time, states,4)  # TODO: allow for semi-batch
        dp['q_rxn'] = self.rxn_heat_prof[:,0]
        dp['q_ht'] = self.rxn_heat_prof[:,1]

        self.outputs = dp
        


        # ---------- Update phases

        self.Solid_1.temp = dp['temp'][-1]
        self.Liquid_1.temp = dp['temp'][-1]

        if type(self) == ReactiveMSMPR:
            vol_slurry = self.Slurry.vol
            vol_liq = (1 - self.Solid_1.kv * dp['mu_n'][-1, 3]) * vol_slurry

            self.Liquid_1.updatePhase(vol=vol_liq,
                                      mole_conc=dp['mole_conc'][-1])
            if self.population_balance_method == '1D-FVM':
                distrib_tilde = dp['distrib'][-1] * vol_slurry
                self.Solid_1.updatePhase(distrib=distrib_tilde)

                self.Slurry = Slurry()

            elif self.population_balance_method == 'moments':
                self.Slurry = Slurry(moments=dp['mu_n'][-1], vol=vol_slurry)

        else:
            vol_liq = dp['vol'][-1]
            if self.basis!='mass_j':
                self.Liquid_1.updatePhase(mass_conc=dp['mass_conc'][-1],
                                    vol=dp['vol'][-1], solvent_pass=True)
            else:
                self.Liquid_1.updatePhase(mass=dp['liq_mass_tot'][-1],
                                          mass_frac=dp['mass_frac'][-1],
                                          solvent_pass=True)
            
            rho_solid = self.Solid_1.getDensity()
            vol_solid = dp['mu_n'][-1, 3] * self.Solid_1.kv
            mass_solid = rho_solid*vol_solid


            vol_slurry = vol_solid + vol_liq
            
            if self.population_balance_method == '1D-FVM':
                distrib_tilde = dp['total_distrib'][-1]
                self.Solid_1.updatePhase(distrib=distrib_tilde,
                                         mass= mass_solid)

                self.Slurry = Slurry()

            elif self.population_balance_method == 'moments':
                pass  # TODO

        self.Slurry.Phases = (self.Solid_1, self.Liquid_1)
        self.elapsed_time = time[-1]

        # ---------- Create output stream
        path = self.Liquid_1.path_data

        solid_comp = np.zeros(self.num_species)
        solid_comp[self.target_ind] = 1

        if type(self) == ReactiveMSMPR:
            liquid_out = LiquidStream(path,
                                      mass_conc=dp['mass_conc'][-1],
                                      temp=dp['temp'][-1], check_input=False)

            solid_out = SolidStream(path, mass_frac=solid_comp)

            if isinstance(inputs['Inlet']['vol_flow'], float):
                vol_flow = inputs['Inlet']['vol_flow']
            else:
                vol_flow = inputs['Inlet']['vol_flow'][-1]

            if self.population_balance_method == '1D-FVM':
                # check cstr for semibatch here
                self.Outlet = SlurryStream(
                    vol_flow=vol_flow,
                    x_distrib=self.x_grid,
                    distrib=dp['distrib'][-1])

            elif self.population_balance_method == 'moments':

                self.Outlet = SlurryStream(
                    vol_flow=vol_flow,
                    moments=dp['mu_n'][-1])


        else:
            liquid_out = copy.deepcopy(self.Liquid_1)
            solid_out = copy.deepcopy(self.Solid_1)

            self.Outlet = Slurry(vol=vol_slurry)

        # self.outputs = y_outputs
        self.Outlet.Phases = (liquid_out, solid_out)

        # check that balances close
        if False:
            self.closure = np.zeros((len(dp['mass_conc']),3)) #mass vol
            for i in range(len(dp['mass_conc'])):
                liquid_out_close = LiquidStream(path,
                                        mass_conc=dp['mass_conc'][i],
                                        temp=dp['temp'][i], check_input=False)

                solid_out_close = SolidStream(path, mass_frac=solid_comp)

                if isinstance(inputs['Inlet']['vol_flow'], float):
                    vol_flow_close = inputs['Inlet']['vol_flow']
                else:
                    vol_flow_close = inputs['Inlet']['vol_flow'][i]
                
                if isinstance(self.Inlet.mass_flow,float):
                    mass_in_close = self.Inlet.mass_flow
                else:
                    mass_in_close = self.Inlet.mass_flow[i]
                


                
                if self.population_balance_method == '1D-FVM':
                    # check cstr for semibatch here
                    Outlet_close = SlurryStream(
                        vol_flow=vol_flow_close,
                        x_distrib=self.x_grid,
                        distrib=dp['distrib'][i])

                elif self.population_balance_method == 'moments':

                    Outlet_close = SlurryStream(
                        vol_flow=vol_flow,
                        moments=dp['mu_n'][i])
                Outlet_close.Phases = (liquid_out_close, solid_out_close)

                # dmapi_liq__dt = self.Inlet['mass_conc'][self.target_ind]*vol_flow_close - dp['mass_conc'][self.target_ind]*vol_flow_close - 
                mass_out = Outlet_close.Solid_1.mass_flow#self.kin_array[time[i]][-1]#Outlet_close.Solid_1.mass*Outlet_close.Solid_1.mass_frac[self.target_ind]
                self.closure[i] = np.array((Outlet_close.mass_flow/mass_in_close, Outlet_close.vol_flow/vol_flow_close, mass_out))

            tot_mass = [sum(self.closure[:i,-1]) for i in range(len(self.closure))]
            dp['m_flow'] = self.closure[:,-1]
            dp['tot_mass_cryst'] = np.array(tot_mass)
            dp['tot_mass_2'] = dp['mu_n'][:,3]*self.Solid_1.kv*vol_flow_close*self.Solid_1.getDensity() # add to fstates

        self.result = DynamicResult(self.states_di, self.fstates_di, **dp) 



    def get_heat_duty(self, time, states,n_components=3):
        q_heat = np.zeros((len(time), n_components))

        if self.params_iter is None:
            merged_params = self.CrystKinetics.concat_params()[self.mask_params_cryst]
        else:
            merged_params = self.params_iter
        if not self.param_estimation_run:
            for ind, row in enumerate(states):
                row = row.copy()
                row[:self.num_distr] *= self.scale  # scale distribution
                q_heat[ind] = self.unit_model(time[ind], row, merged_params,
                                            enrgy_bce=True)

        # q_heat[:, 0] *= -1
        q_cryst, q_ht, flow_term, q_rxn = q_heat.T  # TODO: controlled temperature

        self.heat_prof = q_heat
        self.rxn_heat_prof = np.column_stack((q_rxn,-q_ht,flow_term))
        self.heat_duty = np.array([0, trapezoidal_rule(time, q_ht)])
        self.duty_type = [0, -2]

    

class ReactiveSemibatchCrystallizer(ReactiveMSMPR):
    """
        Assumes:
            constant solid density
            metric units
        
    """
    def __init__(self, target_comp, mask_params_rxn=None,mask_params_cryst=None, temp_ref=298.15, isothermal=True,
                  reset_states=False, controls=None, h_conv=1000, ht_mode='jacket',
                  return_sens=True, state_events=None, method='1D-FVM',
                  scale=1, vol_tank=None, adiabatic=False, rad_zero=0, vol_ht=None,
                  basis='mass_j', jac_type=None, param_wrapper=None, num_interp_points=3, grid_size=500,param_estimation_run=False):
        super().__init__(target_comp,mask_params_rxn,mask_params_cryst, temp_ref, isothermal, reset_states, controls, h_conv, ht_mode, return_sens, state_events, method, scale, vol_tank, adiabatic, rad_zero, vol_ht, basis, jac_type, param_wrapper,param_estimation_run=param_estimation_run)
        self.is_continuous = False
        self.oper_mode = 'Semibatch'
        self.vol_ht = vol_tank*0.15 # this is hard-coded geometry
        self._Inlet = None
        # self.vol_offset = 0.75
        # self.num_interp_points = num_interp_points
        self.mydistrib = np.zeros(grid_size)
        # self.checker = progress_checker(flag='Solver', minh=0.01)
        self.kin_array = {}
 
   
    def quickflow(self,val,alt=0):
        if not self._zero_flow:
            return val
        else:
            return alt
    def material_balances(self, time, params, u_inputs, rhos, mu_n,
                          distrib, mass_j, temp, temp_ht, phi_in):

        rho_susp, rho_in = rhos 
        # self.checker.check(time)
        rho_liq, rho_sol = rho_susp
        rho_in_liq, _ = rho_in
        
        
        input_flow = u_inputs['Inlet']['vol_flow']
        input_flow = self.quickflow(np.max([eps, input_flow]))

        # TODO: generalize dictionary iteration ('Inlet', 'Liquid_1', ...)?
        input_distrib = self.quickflow(u_inputs['Inlet']['distrib'] * self.scale, np.zeros_like(distrib))
        input_conc = self.quickflow(u_inputs['Liquid_1']['mass_conc'], np.zeros_like(u_inputs['Liquid_1']['mass_conc']))

        # print('time = %.2f, vol = %.2e, flowrate = %.2e' % (time, vol, input_flow))

        vol_solid = mu_n[3] * self.Solid_1.kv  # mu_3 is total, not by volume
        mole_j = mass_j/self.Liquid_1.mw*1000 # moles
        tot_mass = sum(mass_j)
        tot_moles = sum(mole_j)
        mass_frac = mass_j/tot_mass
        mole_frac = mole_j/tot_moles
        ## for debugging
        # raw_vol = np.sum(mass_j / self.Liquid_1.getDensityPure()[0])
        # raw_rho =mass_j.sum() / np.sum(mass_j / self.Liquid_1.getDensityPure()[0])

        # print(time)
        # print(vol*rho_liq, vol, rho_liq, sum(mass_conc)*vol)
        # print(self.Liquid_1.mass,self.Liquid_1.vol, self.Liquid_1.getDensity(),sum(self.Liquid_1.mass_frac)*self.Liquid_1.vol)
        # print()
        self.Liquid_1.updatePhase(mass_frac=mass_frac,mass=tot_mass,solvent_pass=True)
        mole_conc = mole_j/self.Liquid_1.vol/1000 #kgmol/m3  mole/L
        vol_slurry = self.Liquid_1.vol + vol_solid
        mass_conc_per_solvent_vol = mass_j/mass_j[self.Liquid_1.ind_solv]*self.Liquid_1.getDensityPure()[0][self.Liquid_1.ind_solv]
        # note that transf is mass of liquid transferred in kg
        if self.population_balance_method == 'moments':
            ddistr_dt, transf = self.method_of_moments(distrib, mass_conc_per_solvent_vol, temp,
                                                       params, rho_sol,
                                                       vol=vol_slurry)

        elif self.population_balance_method == '1D-FVM':
            ddistr_dt, transf = self.fvm_method(distrib, mu_n, mass_conc_per_solvent_vol, temp,
                                                params, rho_sol,
                                                vol=vol_slurry)

        # ---------- Add flow terms
        # Distribution
        flow_distrib = input_flow * input_distrib

        ddistr_dt = ddistr_dt + flow_distrib

        # Liquid phase
        # c_tank = mass_conc
        # mole_conc = mass_conc /self.Liquid_1.mw
        #### From Reactor ####
        if self.RxnKinetics.keq_params is None:
            rate = self.RxnKinetics.get_rxn_rates(mole_conc[self.mask_species], temp) #rates in kmol/m^3_reactor_liquid negative=consumed 
        else:
            concentr = np.zeros(len(self.name_species))
            concentr[self.mask_species] = mole_conc[self.mask_species]
            concentr[~self.mask_species] = self.conc_inert
            deltah_rxn = self.Liquid_1.getHeatOfRxn(
                self.RxnKinetics.stoich_matrix, temp, self.mask_species,
                self.RxnKinetics.delta_hrxn, self.RxnKinetics.tref_hrxn)

            rate = self.RxnKinetics.get_rxn_rates(mole_conc[self.mask_species], temp,
                                               delta_hrxn=deltah_rxn)
        rates = np.zeros_like(mole_conc)
        rates[self.mask_species] = rate
        rates*=self.Liquid_1.mw #from kmol_j/m^3_reactor_liquid to kg_j/m^3_reactor_liquid
        #### End from Reactor ####
        ####
        rxnterm= rates*self.Liquid_1.vol
        dmj_dt = phi_in[0]*input_flow*input_conc + rxnterm - transf*self.kron_jtg
        

        dmaterial_dt = np.concatenate((ddistr_dt, dmj_dt))

        return dmaterial_dt, transf

    def energy_balances(self, time, params, cryst_rate, u_inputs, rhos,
                        distrib, mass_j, temp, temp_ht, mu_n, h_in,heat_prof=False):

        rho_susp, rho_in = rhos

        # Input properties
        vol_solid = mu_n[3] * self.Solid_1.kv  # mu_3 is total, not by volume
        mole_j = mass_j/self.Liquid_1.mw*1000 # moles
        tot_mass = sum(mass_j)
        tot_moles = sum(mole_j)
        mass_frac = mass_j/tot_mass
        mole_frac = mole_j/tot_moles
        self.Liquid_1.updatePhase(mass_frac=mass_frac,mass=tot_mass,solvent_pass=True)
        
        input_flow = u_inputs['Inlet']['vol_flow']
        input_flow = self.quickflow(np.max([eps, input_flow]))
        vol = self.Liquid_1.vol
        vol_solid = mu_n[3] * self.Solid_1.kv  # mu_3 is total, not by volume
        vol_total = vol + vol_solid

        phi = vol / vol_total
        phis = [phi, 1 - phi]
        dens_slurry = np.dot(rho_susp, phis)

        # Suspension properties
        capacitance = self.Slurry.getCp(temp, phis, rho_susp,
                                        times_vliq=True)
        h_sp = self.Slurry.getEnthalpy(temp, phis, rho_susp)

        # Renaming
        dh_cryst = -1.46e4  # J/kg
        ##### From Reaction#####
        delta_href = self.RxnKinetics.delta_hrxn
        stoich = self.RxnKinetics.stoich_matrix
        tref_hrxn = self.RxnKinetics.tref_hrxn
        deltah_rxn = self.Liquid_1.getHeatOfRxn(
            stoich, temp, self.mask_species, delta_href, tref_hrxn)  # J/mol

        rxn_rates = self.RxnKinetics.get_rxn_rates(self.Liquid_1.mole_conc[self.mask_species], temp,
                                            overall_rates=False,
                                            delta_hrxn=deltah_rxn)
        rxn_term = -(deltah_rxn * rxn_rates).sum() *vol*1000 # TODO check if need mult by 1000
        ###### End Reaction terms #####
        # dh_cryst = -self.Liquid_1.delta_fus[self.target_ind] / \
        #     self.Liquid_1.mw[self.target_ind] * 1000  # J/kg

        # Terms
        dens_in_liq = rho_in[0]
        dmass_dt = input_flow * dens_in_liq

        accum_term = dmass_dt * h_sp/dens_slurry
        flow_term = input_flow * h_in

        cryst_term = dh_cryst * cryst_rate*vol

        height_liq = vol / (np.pi/4 * self.diam_tank**2)
        area_ht = np.pi * self.diam_tank * height_liq + self.area_base  # m**2

        if 'temp' in self.controls.keys():
            ht_term = capacitance * vol  # return capacitance TODO check if works with RC
        elif 'temp' in self.solver_states:
            ht_term = self.u_ht*area_ht*(temp - temp_ht) if not self.isothermal else (flow_term + rxn_term - cryst_term)
        if heat_prof:
            heat_components = np.hstack([cryst_term, ht_term, flow_term, rxn_term])
            return heat_components
        else:
            # Balance inside the tank
            dtemp_dt = (flow_term + rxn_term - cryst_term - ht_term) / vol / capacitance

            # Balance in the jacket
            ht_media = self.Utility.get_inputs(time)
            flow_ht = ht_media['vol_flow']
            tht_in = ht_media['temp_in']

            cp_ht = self.Utility.cp
            rho_ht = self.Utility.rho

            vol_ht = self.vol_tank*0.14  # m**3

            dtht_dt = flow_ht / vol_ht * (tht_in - temp_ht) - \
                self.u_ht*area_ht*(temp_ht - temp) / rho_ht/vol_ht/cp_ht

            return dtemp_dt, dtht_dt
        

if __name__=='__main__':
    a = MultiPhaseVessel('C',300,False,True,{},0,True,{},TransferMechanism,1,False,'AD')
    dpath = r"C:\Users\zhillma\OneDrivePZH\Documents\Documents\_Grad_School\mypharmadev\PharmaPy\tests\Flowsheet\data\compound_database.json"
    liquid1 = LiquidPhase(dpath,mass=1,mass_frac=[.3,.7,0,0,0])
    solid1 = SolidPhase(dpath,mass=0,mass_frac=[0,0,1,0,0])
    a.Phases = [liquid1,solid1]
    a.phase_connections
