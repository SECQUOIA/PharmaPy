
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
from typing import Optional, Sequence
from collections import OrderedDict


eps = np.finfo(float).eps
# gas_ct = 8.314  # J/mol/K
class TransferMechanism:
    def add_state_variables(self,collection):
        pass
    def add_output_state_variables(self, outputs):
        pass
    def transfer(self,mass_j):
        return mass_j
class DirectTransfer(TransferMechanism):
    pass
def build_transfer_vectors(paths):

    vectors = []

    for path in paths:
        vec = np.asarray(path.species_weights)

        vec = vec / vec.sum()

        vectors.append(vec)

    return vectors
@dataclass(frozen=True)
class PhaseRef:
    phase_type: str
    index: int

@dataclass
class PhaseConnection:
    #TODO move this to connections when done
    #active_condition checks the mass_j and temp of the source_phase and must return a boolean
    source_phase: PhaseRef
    sink_phase: PhaseRef
    kinetics:pk.CrystKinetics|pk.RxnKinetics
    species_weights: np.ndarray | None = None
    active_condition: callable=lambda mass_j,temp:True
    mechanism:TransferMechanism|callable=DirectTransfer()

@dataclass
class ReactionRegion:
    phase:PhaseRef
    kinetics:pk.RxnKinetics

@dataclass
class StateVariable:
    name: str
    dim: int
    units: str
    state_type: str = "diff"
    index: Optional[Sequence] = None
    depends_on: tuple = ("time",)
    stream:str|None=None
    phase: PhaseRef | None = None

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

@dataclass
class StateCollection:
    states: dict[str, StateVariable] = field(default_factory=dict)

    def add(self, state: StateVariable):
        if state.name in self.states:
            raise ValueError(f"State {state.name} already exists")
        self.states[state.name] = state

    def names(self):
        return list(self.states.keys())

    def dims(self):
        return [state.dim for state in self.states.values()]

    def __contains__(self, name):
        return name in self.states

@dataclass
class PhaseStateVariable:
    phase: PhaseRef
    state: StateVariable

@dataclass
class PhaseStateCollection:
    states: dict[PhaseRef, StateCollection] = field(default_factory=dict)

    def add(self, phase: PhaseRef, state: StateVariable):
        if phase not in self.states:
            self.states[phase] = StateCollection()

        self.states[phase].add(state)

    def __getitem__(self, phase):
        return self.states[phase]
    def __iter__(self):
        for phase, collection in self.states.items():
            for state in collection.states.values():
                yield PhaseStateVariable(phase, state)

class MultiPhaseVessel():
    def __init__(self,target_comp, temp_ref,
     isothermal, reset_states, controls, h_conv, ht_mode,
      state_events,population_balance_method,scale,
      adiabatic,jac_type,
      basis='mass_j'):
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
        self._initialize_states(reset_states)
        
        # Phase init
        self._Phases = None
        self.Slurry = None
        self.phase_connections = []

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

        #Reaction
        self.temp_ref = temp_ref
        self._RxnKinetics = None
        self._reaction_regions = []
        
        #State events
        if state_events is None:
            state_events = []
        self.state_event_list = state_events


    @property
    def Phases(self):
        return self._Phases
    
    @Phases.setter
    def Phases(self, phases):
        self._normalize_phases(phases)
        self.__original_phase_dict__ = [copy.deepcopy(phase.__dict__) for phase in self._Phases]
        self._post_set_phases()
        
    def _normalize_phases(self,phases):
        if isinstance(phases, (list, tuple)): #TODO ZZ _phases should now always be
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
    # def _post_set_phases(self):
        # pass
    def _basis_units(self):

        units = {
            "mass_j": "kg",
            "mass_conc": "kg/m3",
            "mole_j": "kmol",
            "mole_conc": "kmol/m3"
        }

        return units[self.basis]
    def _material_state_definition(self):

        return StateVariable(
            name=self.basis,
            dim=self.num_species,
            units=self._basis_units(),
            index=self.name_species
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
        self._initialize_phase_states()
        # # Input defaults TODO delete this if not necessary
        # self.input_defaults = {
        #     'distrib': np.zeros_like(self.Solid_1.distrib)}
        # Species
        self.define_material_states()
        self.nomenclature() 

    def _initialize_phase_states(self):

        self.phase_states = PhaseStateCollection()

        # Stream variables are not phases. They describe inlet/outlet connections.
        self.stream_states = StateCollection()

        # States exposed to solver
        self.solver_state_collection = StateCollection()

        # States exposed as outputs
        self.output_state_collection = StateCollection()

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
                                             active_condition=lambda mass_j,temp:True,
                                             mechanism=self.Solid_1.transfer_mechanism
                                             )
                connections.append(connection)
            if ck.supports('dissolution'):
                #solid to liquid because dissolution
                connection = PhaseConnection(source_phase=PhaseRef("solid",0),
                                             sink_phase=PhaseRef('liquid',0),
                                             kinetics=ck,
                                             species_weights=weights,
                                             active_condition=lambda mass_j,temp:True,
                                             mechanism=DirectTransfer()

                                             )
                connections.append(connection)
        

    @property
    def reaction_regions(self):
        return self._reaction_regions
    @reaction_regions.setter
    def reaction_regions(self,regions):
        if not all([isinstance(r,ReactionRegion) for r in regions]):
            raise TypeError("reaction_regions should all be ReactionRegion objects")
        if not isinstance(regions,list):
            raise TypeError("reaction_regions is expected to be a list")
        self._reaction_regions = regions
        if len(self.reaction_regions) == 0:
            self.mask_species = np.ones(
                self.num_species,
                dtype=bool
            )
            return

        participating_species = set()

        for region in self.reaction_regions:

            rk = region.kinetics

            if hasattr(rk, "partic_species"):
                participating_species.update(
                    rk.partic_species
                )

        self.mask_species = np.array(
            [
                species in participating_species
                for species in self.name_species
            ]
        )

        self.conc_inert = np.zeros(self.num_species)

        for i, active in enumerate(self.mask_species):
            if not active:
                self.conc_inert[i] = getattr(
                    self.Phases[0],
                    self.basis
                )[i]
        
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

        self.solver_states = []

        self.elapsed_time = 0

        self.outputs = None
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
    def define_solver_states(self):

        
        for phase_state in self.phase_states:

            self.solver_state_collection.add(
                copy.deepcopy(phase_state.state)
            )
        
        if self.adiabatic:

            self.solver_state_collection.add(
                StateVariable(
                    name="temp",
                    dim=1,
                    units="K"
                )
            )
            
        elif "temp" not in self.controls:

            self.solver_state_collection.add(
                StateVariable(
                    name="temp",
                    dim=1,
                    units="K"
                )
            )
            self.solver_state_collection.add(
                StateVariable(
                    name="temp_ht",
                    dim=1,
                    units="K"
                )
            )
    def define_output_states(self):

        self.output_state_collection.add(
            StateVariable(
                name="q_rxn",
                dim=1,
                units="W",
                state_type="alg"
            )
        )

        self.output_state_collection.add(
            StateVariable(
                name="q_ht",
                dim=1,
                units="W",
                state_type="alg"
            )
        )

        self.output_state_collection.add(
            StateVariable(
                name="m_flow",
                dim=1,
                units="kg/s",
                state_type="alg"
            )
        )

        self.output_state_collection.add(
            StateVariable(
                name="tot_mass_cryst",
                dim=1,
                units="kg",
                state_type="alg"
            )
        )

        for conn in self.phase_connections:

            if conn.mechanism is not None:
                conn.mechanism.add_output_state_variables(
                    self.output_state_collection
                )
    def nomenclature(self):

        self.solver_state_collection = StateCollection()
        self.output_state_collection = StateCollection()
        self.define_solver_states()
        self.define_output_states()
        self.define_stream_states()
        # ------------------------------
        # mechanism states
        # ------------------------------

        for connection in self.phase_connections:

            if connection.mechanism is not None:
                connection.mechanism.add_solver_state_variables(
                    self.solver_state_collection
                )

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
    
    def get_inputs(self, time,solvent_pass=False):

        inlet = getattr(self, "Inlet", None)

        if inlet is None:
            return {}

        inputs = {}

        if hasattr(inlet, "add_stream_state_variables"):

            temp_collection = StateCollection()

            inlet.add_stream_state_variables(
                temp_collection
            )

            for state in temp_collection.states.values():

                value = getattr(inlet, state.name)

                if callable(value):
                    value = value(time)

                inputs[state.name] = np.asarray(value)

        inputs.update(
            get_inputs_new(
                time,
                inlet,
                self.states_in_dict,
                solvent_pass=solvent_pass
                )
            )
        return inputs
        
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
    def method_of_moments(self, mu, conc, temp, params, rho_cry, vol=1):
        kv = self.Solid_1.kv # shape factor

        # Kinetics
        if self.basis == 'mass_frac':
            rho_liq = self.Liquid_1.getDensity()
            comp_kin = conc / rho_liq
        else:
            comp_kin = conc

        # Kinetic terms
        mu_susp = mu*(1e-6)**np.arange(self.num_distr) / vol  # m**n/m**3_susp
        nucl, growth, dissol = self.CrystKinetics.get_kinetics(comp_kin, temp, kv,
                                                          mu_susp)

        growth = growth * self.CrystKinetics.alpha_fn(conc)

        ind_mom = np.arange(1, len(mu))

        # Model
        dmu_zero_dt = np.atleast_1d(nucl * vol)
        dmu_1on_dt = ind_mom * (growth + dissol) * mu[:-1] + \
            nucl * self.rad**ind_mom
        dmu_dt = np.concatenate((dmu_zero_dt, dmu_1on_dt))

        # Material balance in kg_API/s --> G in um, u_2 in um**2 (or m**2/m**3)
        mass_transf = np.atleast_1d(rho_cry * kv * (
            3*(growth + dissol)*mu[2] + nucl*self.rad**3)) * (1e-6)**3

        return dmu_dt, mass_transf

    def fvm_method(self, csd, moms, conc, temp, params, rho_cry,
                   output='dstates', vol=1):

        mu_2 = moms[2]
        #assumes solid1 is target
        kv_cry = self.Solid_1.kv # volumetric shape factor

        # Kinetic terms
        if self.basis == 'mass_frac':
            rho_liq = self.Liquid_1.getDensity()
            comp_kin = conc / rho_liq
        else:
            comp_kin = conc

        nucl, growth, dissol = self.CrystKinetics.get_kinetics(comp_kin, temp,
                                                          kv_cry, moms)

        nucl = nucl * self.scale * vol 

        impurity_factor = self.CrystKinetics.alpha_fn(conc)
        growth = growth * impurity_factor  # um/s 
        gparams = self.CrystKinetics.params['growth']
        

        # dissol = dissol  # um/s
        boundary_cond = nucl / np.maximum(growth, eps) # num/um or num/um/m**3 initial
        f_aug = np.concatenate(([boundary_cond]*2, csd, [csd[-1]])) # TODO adjust for reaction or handled by concentration? 

        # Flux source terms
        f_diff = np.diff(f_aug)
        
        # f_diff[f_diff == 0] = eps  # avoid division by zero for theta

        if growth > 0:
            theta = f_diff[:-1] / (f_diff[1:] + eps*10)
            # theta = f_diff[:-1] / (f_diff[1:] + eps)
            # theta = f_diff[:-1] / f_diff[1:]
        else:
            theta = f_diff[1:] / (f_diff[:-1] + eps*10)
            # theta = f_diff[:-1] / (f_diff[1:] + eps)
            # theta = f_diff[:-1] / f_diff[1:]
        # Van-Leer limiter
        limiter = np.zeros_like(f_diff)
        limiter[:-1] = (np.abs(theta) + theta) / (1 + np.abs(theta))
        if len(gparams)==3:
        
            growth_term = growth * (f_aug[1:-1] + 0.5 * f_diff[1:] * limiter[:-1])
            dissol_term = dissol * (f_aug[2:] - 0.5 * f_diff[1:] * limiter[1:])
        else:
            growth_dependent = growth * (1 + self.x_grid * gparams[4])**gparams[3]
            dissol_dependent = dissol * (1 + self.x_grid * 0) # TODO add size-dependent dissol params
            growth_pad = np.append(growth_dependent,growth_dependent[-1])
            dissol_pad = np.append(dissol_dependent, dissol_dependent[-1])
            growth_term = growth_pad * (f_aug[1:-1] + 0.5 * f_diff[1:] * limiter[:-1])
            dissol_term = dissol_pad * (f_aug[2:] - 0.5 * f_diff[1:] * limiter[1:])
        flux = growth_term + dissol_term

         
        if output == 'flux':
            return flux  # TODO: isn't it necessary to divide by dx?
        elif output=='dstates':
            dcsd_dt = -np.diff(flux) / self.dx

            # Material bce in kg_API/s --> G in um, mu_2 in m**2 (or m**2/m**3)
            # AKA R_v (rho_c*kv*d_mu3_d_t)
            # Handle stoich in material balance
            if len(gparams)==3:
                mass_transfer = rho_cry * kv_cry * (
                    3*(growth + dissol)*mu_2 + nucl*self.rad**3) * (1e-6)
            else:
                r_m = self.x_grid
                mass_transfer_growth = np.trapezoid(growth_dependent*csd*r_m**2,r_m)
                mass_transfer_dissol = np.trapezoid(dissol_dependent*csd*r_m**2,r_m)
                mass_transfer_nucl = nucl*self.rad**3
                mass_transfer = rho_cry*kv_cry*3*(mass_transfer_dissol+mass_transfer_growth+mass_transfer_nucl)*1e-18
            return dcsd_dt, np.array(mass_transfer)
        
    def unit_model(self, time, states, params=None, sw=None,
                    mat_bce=False, enrgy_bce=False):
        # TODO reconcile with RC
        di_states = unpack_states(states, self.dim_states, self.name_states)
        if self.basis!='mass_j':
            di_states['mole_conc'] = di_states['mass_conc']/self.Liquid_1.mw 
            di_states = complete_dict_states(time, di_states,
                                            ('temp', 'temp_ht', 'vol'),
                                            self.Slurry, self.controls) # this is used to update when there are controls
                # ---------- Physical properties
            # self.Liquid_1.updatePhase(vol=di_states['vol'])
            self.Liquid_1.updatePhase(mass_conc=di_states['mass_conc'],vol=di_states['vol'], solvent_pass=True)
            self.Liquid_1.temp = di_states['temp']
            self.Solid_1.temp = di_states['temp']

            rhos_susp = self.Slurry.getDensity(temp=di_states['temp'])
            rhos_susp[0] = sum(di_states['mass_conc'])
        else:
            di_states = complete_dict_states(time, di_states,
                                            ('temp', 'temp_ht'),
                                            self.Slurry, self.controls) # this is used to update when there are controls
            # ---------- Physical properties
            # self.Liquid_1.updatePhase(vol=di_states['vol'])
            tot_mass = sum(di_states['mass_j'])
            mass_frac = di_states['mass_j']/tot_mass
            self.Liquid_1.updatePhase(mass_frac=mass_frac,mass=tot_mass, solvent_pass=True)
            self.Liquid_1.temp = di_states['temp']
            self.Solid_1.temp = di_states['temp']
            # di_states['mole_conc']

            rhos_susp = self.Slurry.getDensity(temp=di_states['temp'])
            # rhos_susp[0] = sum(di_states['mass_conc'])
        # Inputs
        u_input = self.get_inputs(time,solvent_pass=True)


        
        name_unit = self.__class__.__name__

        if self.population_balance_method == 'moments':
            di_states['distrib'] = di_states['mu_n']
            moms = di_states['mu_n'] * \
                (1e-6)**np.arange(self.states_di['mu_n']['dim']) ###units

        else:
            moms = self.Solid_1.getMoments(
                distrib=di_states['distrib']/self.scale)  # m**n

        di_states['mu_n'] = moms
        batched = ('batch' in name_unit.lower() and 'semi' not in name_unit.lower())
        try:
            self._zero_flow = True if u_input['Inlet']['vol_flow'] == 0 else False
        except KeyError:
            self._zero_flow = True if u_input['Inlet']['mass_flow'] == 0 else False
        if batched: 
            rhos = rhos_susp
            h_in = None
            phis_in = None
        elif 'semi' in name_unit.lower() and self._zero_flow:
            rhos = [rhos_susp, np.zeros_like(self.Inlet.getDensity(temp=di_states['temp']))]
            mom_in = np.zeros(1)
            phi_in = 1- self.Inlet.Solid_1.kv * mom_in
            phis_in = np.concatenate([phi_in, 1 - phi_in]) # TODO assumes only two pahses

            h_in = 0
        elif 'semi' in name_unit.lower() or 'msmpr' in name_unit.lower():
            inlet_temp = u_input['Inlet']['temp']

            if self.Inlet.__module__ == 'PharmaPy.MixedPhases':
                rhos_in = self.Inlet.getDensity(temp=di_states['temp'])

                if 'distrib' in u_input['Inlet']:

                    inlet_distr = u_input['Inlet']['distrib']

                    mom_in = self.Inlet.Solid_1.getMoments(distrib=inlet_distr,
                                                            mom_num=3)
                elif 'mu_n' in u_input['Inlet']:

                    mom_in = np.array([u_input['Inlet']['mu_n'][3]])


                phi_in = 1 - self.Inlet.Solid_1.kv * mom_in
                phis_in = np.concatenate([phi_in, 1 - phi_in]) # TODO assumes only two pahses

                h_in = self.Inlet.getEnthalpy(inlet_temp, phis_in, rhos_in)
            else:
                rho_liq_in = self.Inlet.getDensity(temp=inlet_temp)
                rho_sol_in = None

                rhos_in = np.array([rho_liq_in, rho_sol_in])
                h_in = self.Inlet.getEnthalpy(temp=inlet_temp)

                phis_in = [1, 0]

            rhos = [rhos_susp, rhos_in]

        # Balances
        material_bces, cryst_rate = self.material_balances(
            time, params, u_input, rhos, **di_states, phi_in=phis_in)
        

        if mat_bce:
            return material_bces
        elif enrgy_bce:
            energy_bce = self.energy_balances(
                time, params, cryst_rate, u_input, rhos, **di_states,
                h_in=h_in, heat_prof=True)

            return energy_bce
        #equal here with ~R386
        else:

            if 'temp' in self.name_states:
                energy_bce = self.energy_balances(
                    time, params, cryst_rate, u_input, rhos, **di_states,
                    h_in=h_in)

                balances = np.append(material_bces, energy_bce)
            else:
                balances = material_bces

            self.derivatives = balances

            
            return balances
        # TODO jsut fix multiple states if needed else good
    

    def unit_jacobians(self, time, states, sens, params, fy, v_vector):
        if sens is not None:
            jac_states = self.jac_states_fun(time, states, params)
            jac_params = self.jac_params_fun(time, states, params)

            dsens_dt = np.dot(jac_states, sens) + jac_params

            if not isinstance(dsens_dt, np.ndarray):
                dsens_dt = dsens_dt._value

            return dsens_dt
        elif v_vector is not None:
            _, jac_v = self.jac_states_fun(time, states, params)(v_vector)

            return jac_v
        else:
            jac_states = self.jac_states_fun(time, states, params)

            if not isinstance(jac_states, np.ndarray):
                jac_states = jac_states._value

            return jac_states

    def jac_states_numerical(self, time, states, params, return_only=True):
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
        states_init, merged_params = self._build_initial_state_and_params()

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
        if eval_sens:
            problem = Explicit_Problem(self.unit_model, states_init,
                                       t0=self.elapsed_time,
                                       p0=params_mergd)

            if self.jac_type == 'finite_diff':
                self.jac_states_fn = self.jac_states_numerical
                self.jac_params_fn = self.jac_params_numerical

                problem.jac = self.jac_states_fn
                problem.rhs_sens = self.rhs_sensitivity

            elif self.jac_type == 'AD':
                self.jac_states_fn = self.jac_states_ad
                self.jac_params_fn = self.jac_params_ad

                problem.jac = self.jac_states_fn
                problem.rhs_sens = self.rhs_sensitivity

            elif self.jac_type == 'analytical':
                self.jac_states_fn = self.jac_states
                self.jac_params_fn = self.jac_params

                problem.jac = self.jac_states_fn
                problem.rhs_sens = self.rhs_sensitivity

            elif self.jac_type is None:
                pass
            else:
                raise NameError("Bad string value for the 'jac_type' argument")

        else:
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
        self.set_names()

        if self.population_balance_method == 'moments':
            pass  # TODO: MSMPR MoM should be addressed?
        else:
            x_distr = getattr(self.Solid_1, 'x_distrib', [])
            self.states_in_dict['Inlet']['distrib'] = len(x_distr)

        self.CrystKinetics.target_idx = self.target_ind

        # ---------- Solid phase states
        if 'vol' in self.solver_states or 'mass_j' in self.solver_states:
            if self.population_balance_method == 'moments':
                init_solid = self.Solid_1.moments
                # exp = np.arange(0, self.Solid_1.num_mom) # TODO: problematic line for seeded crystallization.
                # init_solid = init_solid * (1e6)**exp

            elif self.population_balance_method == '1D-FVM':
                x_grid = self.Solid_1.x_distrib
                init_solid = self.Solid_1.distrib * self.scale

        else:
            if self.population_balance_method == 'moments':
                init_solid = self.Slurry.moments
                # exp = np.arange(0, self.Solid_1.num_mom) # TODO
                # init_solid = init_solid * (1e6)**exp

            elif self.population_balance_method == '1D-FVM':
                x_grid = self.Slurry.x_distrib
                init_solid = self.Slurry.distrib * self.scale

        self.dx = self.Slurry.dx
        self.x_grid = self.Slurry.x_distrib

        # ---------- Liquid phase states
        init_liquid = self.Liquid_1.mass_conc.copy()
        if self.basis == 'mass_j':
            init_liquid*= self.Liquid_1.vol

        self.num_species = len(init_liquid)

        self.len_states = [self.num_distr, self.num_species]  # TODO: not neces

        if 'vol' in self.solver_states:  # Batch or semibatch
            vol_init = self.Slurry.getTotalVol()
            init_susp = np.append(init_liquid, vol_init)

            self.len_states.append(1)
        else:
            init_susp = init_liquid

        if self.reset_states:
            self.reset()

        # ---------- Read time
        # if runtime is not None:
        #     final_time = runtime + self.elapsed_time

        # if time_grid is not None:
        #     final_time = time_grid[-1]

        

        states_init = np.append(init_solid, init_susp) # solid=fvm bins/moments susp = liquid concentrations/volume balance

        # if self.vol_tank is None:
        #     if isinstance(self, ReactiveSemibatchCrystallizer):
        #         time_vec = np.linspace(self.elapsed_time, final_time)
        #         vol_flow = self.get_inputs(time_vec)['Inlet']['vol_flow']

        #         self.vol_tank = trapezoidal_rule(time_vec, vol_flow)

        #     else:
        #         self.vol_tank = self.Slurry.vol

        self.diam_tank = (4/np.pi * self.vol_tank)**(1/3) # TODO ensure redefinition is fine Z
        self.area_base = np.pi/4 * self.diam_tank**2
        self.vol_tank *= 1 / self.vol_offset

        if 'temp_ht' in self.solver_states:

            if len(self.profiles_runs) == 0:
                temp_ht = self.Utility.evaluate_inputs(0)['temp_in']
            else:
                temp_ht = self.profiles_runs[-1]['temp_ht'][-1]

            states_init = np.concatenate(
                (states_init, [self.Liquid_1.temp, temp_ht]))

            self.len_states += [1, 1]
        elif 'temp' in self.solver_states:
            states_init = np.append(states_init, self.Liquid_1.temp)
            self.len_states += [1]

        merged_params = self.CrystKinetics.concat_params()[self.mask_params_cryst]
        return states_init, merged_params
    
    def _rebase_initial(self, state_vec):
        """
        Rebase all internal phase objects to match a provided state vector.
        This allows restarting simulations consistently from saved states.
        """

        idx = 0

        # --------------------------------------------------
        # 1) Solid distribution / moments
        # --------------------------------------------------
        nd = self.num_distr
        solid_part = state_vec[idx:idx + nd]
        idx += nd

        if self.population_balance_method == 'moments':
            # Directly update moments
            if 'vol' in self.solver_states or 'mass_j' in self.solver_states:
                self.Solid_1.moments = solid_part.copy()
            else:
                self.Slurry.moments = solid_part.copy()

        elif self.population_balance_method == '1D-FVM':
            # Unscale if needed
            solid_part = solid_part / self.scale

            if 'vol' in self.solver_states or 'mass_j' in self.solver_states:
                self.Solid_1.distrib = solid_part.copy()
            else:
                self.Slurry.distrib = solid_part.copy()

        # --------------------------------------------------
        # 2) Liquid species
        # --------------------------------------------------
        ns = self.num_species
        liquid_part = state_vec[idx:idx + ns]
        idx += ns

        if self.basis == 'mass_j':
            # state contains total mass of each species
            total_mass = liquid_part.sum()
            mass_frac = liquid_part / total_mass

            self.Liquid_1.updatePhase(
                mass_frac=mass_frac,
                mass=total_mass,
                solvent_pass=True
            )
        else:
            # state contains concentration
            self.Liquid_1.mass_conc = liquid_part.copy()

        # --------------------------------------------------
        # 3) Volume state (if dynamic)
        # --------------------------------------------------
        if 'vol' in self.solver_states:
            vol_state = state_vec[idx]
            idx += 1

            self.vol_tank = vol_state * self.vol_offset
            self.Liquid_1.vol = vol_state

            # Update geometry
            self.diam_tank = (4 / np.pi * self.vol_tank)**(1/3)
            self.area_base = np.pi / 4 * self.diam_tank**2

        # --------------------------------------------------
        # 4) Temperature states
        # --------------------------------------------------
        if 'temp_ht' in self.solver_states:
            temp_liq = state_vec[idx]
            temp_ht = state_vec[idx + 1]
            idx += 2

            self.Liquid_1.temp = temp_liq
            # do NOT overwrite Utility internals permanently
            # only store current value
            self._temp_ht_current = temp_ht

        elif 'temp' in self.solver_states:
            temp_liq = state_vec[idx]
            idx += 1
            self.Liquid_1.temp = temp_liq

        # --------------------------------------------------
        # 5) Ensure slurry references are consistent
        # --------------------------------------------------
        self.dx = self.Slurry.dx
        self.x_grid = self.Slurry.x_distrib


    def _fast_solve(self, runtime, time_grid,
                eval_sens, jac_v_prod,
                verbose, test,
                sundials_opts, any_event):

        
        
        states_init, merged_params = self._build_initial_state_and_params()

        if self._compiled_val_sens:
            self._problem.p = merged_params

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
                                    eval_sens, jac_v_prod,
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

        if self.vol_tank is None:
            self.vol_tank = self.Slurry.vol      

        if self.vol_tank is None:
            if isinstance(self, ReactiveSemibatchCrystallizer):
                time_vec = np.linspace(self.elapsed_time, final_time)
                vol_flow = self.get_inputs(time_vec)['Inlet']['vol_flow']

                self.vol_tank = trapezoidal_rule(time_vec, vol_flow)

            else:
                self.vol_tank = self.Slurry.vol

        self.sundials_opt = self._solver.get_options()
        if runtime is not None:
            final_time = runtime + self.elapsed_time
        if time_grid is not None:
            final_time = time_grid[-1]
        


       
        time, states= self._solver.simulate(final_time, ncp_list=time_grid)
        
        self.retrieve_results(time, states)

        # ---------- Organize sensitivity
        if eval_sens:
            sensit = []
            for elem in self._solver.p_sol:
                sens = np.array(elem)
                sens[0] = 0  # correct NaN's at t = 0 for sensitivities
                sensit.append(sens)

            self.sensit = sensit

            return time, states, sensit
        else:
            return time, states

    def flatten_states(self):
        out = flatten_states(self.profiles_runs)

        return out
    
        

    

class ReactiveMSMPR(_BaseReactiveCryst):
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