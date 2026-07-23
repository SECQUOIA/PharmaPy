import numpy as np
import copy
    
class Mechanism:
    """
    Base class for all physical mechanisms.

    Responsibilities
    ----------------
    • Store internal solver/output state variables.
    • Compute material source terms.
    • Cache intermediate quantities in aux.
    • Compute energy source terms from aux.

    Contract
    --------
    material_rate, aux = get_material_rate(...)

    qdot = get_heat_generation(
        aux,
        completed_state,
        time,
    )

    The returned aux object is opaque to the vessel and should contain all
    intermediate quantities required by get_heat_generation().
    """

    solver_states = ()
    output_states = ()

    def update_state(self,state):

        for variable in self.solver_states:

            if variable.name in state:
                setattr(
                    self,
                    variable.name,
                    state[variable.name]
                )

    def add_solver_state_variables(
            self,
            state_collection,
            overwrite=False
        ):

        for state in self.solver_states:

            state_collection.add(
                copy.deepcopy(state),
                overwrite
            )


    def add_output_state_variables(
            self,
            state_collection,
            overwrite=False
        ):

        for state in self.output_states:

            state_collection.add(
                copy.deepcopy(state),
                overwrite
            )
    def get_heat_generation(self,aux:dict,completed_state:dict,time:np.ndarray)->float:
        raise NotImplementedError
        return float
class TransferMechanism(Mechanism):
    """
    Computes material exchanged between two phases.

    Returns
    -------
    transfer_rate : ndarray
        Species mass transfer rates.

    aux : dict
        Cached information required by get_heat_generation().
    """
    def get_material_transfer_rate(
            source_phase,
            sink_phase,
            connection,
            completed_state,
            time,
            params
        ):

        raise NotImplementedError
    
    
class ReactionMechanism(Mechanism):
    """
    Computes intraphase material generation or consumption.

    Returns
    -------
    species_mass_rates : ndarray
        Species generation/consumption rates.

    aux : dict
        Cached reaction information required by
        get_heat_generation().
    """
    def __init__(self, kinetics):
        self.kinetics=kinetics
    def get_material_rate(
            self,
            process,
            time,
            completed_state,
            params,
            u_input,
            molarity_molPerLiter=True
        ):
        "aux must have process field that stores process"
        mole_adjust = 1000 if molarity_molPerLiter else 1
        rk = self.kinetics
        phase = process.phase
        

        temp = phase.temp

        mask = np.array([
            species in rk.partic_species
            for species in phase.name_species
        ])

        conc = phase.mole_conc

        deltah_rxn = None

        if rk.keq_params is not None:
            deltah_rxn = phase.getHeatOfRxn(
                rk.stoich_matrix,
                temp,
                mask,
                rk.delta_hrxn,
                rk.tref_hrxn
            )

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

        aux = {
                "phase": phase,
                "rxn_rates": reaction_rates,
                "process":process
                }

        return species_mass_rates, aux
    def get_heat_generation(self, aux, completed_state, time,molarity_molPerLiter=True):
        mole_adjust = 1000 if molarity_molPerLiter else 1
        phase = aux['phase']
        temp = phase.temp
        rk = self.kinetics
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
        q =  -(deltah_rxn * aux["rxn_rates"]).sum()* phase.vol *mole_adjust # molarity here is mol/L, but J/kg is expected for internal consistency
        return q
    
class DirectTransfer(TransferMechanism):

    def get_material_transfer_rate(
            self,
            source_phase,
            sink_phase,
            kinetics,
            **kwargs):

        return kinetics.get_rate(
            source_phase,
            sink_phase,
            **kwargs
        )
    
class PopulationBalance(TransferMechanism):
    
    def add_output_state_variables(self,outputs):
        outputs.add(
                StateVariable(
                    name="supersat",
                    dim=len(self.target_comp),
                    units="kg/m3",
                    state_type="alg"
                )
            )

        outputs.add(
            StateVariable(
                name="solubility",
                dim=len(self.target_comp),
                units="kg/m3",
                state_type="alg"
            )
        )
    def get_material_transfer_rate(self, source, sink, connection, completed_state, time):
        # return super().get_material_rate(source, sink, connection, completed_state, time)
        if connection.species_weights is None:

                species_rate = (transfer_rate* source.mass_frac)

            else:

                species_rate = (transfer_rate* np.asarray(connection.species_weights))

class MomentsPopulationBalance(PopulationBalance):
    
    def add_output_state_variables(self, outputs):
        super().add_output_state_variables(outputs)
        outputs.add(
            StateVariable(
                name="mu_n",
                dim=4,
                index=list(range(4)),
                units="m**n",
                state_type="alg"
            )
        )


class FVMPopulationBalance(PopulationBalance):
    def add_solver_state_variables(self, collection):
        collection.add(
            StateVariable(
                name="distrib",
                ...
            )
        )
    def add_output_state_variables(self, outputs):
        super().add_output_state_variables(outputs)

        outputs.add(
            StateVariable(
                name="mu_n",
                dim=4,
                index=list(range(4)),
                units="m**n",
                state_type="alg"
            )
        )

        outputs.add(
            StateVariable(
                name="vol_distrib",
                dim=self.num_distr,
                index=list(range(self.num_distr)),
                units="m3/m3",
                state_type="alg"
            )
        )
    def get_material_transfer_rate(self, source, sink, connection, completed_state, time):
        transfer_rate = (
            crystal_growth_mass_rate
            + nucleation_mass_rate
            - dissolution_mass_rate
        )
        return transfer_rate




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