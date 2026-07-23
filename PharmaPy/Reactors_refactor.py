from PharmaPy.MultiPhaseVessel import MultiPhaseVessel

class _BaseReactors(MultiPhaseVessel):
    def __init__(self, target_comp, temp_ref, isothermal, reset_states, controls, h_conv, state_events, population_balance_method, scale, adiabatic, jac_type, basis='mass_j', ht_mode="jacket", diam=0, area_base=0):
        super().__init__(target_comp, temp_ref, isothermal, reset_states, controls, h_conv, state_events, population_balance_method, scale, adiabatic, jac_type, basis, ht_mode, diam, area_base):
    
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
            assert len(self.Liquids)>=i-1, "The number of reaction kinetics must match or be less than the number of liquid phases or you must specify reaction_regions manually"
            region = IntraPhaseProcess(phase=PhaseRef("liquid",i),
                                    mechanism=ReactionMechanism(rk))
            processes.append(region)
        self.intraphase_processes = processes

    def _post_RxnKinetics_setter(self,RKs):
        "Place holder in case future children need special behavior"
        pass
