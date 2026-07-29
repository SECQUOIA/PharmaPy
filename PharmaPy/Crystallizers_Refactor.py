from PharmaPy.MultiPhaseVessel import MultiPhaseVessel
from PharmaPy.Mechanisms import *
from PharmaPy.DataClasses import *

class _BaseCrystallizer(MultiPhaseVessel):
    def __init__(self,target_comp):
        super().__init__()
        if isinstance(target_comp, str):
            target_comp = [target_comp]
        self.target_comp = target_comp
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
        if self.target_comp is not None:
            self.target_ind = []
            for tc in self.target_comp:
                name_bool = [name == tc for name in self.name_species] #TODO check that it selects correctly
                self.target_ind.append(np.where(name_bool)[0][0])
        for i,ck in enumerate(self._CrystKinetics):
            weights = np.zeros(self.Phases.num_species)
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

    def _post_CrystKinetics_setter(self):
        "Place holder in case future children need special behavior"
        pass
    def configure_solver(self):
        #Assimulo option, does nothing if not using assimulo backend
        self.integrator.solver.linear_solver = "SPGMR"