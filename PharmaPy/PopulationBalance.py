class PopulationBalance:
    def add_solver_state_variables(self,collection):
        pass
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


