from PharmaPy.Phases import LiquidPhase, SolidPhase
from PharmaPy.Streams import LiquidStream
from PharmaPy.Reactors_refactor import ContinuousReactor
from PharmaPy.IntegratorBackends import AssimuloBackend
from PharmaPy.Kinetics import RxnKinetics
from PharmaPy.Utilities import CoolingWater


import numpy as np


dpath = r"C:\Users\zhillma\OneDrivePZH\Documents\Documents\_Grad_School\mypharmadev\PharmaPy\tests\Flowsheet\data\compound_database.json"


# -----------------------------
# Reactor
# -----------------------------

integrator = AssimuloBackend()

vessel = ContinuousReactor(
    integrator=integrator,
    h_conv=10000,
    diam=.01,

)


# -----------------------------
# Initial phases
# -----------------------------

liquid1 = LiquidPhase(
    dpath,
    mass=1,
    mass_frac=[0.3,0.7,0,0,0],

)

solid1 = SolidPhase(
    dpath,
    mass=0,
    mass_frac=[0,0,1,0,0]
)


vessel.Phases = liquid1


# -----------------------------
# Feed
# -----------------------------

inlet = LiquidStream(
    dpath,
    mass_flow=0.02,
    mass_frac=[0.3,0.7,0,0,0]
)
rxns = ['A + B --> C']
kvals_rxns = np.array([1e7])#, 1e2]) #psuedo instantaneous
ea_vals = np.array([1e2])#,1e4]) #psuedo no activation energy
Rkinetics = RxnKinetics(path=dpath,rxn_list=rxns, k_params=kvals_rxns,ea_params=ea_vals)
Utility = CoolingWater(mass_flow=100, temp_in=273.55)
vessel.Utility = Utility
vessel.RxnKinetics = Rkinetics
vessel.Inlet = inlet


# -----------------------------
# Solve
# -----------------------------

vessel.solve_unit(
    runtime=300
)


# -----------------------------
# Check outlet
# -----------------------------

print("done")

print(
    "Outlet flow:",
    vessel.Outlet.mass_flow
)

print(
    "Outlet composition:",
    vessel.Outlet.mass_frac
)
print("outlet temp:",
      vessel.Outlet.temp)