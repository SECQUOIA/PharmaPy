from PharmaPy.Phases import LiquidPhase, SolidPhase
from PharmaPy.Streams import LiquidStream
from PharmaPy.Reactors_refactor import ContinuousReactor,BatchReactor
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
    mass=1000,
    mass_frac=[0.4,0.6,0,0,0],

)
print("starting vol:", liquid1.vol)
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
    mass_frac=[0.4,0.6,0,0,0]
)
rxns = ['A + B --> C', 'C + A --> D']
kvals_rxns = np.array([1e4,1e4])#, 1e2]) #psuedo instantaneous
ea_vals = np.array([1e2,1e2])#,1e4]) #psuedo no activation energy
Rkinetics = RxnKinetics(path=dpath,rxn_list=rxns, k_params=kvals_rxns,ea_params=ea_vals)
Utility = CoolingWater(mass_flow=100, temp_in=273.55)
vessel.Utility = Utility
vessel.RxnKinetics = Rkinetics
vessel.controller.target_volume=.001
vessel.Inlet = inlet


# -----------------------------
# Solve
# -----------------------------

vessel.solve_unit(
    runtime=3000
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

print(
    "Vessel Final composition:",
    vessel.Phases.mass_frac
)
print(
    "Vessel Final mass:",
    vessel.Phases.mass
)
print(
    "Vessel Final vol:",
    vessel.Phases.vol
)

if True:
    import matplotlib.pyplot as plt
    for mj,spec in zip(vessel.result.mass_j_liquid0.T,['A','B','C','D','Solvent']):
        if spec=='Solvent':continue
        plt.plot(vessel.result.time,mj,label=spec)
    plt.legend()
    plt.show()