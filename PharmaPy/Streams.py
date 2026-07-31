from PharmaPy.Phases import BasePhase
from PharmaPy.Interpolation import NewtonInterpolation
from PharmaPy.Results import DynamicResult
from PharmaPy.DataClasses import StateVariable,StateKey

from scipy.interpolate import CubicSpline
import numpy as np

class BaseStream(BasePhase):
    def __init__(self, path_thermo=None, temp=298.15, pressure=101325, mass_flow=None, mole_flow=None, vol_flow=None, mass_j_flow=None, mass_frac=None, mole_frac=None, mass_conc=None, mole_conc=None, name_solv=None, check_input=True, verbose=True,
                 controls=None,args_control=None,num_interpolation_points=3,verbose=True, check_input=True, **kwargs):
        super().__init__(path_thermo, temp, pressure, mass_flow, mole_flow, vol_flow, mass_j_flow, mass_frac, mole_frac, mass_conc, mole_conc, name_solv, check_input, verbose, **kwargs)
        self._DynamicInlet = None
        self.controllable = ('mass_flow', 'mole_flow', 'vol_flow', 'temp') #Deprecation Warning
        self.time_upstream = None
        self.bipartite = None
        ## Deprecated control logic for compatibility
        # Controls
        if controls is None:
            controls = {}
        else:
            if args_control is None:
                args_control = {key: () for key in controls.keys()}

            update_dict = {}
            for key, fun in controls.items():
                update_dict[key] = fun(0, *args_control[key])

            self.updatePhase(**update_dict)

        self.controls = controls
        self.args_control = args_control

        self.num_interpolation_points = num_interpolation_points
        
    # ====================================================
    # Aliases for flows
    # =====================================================

    @property
    def mass_flow(self):
        return self.mass
    @mass_flow.setter
    def mass_flow(self,value):
        if value is not None:
            self.mass = value
    @property
    def mole_flow(self):
        return self.moles
    @mole_flow.setter
    def mole_flow(self,value):
        if value is not None:
            self.moles = value
    @property
    def vol_flow(self):
        return self.vol
    @vol_flow.setter
    def vol_flow(self,value):
        if value is not None:
            self.vol = value
     @property
    def mass_j_flow(self):
        return self.mass_j
    @mass_j_flow.setter
    def mass_j_flow(self,value):
        if value is not None:
            self.mass_j = value


    # ========================================
    # Change key names for applicable methods
    # ==========================================

    def updatePhase(self,mass_flow=None,mole_flow=None,vol_flow=None,mass_j_flow=None, mole_conc=None, mass_conc=None, mass_frac=None, mole_frac=None, mass_j=None, mass=None, vol=None, moles=None, temp=None, pres=None, **kwargs):
        if any([all([mass_flow,mass]),all([mole_flow,moles]),all(vol_flow,vol), all(mass_j_flow,mass_j)]):
            raise RuntimeError("Only a qunatity or flow can be specified, not both")
        return super().updatePhase(mole_conc, mass_conc, mass_frac, mole_frac, mass_j, mass, vol, moles, temp, pres, **kwargs)

    @property
    def state_dict(self):
        state_dict = super().state_dict()
        state_dict['mass_flow'] = state_dict['mass']
        state_dict['vol_flow'] = state_dict['vol']
        state_dict['mole_flow'] = state_dict['moles']
        state_dict['mass_j_flow'] = state_dict['mass_j']
        return state_dict


    # ========================================
    # Input logic
    # ==========================================
    @property
    def DynamicInlet(self):
        return self._DynamicInlet

    @DynamicInlet.setter
    def DynamicInlet(self, dynamic_object):
        dynamic_object.controllable = self.controllable
        dynamic_object.parent_instance = self

        self._DynamicInlet = dynamic_object

    def InterpolateInputs(self, time):
        if isinstance(time, (float, int)):
            # Assume steady state for extrapolation
            time = min(time, self.time_upstream[-1])

            y_interpol = Interpolation(self.time_upstream, self.y_inlet,
                                        time,
                                        num_points=self.num_interpolation_points)
        else:
            interpol = CubicSpline(self.time_upstream, self.y_inlet)
            flags_interpol = time > self.time_upstream[-1]

            if any(flags_interpol):
                time_interpol = time[~flags_interpol]
                y_interp = interpol(time_interpol)

                y_extrapol = np.tile(y_interp[-1],
                                        (sum(flags_interpol), 1))
                y_interpol = np.vstack((y_interp, y_extrapol))
            else:
                y_interpol = interpol(time)

        return y_interpol
    def evaluate_inputs(self, time):
        if self.DynamicInlet is None:
            inputs = {}
            for attr in self.controllable:
                inputs[attr] = getattr(self, attr)

        else:
            inputs = self.DynamicInlet.evaluate_inputs(time)

        return inputs
    def add_stream_state_variables(self, collection):
        collection.add(
            StateVariable(
                name="vol_flow",
                dim=1,
                units="m3/s",
                stream="inlet"
            )
        )

        collection.add(
            StateVariable(
                name="temp",
                dim=1,
                units="K",
                stream="inlet"
            )
        )

        collection.add(
            StateVariable(
                name="mole_conc",
                dim=len(self.name_species),
                units="kmol/m3",
                stream="inlet"
            )
        )


class LiquidStream(BaseStream):
    def __init__(self, path_thermo=None, temp=298.15, pressure=101325, mass_flow=None, mole_flow=None, vol_flow=None, mass_j_flow=None, mass_frac=None, mole_frac=None, mass_conc=None, mole_conc=None, name_solv=None, check_input=True, verbose=True, controls=None, args_control=None, num_interpolation_points=3, verbose=True, check_input=True, **kwargs):
        self.phase_family='liquid'
        super().__init__(path_thermo, temp, pressure, mass_flow, mole_flow, vol_flow, mass_j_flow, mass_frac, mole_frac, mass_conc, mole_conc, name_solv, check_input, verbose, controls, args_control, num_interpolation_points, verbose, check_input, **kwargs)
        self.cp_liq = np.atleast_2d(self.cp_liq)
        self.p_vap = np.atleast_2d(self.p_vap)
class VaporPhase(BaseStream):
    def __init__(self, path_thermo=None, temp=298.15, pressure=101325, mass_flow=None, mole_flow=None, vol_flow=None, mass_j_flow=None, mass_frac=None, mole_frac=None, mass_conc=None, mole_conc=None, name_solv=None, check_input=True, verbose=True, controls=None, args_control=None, num_interpolation_points=3, verbose=True, check_input=True, **kwargs):
        self.phase_family = 'vapor'
        super().__init__(path_thermo, temp, pressure, mass_flow, mole_flow, vol_flow, mass_j_flow, mass_frac, mole_frac, mass_conc, mole_conc, name_solv, check_input, verbose, controls, args_control, num_interpolation_points, verbose, check_input, **kwargs)

class SolidPhase(BaseStream):
    def __init__(self, path_thermo=None, temp=298.15, pressure=101325, mass_flow=None, mole_flow=None, vol_flow=None, mass_j_flow=None, mass_frac=None, mole_frac=None, mass_conc=None, mole_conc=None, name_solv=None, check_input=True, verbose=True, controls=None, args_control=None, num_interpolation_points=3, verbose=True, check_input=True, **kwargs):
        self.phase_family='solid'
        super().__init__(path_thermo, temp, pressure, mass_flow, mole_flow, vol_flow, mass_j_flow, mass_frac, mole_frac, mass_conc, mole_conc, name_solv, check_input, verbose, controls, args_control, num_interpolation_points, verbose, check_input, **kwargs)

    #============================================================
    # Compatibility API: Deprecation Warning
    # ===========================================================
    def getPorosity(self, distrib=None, diam_filter=1, AR=None,
                    sphericity=None):

        if distrib is None:
            distrib = self.distrib
            mom_zero = self.moments[0]
            mom_one = self.moments[1]
        else:
            mom_zero, mom_one = self.getMoments(mom_num=(0, 1))

        # mom_one *= 1e-6  # m
        x_dist = self.x_distrib * 1e-6  # m

        if AR is None:
            AR = 2

        if sphericity is None:
            sphericity = 0.7

        # Yu, Zou et al (1996) and Yu,Zou, Stnadish (1996) model
        kv = 0.524  # Volumetric shape coefficient
        ks = 3.142  # Surface shape coefficient

        del_x_dist = np.diff(x_dist)
        node_x_dist = (x_dist[:-1] + x_dist[1:]) / 2
        node_CSD = (distrib[:-1] + distrib[1:]) / 2

        # Volume of crystals in each bin
        vol_cry = node_CSD * del_x_dist * (kv * node_x_dist**3)
        frac_vol_cry = vol_cry / (np.sum(vol_cry) + eps)

        vol_particle = kv * node_x_dist**3
        d_part_sphere = (6 * vol_particle / np.pi)**(1/3)
        d_part_equiv_pack = d_part_sphere / (sphericity**2.785 *
                                                np.exp(2.946 * (1 - sphericity)))

        # Initial porosity
        D_mean = mom_one/(mom_zero + eps)
        E_0_Jeschar = 0.375 + 0.34 * D_mean/diam_filter  # average porosity of packing of uniform sized spheres [-]

        initial_porosity = E_0_Jeschar

        V = 1/(1 - initial_porosity) * np.ones_like(node_x_dist)  # Specific Volume for initial porosity

        # Evaluate specific volume using modified linear packing model
        num_x = len(node_x_dist)
        V_T_node = np.zeros(num_x)

        for i in range(num_x):

            r = d_part_equiv_pack[:i] / d_part_equiv_pack[i]
            g_r = (1 - r)**2 + 0.4*r*(1 - r)**3.7
            V_large_j = V[:i] - (V[:i] - 1) * g_r - V[i]
            sum_V_large_term = sum(V_large_j * frac_vol_cry[:i])

            r_inv = r = d_part_equiv_pack[i] / d_part_equiv_pack[i + 1:]
            f_r = (1 - r_inv)**3.3 + 2.8*r_inv*(1 - r_inv)**2.7
            V_small_j = V[i + 1:] * (1 - f_r) - V[i]
            sum_V_small_term = sum(V_small_j * frac_vol_cry[i + 1:])

            V_T_node[i] = V[i] + sum_V_large_term + sum_V_small_term

        V_T = max(V_T_node)

        porosity = 1 - 1/V_T

        return porosity