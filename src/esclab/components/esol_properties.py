from esclab.simulate import *
from eeslib import fluid_properties as fp
from eeslib.functions import convert
import numpy as np

"""
Unified property database for incompressible fluids, including density, viscosity, and specific heat, where available.

Note that some water properties are provided here, and other fluids that are potentially compressible (e.g. air, argon, hydrogen) are included with ideal gas density models. For more detailed and comprehensive property data, users should import **eeslib** which provides an API for CoolProp.

Usage:
    from esclab.components.esol_properties import Incompressible
    props = Incompressible()
    props.viscosity('Air', T=300)           # Viscosity of air at 300 K
    props.specheat('Therminol Oil', T=350)  # Specific heat of Therminol Oil at 350 K
    props.enthalpy('Salt (60 NaNO3, 40 KNO3)', T=400, T_ref=300)  # Enthalpy change of salt from 300 K to 400 K

"""

class Incompressible:
    """Incompressible substance property database with density, viscosity, and specific heat, where available"""
    
    funcmap = {
        "Stainless_AISI316": {
            'density': (lambda T: 8349.38 - 0.341708*T - 0.0000865128*T*T),
            'specheat': (lambda T: 0.368455 + 0.000399548*T - 1.70558E-07*T*T),
        },
        "Water (liquid)": {
            # Use with caution - this is only valid near room temp/pressure
            'density': (lambda T: 1000),
            'specheat': (lambda T: 4.181),
        },
        "Air": {
            'viscosity': (lambda T: max([0.0000010765 + 7.15173E-08*T - 5.03525E-11*T*T + 2.02799E-14*T*T*T, 1.e-6])),
            'specheat': (lambda T: 1.03749 - 0.000305497*T + 7.49335E-07*T*T - 3.39363E-10*T*T*T),
        },
        "Salt (68 KCl, 32 MgCl2)": {
            'density': (lambda T: 1E-10*T*T*T - 3E-07*T*T - 0.4739*T + 2384.2),
            'viscosity': (lambda T: 0.0146*np.exp(2230./T)*0.001),
            'specheat': (lambda T: 1.156),
        },
        "Salt (8 NaF, 92 NaBF4)": {
            'density': (lambda T: 8E-09*T*T*T - 2E-05*T*T - 0.6867*T + 2438.5),
            'viscosity': (lambda T: 0.0877*np.exp(2240./T)*0.001),
            'specheat': (lambda T: 1.507),
        },
        "Salt (25 KF, 75 KBF4)": {
            'density': (lambda T: 2E-08*T*T*T - 6E-05*T*T - 0.7701*T + 2466.1),
            'viscosity': (lambda T: 0.0431*np.exp(3060./T)*0.001),
            'specheat': (lambda T: 1.306),
        },
        "Salt (31 RbF, 69 RbBF4)": {
            'density': (lambda T: -1E-08*T*T*T + 4E-05*T*T - 1.0836*T + 3242.6),
            'viscosity': (lambda T: 0.0009),
            'specheat': (lambda T: 9.127),
        },
        "Salt (46.5 LiF, 11.5 NaF, 42 KF)": {
            'density': (lambda T: -2E-09*T*T*T + 1E-05*T*T - 0.7427*T + 2734.7),
            'viscosity': (lambda T: 0.0400*np.exp(4170./T)*0.001),
            'specheat': (lambda T: 2.010),
        },
        "Salt (49 LiF, 29 NaF, 29 ZrF4)": {
            'density': (lambda T: -2E-11*T*T*T + 1E-07*T*T - 0.5172*T + 3674.3),
            'viscosity': (lambda T: 0.0069),
            'specheat': (lambda T: 1.239),
        },
        "Salt (58 KF, 42 ZrF4)": {
            'density': (lambda T: -6E-10*T*T*T + 4E-06*T*T - 0.8931*T + 3661.3),
            'viscosity': (lambda T: 0.0159*np.exp(3179./T)*0.001),
            'specheat': (lambda T: 1.051),
        },
        "Salt (58 LiCl, 42 RbCl)": {
            'density': (lambda T: -8E-10*T*T*T + 1E-06*T*T - 0.689*T + 2929.5),
            'viscosity': (lambda T: 0.0861*np.exp(2517./T)*0.001),
            'specheat': (lambda T: 8.918),
        },
        "Salt (58 NaCl, 42 MgCl2)": {
            'density': (lambda T: -5E-09*T*T*T + 2E-05*T*T - 0.5298*T + 2444.1),
            'viscosity': (lambda T: 0.0286*np.exp(1441./T)*0.001),
            'specheat': (lambda T: 1.080),
        },
        "Salt (59.5 LiCl, 40.5 KCl)": {
            'density': (lambda T: 1E-09*T*T*T - 5E-06*T*T - 0.864*T + 2112.6),
            'viscosity': (lambda T: 0.0861*np.exp(2517./T)*0.001),
            'specheat': (lambda T: 1.202),
        },
        "Salt (59.5 NaF, 40.5 ZrF4)": {
            'density': (lambda T: -5E-09*T*T*T + 2E-05*T*T - 0.9144*T + 3837.),
            'viscosity': (lambda T: 0.0767*np.exp(3977./T)*0.001),
            'specheat': (lambda T: 1.172),
        },
        "Salt (60 NaNO3, 40 KNO3)": {
            'density': (lambda T: max([-1E-07*T*T*T + 0.0002*T*T - 0.7875*T + 2299.4, 1000.0])),
            'viscosity': (lambda T: max([-1.473302E-10*(T-273.15)**3 + 2.279989E-07*(T-273.15)**2 - 1.199514E-04*(T-273.15) + 2.270616E-02, 0.0001])),
            'specheat': (lambda T: -1E-10*T*T*T + 2E-07*T*T + 5E-06*T + 1.4387),
        },
        "Nitrate Salt": {
            'density': (lambda T: max([2090 - 0.636 * (T-273.15), 1000.0])),
            'viscosity': (lambda T: max([(22.714 - 0.12 * (T-273.15) + 0.0002281 * (T-273.15)**2 - 0.0000001474 * (T-273.15)**3) / 1000, 1.e-6])),
            'specheat': (lambda T: (1443. + 0.172 * (T-273.15))/1000.),
        },
        "Caloria HT 43": {
            'density': (lambda T: max([885 - 0.6617 * (T-273.15) - 0.0001265 * (T-273.15)**2, 100.0])),
            'specheat': (lambda T: (3.88 * (T-273.15) + 1606.0)/1000.),
        },
        "HITEC XL Nitrate Salt,": {
            'density': (lambda T: max([2240 - 0.8266 * (T-273.15), 800.0])),
            'viscosity': (lambda T: 1372000 * (T-273.15)**-3.364),
            'specheat': (lambda T: max([1536 - 0.2624 * (T-273.15) - 0.0001139 * (T-273.15)**2, 1000.])/1000.),
        },
        "Therminol Oil": {
            'density': (lambda T: max([1074.0 - 0.6367 * (T-273.15) - 0.0007762 * (T-273.15)**2, 400.0])),
            'viscosity': (lambda T: 0.001 * (10**0.8703 * max(T-273.15, 20.)**(0.2877 + np.log10(max(T-273.15, 20.)**-0.3638)))),
            'specheat': (lambda T: 1.509 + 0.002496 * (T-273.15) + 0.0000007888 * (T-273.15)**2),
        },
        "HITEC Salt": {
            'density': (lambda T: max([2080 - 0.733 * (T-273.15), 1000.0])),
            'viscosity': (lambda T: max([0.00622 - 0.0000102 * (T-273.15), 1.e-6])),
            'specheat': (lambda T: (1560 - 0.0 * (T-273.15))/1000.),
        },
        "Dowtherm Q": {
            'density': (lambda T: max([-0.757332 * (T-273.15) + 980.787, 100.0])),
            'viscosity': (lambda T: 1 / (132.40658 + 4.36107 * (T-273.15) + 0.0781417 * (T-273.15)**2 - 0.00011035416 * (T-273.15)**3)),
            'specheat': (lambda T: (-0.00053943 * (T-273.15)**2 + 3.2028 * (T-273.15) + 1589.2)/1000.),
        },
        "Dowtherm RP": {
            'density': (lambda T: max([-0.000186495 * (T-273.15)**2 - 0.668337 * (T-273.15) + 1042.11, 200.0])),
            'viscosity': (lambda T: 1 / (4.523003 + 0.39156855 * (T-273.15) + 0.028604206 * (T-273.15)**2)),
            'specheat': (lambda T: (-0.0000031915 * (T-273.15)**2 + 2.977 * (T-273.15) + 1560.8)/1000.),
        },
        "HITEC XL": {
            'density': (lambda T: max([2240 - 0.8266 * (T-273.15), 800.0])),
            'viscosity': (lambda T: 1372000 * (T-273.15)**-3.364),
            'specheat': (lambda T: max([1536 - 0.2624 * (T-273.15) - 0.0001139 * (T-273.15)**2, 1000.])/1000.),
        },
        "T-91 Steel": {
            'density': (lambda T: -0.3289*(T-273.15) + 7742.5),
            'specheat': (lambda T: 0.0004*(T-273.15)**2 + 0.2473*(T-273.15) + 450.08),
        },
        "Therminol 66": {
            'density': (lambda T: -0.7146*(T-273.15) + 1024.8),
            'viscosity': (lambda T: (1.31959963 - 0.171204729*(T-273.15) + 0.0100351594*(T-273.15)**2 - 0.000313556341*(T-273.15)**3 + 0.0000053430666*(T-273.15)**4 - 4.66597650E-08*(T-273.15)**5 + 1.63046296E-10*(T-273.15)**6) if (T-273.15) < 80. else (0.0490075884 - 0.00120478233*(T-273.15) + 0.0000130162082*(T-273.15)**2 - 7.58913847E-08*(T-273.15)**3 + 2.47856063E-10*(T-273.15)**4 - 4.26872345E-13*(T-273.15)**5 + 3.01949160E-16*(T-273.15)**6)),
            'specheat': (lambda T: 0.0036*(T-273.15) + 1.4801),
        },
        "Therminol 59": {
            'density': (lambda T: -0.0003*(T-273.15)**2 - 0.6963*(T-273.15) + 988.44),
            'viscosity': (lambda T: (0.0137267822 - 0.000218740224*(T-273.15) + 0.0000759248815*(T-273.15)**2 - 0.00000473464744*(T-273.15)**3 - 1.97083667E-07*(T-273.15)**4 + 4.35487179E-09*(T-273.15)**5 + 2.40243056E-10*(T-273.15)**6) if (T-273.15) < 25. else (0.0114608807 - 0.000313431056*(T-273.15) + 0.00000416778121*(T-273.15)**2 - 3.04668508E-08*(T-273.15)**3 + 1.23719006E-10*(T-273.15)**4 - 2.60834697E-13*(T-273.15)**5 + 2.22227675E-16*(T-273.15)**6)),
            'specheat': (lambda T: 0.0033*(T-273.15) + 1.6132),
        },
        "Argon": {
            'density': (lambda T, P: max([P/(208.13*T), 1.e-10])),
            'viscosity': (lambda T: 4.4997e-6 + 6.38920E-08*T - 1.24550E-11*T*T),
            'specheat': (lambda T: 0.5203),
        },
        "Hydrogen": {
            'density': (lambda T, P: max([P/(4124.*T), 1.e-10])),
            'viscosity': (lambda T: 0.00000231 + 2.37842E-08*T - 5.73624E-12*T*T),
            'specheat': (lambda T: min([max([-45.4022 + 0.690156*T - 0.00327354*T*T + 0.00000817326*T**3 - 1.13234E-08*T**4 + 8.24995E-12*T**5 - 2.46804E-15*T**6, 11.3]), 14.7])),
        },
        "Dowtherm A": {
            # Density polynomial fit to Dow Chemical published data (kg/m³), T in K.
            # Fitted to: T_C=15→1061, 50→1031, 100→985, 150→936, 200→884, 250→830, 300→772, 350→712 kg/m³
            'density': (lambda T: max([1071.0 - 0.803*(T-273.15) - 0.000644*(T-273.15)**2, 100.0])),
            'viscosity': (lambda T: 0.786512*max(T-273.15, 20.)**-1.44263),
            'specheat': (lambda T: 1.47524 + 0.00368606*((T-273.15)-273) - 0.00000516458*((T-273.15)-273)**2 + 8.99399E-09*((T-273.15)-273)**3),
        },
    }

    def enthalpy(self, fluid: str, T: float, T_ref: float = 273.15, P: float = float('nan')):
        """Calculate specific enthalpy [J/kg] from temperature T [K] relative to T_ref [K] (default 0°C).
        Computed by integrating specific heat [J/(kg·K)] from T_ref to T.
        Only enthalpy differences are physically meaningful; the reference temperature cancels in all
        energy balance applications."""
        from scipy.integrate import quad
        assert fluid in self.funcmap.keys(), f"Fluid '{fluid}' not found in database"
        assert 'specheat' in self.funcmap[fluid], f"No specific heat data available for '{fluid}'"
        cp_func = self.funcmap[fluid]['specheat']
        h_kJ_per_kg, _ = quad(cp_func, T_ref, T)
        return h_kJ_per_kg  # Convert kJ/kg → J/kg

    def density(self, fluid: str, T: float, P: float = float('nan')):
        """Calculate density [kg/m^3] from temperature [K] and optional pressure [Pa].
        Argon and Hydrogen require pressure."""
        assert fluid in self.funcmap.keys(), f"Fluid '{fluid}' not found in database"
        assert 'density' in self.funcmap[fluid], f"No density data available for '{fluid}'"
        
        if fluid in ['Argon', 'Hydrogen']:
            assert P == P, "Pressure required for Argon and Hydrogen"
            return self.funcmap[fluid]['density'](T, P)
        else:
            return self.funcmap[fluid]['density'](T)

    def viscosity(self, fluid: str, T: float, P: float = float('nan')):
        """Calculate dynamic viscosity [Pa-s] from temperature [K] and optional pressure [Pa]."""
        assert fluid in self.funcmap.keys(), f"Fluid '{fluid}' not found in database"
        assert 'viscosity' in self.funcmap[fluid], f"No viscosity data available for '{fluid}'"
        return self.funcmap[fluid]['viscosity'](T)

    def specheat(self, fluid: str, T: float, P: float = float('nan')):
        """Calculate specific heat [J/kg-K] from temperature [K] and optional pressure [Pa]."""
        assert fluid in self.funcmap.keys(), f"Fluid '{fluid}' not found in database"
        assert 'specheat' in self.funcmap[fluid], f"No specific heat data available for '{fluid}'"
        return self.funcmap[fluid]['specheat'](T) * 1000. # Convert kJ/kg-K → J/kg-K
    
    def print_fluids(self):
        """Print list of fluids in database"""
        for fluid in self.funcmap.keys():
            available_properties = self.funcmap[fluid].keys()
            abbrev_properties = {'density': 'ρ', 'viscosity': 'μ', 'specheat': 'c'}
            abbrev_available = [abbrev_properties[prop] for prop in available_properties]
            print(fluid + " → (" + ", ".join(abbrev_available) + ")")

if __name__ == "__main__":
    props = Incompressible()
    props.print_fluids()