import numpy as np
from CoolProp.CoolProp import PropsSI

# this is a new comment
# test

in2m    = 0.0254
Cv2Kv   = 0.862
pa2psi  = 1/6895
g       = 9.81 #m/s2

# Fluid System Component Framework
class Line:
    type = "line"

    def __init__(self, common_name, OD, ID, length, material):
        self.common_name = common_name
        self.OD = OD
        self.ID = ID
        self.length = length
        self.material = material

    def __repr__(self):
        return f"<Line | {self.common_name}>"


class Component:
    type = "component"

    def __init__(self, common_name, Cv, mass, material, ID):
        self.common_name = common_name
        self.Cv = Cv
        self.mass = mass
        self.material = material
        self.ID = ID

    def __repr__(self):
        return f"<Component | {self.common_name}>"


class MinorLoss:
    type = "minor_loss"

    def __init__(self, common_name, ID, k_factor):
        self.common_name = common_name
        self.ID = ID
        self.k_factor = k_factor

    def __repr__(self):
        return f"<MinorLoss | {self.common_name}>"

class Injector:
    type = "injector"

    def __init__(self, common_name, Cd, area, mass, material):
        self.common_name = common_name
        self.Cd = Cd
        self.area = area
        self.mass = mass
        self.material = material


    def __repr__(self):
        return f"<Injector | {self.common_name}>"

# Materials Class
class Material:
    def __init__(self, name, roughness, Cp, rho):
        self.name = name
        self.roughness = roughness  # [m]
        self.Cp = Cp                # [J/kg*K]
        self.rho = rho              # [kg/m^3]

    def __call__(self, prop):
        return getattr(self, prop)

    def __repr__(self):
        return f"<Material | {self.name}>"


# Material Properties dimensionless roughness, Cp [J/kg K], rho [kg/m3]

SS316  = Material("SS316",  roughness=4.6e-5, Cp=500,  rho=8000)
Brass  = Material("Brass",  roughness=1.5e-6, Cp=380,  rho=8500)
AL6061 = Material("AL6061", roughness=1.6e-6, Cp=896, rho=2700)

"""
Input Initial Conditions
"""

fluid       = "OXYGEN"
bottle_high = 2600 # psi
bottle_low  = 1400 # psi
p_0         = 650 / pa2psi # Pa
T_0         = 90  # K
m_dot       = 1.8534 #kg/s
time_soak   = 2 #seconds
t_ambient   = 298.15 # K

System_Arch = [
    Line("LINE1",OD=.75,ID=.694,length=6,material=SS316),
    Component("O-FLOW-BV6/NC",Cv=21,mass=1,material=SS316,ID=.694),
    Line("LINE2", OD=.75, ID=.694, length= 6,material= SS316),
    Component("O-ISOLATE-BV7", Cv=28, mass=1, material=SS316, ID=.694),
    Line("LINE3", OD=.75, ID=.694, length=6, material=SS316),
    Component("O-CHECK-CV2", Cv=10, mass=0.5, material=SS316, ID=.694),
    MinorLoss("90_ELBOW",ID=.694, k_factor=0.75),
    Line("LINE4", OD=.75, ID=.694, length=24, material=SS316),
    Injector("Injector (Ox side)", Cd=0.65,area=55, mass=1.1,material=AL6061)


]

rho_0   = PropsSI("Dmass", "T", T_0, "P", p_0, fluid)  # kg/m^3
sg_0    = rho_0/1000
q_0     = m_dot / rho_0 #m3/s
mu_0    = 300e-6 + (100e-6 - 300e-6) * (T_0 - 200) / 100  # kg/m*s

print(f"Starting analysis of {fluid}.\n"
      f"Initial Pressure:{p_0*pa2psi:7.3f} Psi\n"
      f"        Temperature: {T_0} k\n "
      f"        Density: {rho_0:10.3f} kg/m^3\n"
      f"        Specific Gravity: {sg_0}\n"
      f"        Dynamic Viscosity: {mu_0:10.7f} kg/m*s\n"
      f"        Mass Flow Rate: {m_dot:10.3f} kg/s\n"
      f"        Volumetric Flow Rate: {q_0:10.5f} m3/s\n")

p_current = p_0
t_current = T_0
rho_current = rho_0
q_current = q_0
mu_current = mu_0
sg_current = sg_0
dP_from_lines = 0
dP_from_component = 0
dP_from_minor_loss = 0
dP_from_injector = 0
dT_total = 0

# find ullage requirements
M_N2 = 0.028 #molar mass n2 in kg/mol
R = 8.314 # gas constant J/molK
m_dot_ullage = M_N2 * p_0/(R*T_0) * (m_dot/rho_0)



for item in System_Arch:
    dP = 0
    dT = 0

    # compute velocity and Reynolds number for items with a port diameter
    if hasattr(item, 'ID'):
        v = q_current / (np.pi/4 * (item.ID*in2m)**2)
        Re = rho_current * v * item.ID*in2m / mu_current

    # --- Pressure Drop ---
    if item.type == "line":
        # darcy weisbach friction factor from colebrook equation
        def colebrook(f, e):
            return 1 / f ** 0.5 + 2 * np.log10(e / (3.7 * item.ID * in2m) + 2.51 / (Re * f ** 0.5))

        a, b = 0.001, 0.1
        for _ in range(64): # bisection method
            mid = (a + b) / 2
            if colebrook(a,item.material("roughness")) * colebrook(mid,item.material("roughness")) < 0:
                b = mid
            else:
                a = mid
        f = (a + b) / 2
        dP = f * (item.length)/item.ID * (rho_current * v**2 )/2 #dP in Pa
        dP_from_lines = dP_from_lines + dP

    if item.type == "component":
        dP = sg_current * (q_current*3600/(item.Cv*Cv2Kv))**2 #dP in bar
        dP = dP *100*1000
        dP_from_component = dP_from_component + dP

    if item.type == "minor_loss":
        hm = item.k_factor * v**2 / (2 * g) #dP in m head
        dP = hm * rho_current * g #dP in Pa
        dP_from_minor_loss = dP_from_minor_loss + dP

    if item.type == "injector":
        Cd = item.Cd
        area = item.area/ 1e+6
        dP = 1/(2*rho_current) * (m_dot/(Cd * area))**2 #dp in Pa from standard orifice eq
        dP_from_injector = dP_from_injector + dP

    # --- Thermal Update (Level 2: Convective + Lumped Capacitance) ---
    try:
        k_fluid = PropsSI("L", "T", t_current, "P", p_current, fluid)
    except ValueError:
        k_fluid = 0.12 + (0.08 - 0.12) * (t_current - 200) / 100  # W/mK approx for liquid N2O
    cp_fluid = PropsSI("Cpmass", "T", t_current, "P", p_current, fluid)
    Pr = mu_current * cp_fluid / k_fluid

    if item.type == "line":
        Nu = 0.023 * Re**0.8 * Pr**0.4
        h_conv = Nu * k_fluid / (item.ID * in2m)
        A_contact = np.pi * item.ID * in2m * item.length * in2m
        m_metal = item.material.rho * (np.pi/4 * ((item.OD*in2m)**2 - (item.ID*in2m)**2)) * item.length * in2m
        tau = m_metal * item.material.Cp / (h_conv * A_contact)
        Q = m_metal * item.material.Cp * (t_ambient - t_current) * (1 - np.exp(-time_soak / tau))
        dT = Q / (m_dot * time_soak * cp_fluid)
        t_current += dT

    elif item.type == "component":
        Nu = 0.023 * Re**0.8 * Pr**0.4
        h_conv = Nu * k_fluid / (item.ID * in2m)
        effective_length = 3 * item.ID * in2m
        A_contact = np.pi * item.ID * in2m * effective_length
        tau = item.mass * item.material.Cp / (h_conv * A_contact)
        Q = item.mass * item.material.Cp * (t_ambient - t_current) * (1 - np.exp(-time_soak / tau))
        dT = Q / (m_dot * time_soak * cp_fluid)
        t_current += dT

    elif item.type == "injector":
        d_eff = np.sqrt(4 * item.area / 1e6 / np.pi)
        v_inj = q_current / (np.pi/4 * d_eff**2)
        Re_inj = rho_current * v_inj * d_eff / mu_current
        Nu = 0.023 * Re_inj**0.8 * Pr**0.4
        h_conv = Nu * k_fluid / d_eff
        effective_length = 5 * d_eff
        A_contact = np.pi * d_eff * effective_length
        tau = item.mass * item.material.Cp / (h_conv * A_contact)
        Q = item.mass * item.material.Cp * (t_ambient - t_current) * (1 - np.exp(-time_soak / tau))
        dT = Q / (m_dot * time_soak * cp_fluid)
        t_current += dT

    dT_total += dT

    # --- State Update ---
    p_current = p_current - dP
    rho_current = PropsSI("Dmass", "T", t_current, "P", p_current, fluid)  # kg/m^3
    sg_current = rho_current / 1000
    q_current = m_dot / rho_current  # m3/s
    mu_current = 300e-6 + (100e-6 - 300e-6) * (t_current - 200) / 100

    # --- Phase Check ---
    try:
        p_sat = PropsSI("P", "T", t_current, "Q", 0, fluid)
        if p_sat > p_current:
            print(f"  *** WARNING: Two-phase condition after {item.common_name}! "
                  f"P_sat={p_sat*pa2psi:.1f} psi > P_local={p_current*pa2psi:.1f} psi ***")
    except Exception:
        pass  # above critical temp, no saturation defined

    print(f"After {item.common_name}.\n"
          f"        Pressure: {p_current * pa2psi:7.3f} Psi\n"
          f"        Temperature: {t_current:.2f} K (dT = +{dT:.2f} K)\n "
          f"        Density: {rho_current:10.3f} kg/m^3\n"
          f"        Specific Gravity: {sg_current}\n"
          f"        Dynamic Viscosity: {mu_current:10.7f} kg/m*s\n"
          f"        Mass Flow Rate: {m_dot:10.3f} kg/s\n"
          f"        Volumetric Flow Rate: {q_current:10.5f} m3/s\n")

#end of per item loop
total_DP = (p_0 - p_current) * pa2psi

# Orifice sizing
gamma_n2 = 1.4
R_spec_n2 = 296.95
Cd_orifice = 0.85
choked_term = np.sqrt(gamma_n2 / (R_spec_n2 * T_0) * (2/2.4)**((2.4)/0.4))

# Size orifice at low bottle pressure for 1.5x ullage demand
A_orifice = m_dot_ullage * 1.5 / (Cd_orifice * (bottle_low/pa2psi + 101325) * choked_term)
d_orifice = np.sqrt(4 * A_orifice / np.pi)

# Max flow through that same orifice at high bottle pressure = worst case for burst disk
m_dot_max = Cd_orifice * A_orifice * (bottle_high/pa2psi + 101325) * choked_term
rho_n2_std = PropsSI("Dmass", "T", 293.15, "P", 101325, "N2")
max_vent_scfm = m_dot_max / rho_n2_std * 2118.88



print(f"Final Properties:\n"
          f"        Pressure: {p_current * pa2psi:.3f} Psi\n"
          f"        Temperature: {t_current:.2f} K (total dT = +{dT_total:.2f} K from {T_0:.2f} K)\n"
          f"        Density: {rho_current:.3f} kg/m^3\n"
          f"        Specific Gravity: {sg_current}\n"
          f"        Dynamic Viscosity: {mu_current:.7f} kg/m*s\n"
          f"        Mass Flow Rate: {m_dot:.3f} kg/s\n"
          f"        Volumetric Flow Rate: {q_current:.5f} m3/s\n"
          f"        dP from lines: {dP_from_lines*pa2psi:.3f} \n"
          f"        dP from minor losses: {dP_from_minor_loss*pa2psi:.3f} \n"
          f"        dP from components: {dP_from_component*pa2psi:.3f} \n"
          f"        dP from injector: {dP_from_injector*pa2psi:.3f} \n"
          f"        Percent dP at injector: {dP_from_injector/(p_0-dP_from_component-dP_from_lines-dP_from_minor_loss) *100:.2f}%\n"
          f"        Required Ullage Nitrogen: {m_dot_ullage:.3f} kg/s\n"
          f"        Required Ullage Nitrogen: {m_dot_ullage/PropsSI("Dmass", "T", 293.15, "P", 101000, "N2") * 2118.88:.3f} SCFM\n"
          f"        Required Flow Control Orifice Diameter: {d_orifice / in2m:.4f} in @ 1.5x max flow @ {bottle_low} psi bottle \n"
          f"        Required Maximum Venting Capability  {max_vent_scfm:.1f} SCFM @ {bottle_high} psi bottle")