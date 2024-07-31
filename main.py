from numpy import longdouble, pi, arange, fmod
import time
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm

from vpython import graph, gcurve, color


# Composite propellant rocket motor simulator
# John Nguyen <knightofthealtar64@gmail.com>
# MIT License

# Defines the mechanical properties of the pressure case.

class PressureCase:
    # Safety factors for each component
    CASE_PRESSURE_SAFE_FACTOR = 3
    PIN_PRESSURE_SAFE_FACTOR = 3

    def __init__(self, yield_strength: float, thickness: float, outer_radius: float, length: float, T_allowable: float, density: float, thermal_conductivity: float, specific_heat: float):
         self.T_allowable = T_allowable # deg F, case allowable temperature
         self.yield_strength = yield_strength # ksi, yield strength
         self.thickness = thickness # in, case thickness
         self.outer_radius = outer_radius # in, case outer radius
         self.inner_radius = outer_radius - thickness # in, case inner radius
         self.length = length # pressure case length
         self.density = density # lb / in^3, case density
         self.mass = density * pi * (self.outer_radius**2 - self.inner_radius**2) * self.length # lb, case mass
         self.thermal_conductivity = thermal_conductivity # BTU / in / s thermal conductivity
         self.specific_heat = specific_heat # BTU / lbm / R specific heat capacity
         self.inner_area = 2 * pi * self.inner_radius * self.length # in^2, internal area
         self.yield_pressure = 2 * yield_strength * thickness / (outer_radius*2) # ksi, yield pressure
         self.allowable_pressure = self.yield_pressure / PressureCase.CASE_PRESSURE_SAFE_FACTOR # ksi, case allowable pressure

    def insulate(self, insulation_thickness: float, insulation_decomposition_temp: float, insulation_decomposition_rate: float, insulation_density: float, insulation_thermal_conductivity: float, insulation_specific_heat: float):
         self.insulation_thickness = insulation_thickness # in, thickness of insulation
         self.insulation_decomposition_temp = insulation_decomposition_temp # degrees F, insulation decomposition temperature
         self.insulation_decomposition_rate = insulation_decomposition_rate # in / (thickness * deg F) insulation decomposition rate
         self.insulation_density = insulation_density # lb / in^3 insulation density
         self.insulation_thermal_conductivity = insulation_thermal_conductivity / 1055.06 * 5 / 9 / 39.37 / 3600 # BTU / (thickness * deg F * in) insulation conductivity
         self.insulation_specific_heat = insulation_specific_heat # BTU / (lb * deg F) insulation specific heat capacity

class TempMatrixElement:
    def __init__(self, inner_radius, outer_radius, length, density, temperature, thickness, specific_heat, conduction_htc):
        self.mass = pi * (outer_radius**2 - inner_radius**2) * length * density # lbm
        self.temperature = float(temperature) # R
        self.density = density # lbm / in^3
        self.inner_radius = inner_radius # in
        self.outer_radius = outer_radius # in
        self.length = length # in
        self.outer_area = 2 * pi * outer_radius * length # in^2
        self.inner_area = 2 * pi * inner_radius * length # in^2
        self.thickness = thickness # in
        self.specific_heat = specific_heat # BTU / R lbm
        self.conduction_htc = conduction_htc # Conductive heat transfer coefficient

    def erode(self, Td, dr, dt):
        if self.temperature - Td > 0 and self.thickness > 0:
            self.thickness -= (self.temperature - Td) * dr * dt 
            if self.thickness < 0:
                self.thickness = 0
            self.mass = pi * (self.outer_radius**2 - (self.outer_radius - self.thickness)**2) * self.length * self.density
            self.inner_area = 2 * pi * (self.outer_radius - self.thickness) / 2 * self.length

    def conductiveHeatTransfer(self, temperature_difference, dt):
        print(f"dT = {temperature_difference}, thick = {self.thickness}, specific_heat = {self.specific_heat}, mass = {self.mass}")
        dT = longdouble(self.conduction_htc / self.thickness * self.inner_area * temperature_difference / (self.specific_heat * self.mass) )
        self.temperature += dT * dt
        return dT * dt

def convectiveHeatTransfer(heat_transfer_coefficient, area, temperature_difference, specific_heat, mass):
        return longdouble(heat_transfer_coefficient * area * temperature_difference / (specific_heat * mass))

steel_sugar_case = dict(yield_strength=36, thickness=0.049, outer_radius=1.5, length=7.3, T_allowable=1000+459.7, density=0.283, thermal_conductivity=50/12/3600, specific_heat=0.117)
case = PressureCase(**steel_sugar_case)
test_insulation = dict(insulation_thickness=0.025, insulation_decomposition_temp=300+459.7, insulation_decomposition_rate=7e-6, insulation_density=0.0361, insulation_thermal_conductivity=0.144/12/3600, insulation_specific_heat=2.38e-4)
case.insulate(**test_insulation)

def heatFlowSim():
    t = 0
    dt = longdouble(1e-4)
    T = 525 # Rankine, ambient temperature
    ds = longdouble(3e-3) # diameter step
    Tc = 1450*9/5+32+459.7 # Rankine, chamber temperature.  This value of 1450 degrees Celcius corresponds to the adiabatic flame temperature of 65/35 KNO3/sucrose "rocket candy".

    Tmax = T # maximum temperature within the array

    convective_heat_transfer_coefficient = longdouble(1000 / 0.454 / 1055.06 * 5 / 9 / 39.37**2) # BTU / (in^2 * s * deg F) Convection heat transfer coefficient b/w fluid and casing/insulation

    elements = []
    insulation = None
    for idx, i in enumerate(arange(case.inner_radius, case.outer_radius, ds)):
        elements.append(TempMatrixElement(i, i+ds, case.length, case.density, T, ds, case.specific_heat, case.thermal_conductivity))
    if case.insulation_thickness:
        insulation = TempMatrixElement(case.inner_radius - case.insulation_thickness, case.inner_radius, case.length, case.insulation_density, T, case.insulation_thickness, case.insulation_specific_heat, case.insulation_thermal_conductivity)

    print(f"# of elements: {len(elements)}")

    tgraph = graph(title="Maximum temperature", xtitle="t", ytitle="T (Rankine)", fast=True)
    max_temp_curve = gcurve(label="Max casing temperature", color=color.red)
    insulation_temp_curve = gcurve(label="Insulation temperature", color=color.green)

    section_graph = graph(title="Slice temperature", xtitle="inches from inside of casing", ytitle="T (Rankine)", fast=False)
    slice_temperature_curve = gcurve(label="T at thickness", color=color.red)
    initial_temperature_curve = gcurve(label="initial T at thickness", color=color.blue)
    final_temperature_curve = gcurve(label="final T at thickness", color=color.green)

    insgraph = graph(title="Insulation thickness", xtitle="t", ytitle="in", fast=True)
    ins_s_curve = gcurve(label="Insulation thickness", color=color.blue)

    temp_record_z = []
    thick_record_x = []
    time_record_y = []

    max_idx = 0
    for idx, i in enumerate(elements):
        initial_temperature_curve.plot(i.thickness, i.temperature)
        max_idx = idx

    rt_0 = time.perf_counter()
    while elements[0].temperature <= case.T_allowable and t <= 1.5:
        disp = False
        if fmod(t,0.01) < 2*dt:
            disp = True
        for idx, i in enumerate(elements):
            if 0 < idx < max_idx:
                Tr = longdouble(elements[idx - 1].temperature)
                Td = longdouble(i.temperature)
                dT = i.conductiveHeatTransfer(Tr - Td, dt)
                elements[idx-1].temperature -= dT
            elif idx == 0 and not insulation: # innermost element receives heat from combustion gasses or the insulation, not a previous element
                Tr = longdouble(Tc) # high temperature
                Td = longdouble(i.temperature) # low temperature
                i.conductiveHeatTransfer(Tr - Td, dt)
            elif idx == 0:
                insulation.erode(case.insulation_decomposition_temp, case.insulation_decomposition_rate, dt)
                if disp:
                    ins_s_curve.plot(t, insulation.thickness)
                if insulation.thickness > 0:
                    insulation.temperature += convectiveHeatTransfer(convective_heat_transfer_coefficient, insulation.inner_area, Tc - insulation.temperature, case.insulation_specific_heat, insulation.mass)
                    if disp:
                        insulation_temp_curve.plot(t, insulation.temperature)
                    dT = i.conductiveHeatTransfer(insulation.temperature - i.temperature, dt)
                    insulation.temperature -= dT
                else:
                    # If insulation has been eroded away, heat transfers from convection
                    i.temperature += convectiveHeatTransfer(convective_heat_transfer_coefficient, i.inner_area, Tc - i.temperature, case.specific_heat, i.mass) * dt
            else: # outermost element pushes heat to the outside as well as receiving heat from within
                thermal_conductivity = longdouble(case.thermal_conductivity)
                Tr = longdouble(elements[idx - 1].temperature)
                Td = longdouble(i.temperature)
                i.conductiveHeatTransfer(elements[idx - 1].temperature - i.temperature, dt) * dt # heat gained from previous slice
                Tr = longdouble(i.temperature)
                Td = longdouble(T)
                i.conductiveHeatTransfer(Tr - Td, dt) # heat lost to surroundings
            temp_record_z.append(i.temperature)
            thick_record_x.append(i.inner_radius)
            time_record_y.append(t)
        if disp:
            max_temp_curve.plot(t, elements[0].temperature)
        print(f"t={t:.6f}, Tmax={elements[0].temperature:.3f}, realtime={time.perf_counter()-rt_0:.2f}")
        t+=dt

    for idx, i in enumerate(elements):
        final_temperature_curve.plot(i.thickness, i.temperature)
    # Termination of heat flow simulation

    thick_record_x = np.matrix(thick_record_x)
    print(thick_record_x)
    time_record_y = np.matrix(time_record_y)
    print(time_record_y)
    fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
    temp_record_z = np.matrix(temp_record_z)
    print(temp_record_z)
    plt.style.use('_mpl-gallery')
    ax.plot_surface(thick_record_x, time_record_y, temp_record_z, norm=False)
    ax.set(xticklabels=[],yticklabels=[],zticklabels=[])
    plt.show()

def main():
    heatFlowSim()

if __name__ == "__main__":
    main()
