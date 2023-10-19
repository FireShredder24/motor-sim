from numpy import longdouble, pi, arange, fmod
import time

from vpython import graph, gcurve, color


# Composite propellant rocket motor simulator
# John Nguyen <knightofthealtar64@gmail.com>
# MIT License

# Defines the mechanical properties of the pressure case.

class PressureCase:
    # Safety factors for each component
    CASE_PRESSURE_SAFE_FACTOR = 3
    PIN_PRESSURE_SAFE_FACTOR = 3

    def __init__(self, yield_strength: float, thickness: float, outer_radius: float, length: float, T_allowable: float, density: float, thermal_conductivity_btu_hr_ft: float, specific_heat: float):
         self.T_allowable = T_allowable # deg F, case allowable temperature
         self.yield_strength = yield_strength # ksi, yield strength
         self.thickness = thickness # in, case thickness
         self.outer_radius = outer_radius # in, case outer radius
         self.inner_radius = outer_radius - thickness # in, case inner radius
         self.length = length # pressure case length
         self.density = density # lb / in^3, case density
         self.mass = density * pi * (self.outer_radius**2 - self.inner_radius**2) * self.length # lb, case mass
         self.thermal_conductivity_btu_hr_ft = thermal_conductivity_btu_hr_ft # BTU / (hr * deg F * ft) thermal conductivity
         self.thermal_conductivity = thermal_conductivity_btu_hr_ft * 1 / 3600 * 1 / 12 # BTU / (thickness * deg F * in) thermal conductivity
         self.specific_heat = specific_heat # BTU / (lb * deg F) specific heat capacity
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
    def __init__(self, inner_radius, outer_radius, length, density, temperature, thickness):
        self.mass = pi * (outer_radius**2 - inner_radius**2) * length * density
        self.temperature = float(temperature)
        self.density = density
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.length = length
        self.outer_area = 2 * pi * outer_radius * length
        self.inner_area = 2 * pi * inner_radius * length
        self.thickness = thickness

    def erode(self, Td, dr, dt):
        if self.temperature - Td > 0 and self.thickness > 0:
            self.thickness -= (self.temperature - Td) * dr * dt 
            if self.thickness < 0:
                self.thickness = 0
            self.mass = pi * (self.outer_radius**2 - (self.outer_radius - self.thickness)**2) * self.length * self.density
            self.inner_area = 2 * pi * (self.outer_radius - self.thickness) / 2 * self.length

def conductiveHeatTransfer(heat_transfer_coefficient, thickness, area, temperature_difference, specific_heat, mass):
    return longdouble(heat_transfer_coefficient / thickness * area * temperature_difference / (specific_heat * mass) )

def convectiveHeatTransfer(heat_transfer_coefficient, area, temperature_difference, specific_heat, mass):
    return longdouble(heat_transfer_coefficient * area * temperature_difference / (specific_heat * mass))

case = PressureCase(36, 0.049, 1.5, 7.3, 1000 + 459.7, 2.83e-1, 23.116, 1.17e-1)
case.insulate(0.025, 300 + 459.7, 7e-6, 0.0361, 1, 2.38e-4)

def heatFlowSim():
    t = 0
    dt = longdouble(1e-5)
    T = 9/5*20 + 32 + 459.7# Rankine, ambient temperature
    ds = longdouble(1e-3) # diameter step
    Tc = 1450*9/5 + 32 + 459.7 # Rankine, chamber temperature

    Tmax = T # maximum temperature within the array

    convective_heat_transfer_coefficient = longdouble(1000 / 0.454 / 1055.06 * 5 / 9 / 39.37**2) # BTU / (in^2 * s * deg F) Convection heat transfer coefficient b/w fluid and casing/insulation

    elements = []
    insulation = None
    for idx, i in enumerate(arange(case.inner_radius, case.outer_radius, ds)):
        elements.append(TempMatrixElement(i, i+ds, case.length, case.density, T, idx*ds))
    if case.insulation_thickness:
        insulation = TempMatrixElement(case.inner_radius - case.insulation_thickness, case.inner_radius, case.length, case.insulation_density, T, case.insulation_thickness)

    print(f"# of elements: {len(elements)}")

    tgraph = graph(title="Maximum temperature", xtitle="t", ytitle="T (deg F)", fast=True)
    max_temp_curve = gcurve(label="Max casing temperature", color=color.red)
    insulation_temp_curve = gcurve(label="Insulation temperature", color=color.green)

    section_graph = graph(title="Slice temperature", xtitle="inches from inside of casing", ytitle="T (deg F)", fast=False)
    slice_temperature_curve = gcurve(label="T at thickness", color=color.red)
    initial_temperature_curve = gcurve(label="initial T at thickness", color=color.blue)
    final_temperature_curve = gcurve(label="final T at thickness", color=color.green)

    insgraph = graph(title="Insulation properties", xtitle="t", ytitle="T deg F, thickness in", fast=True)
    ins_s_curve = gcurve(label="Insulation thickness", color=color.blue)

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
                thermal_conductivity = longdouble(case.thermal_conductivity)
                Tr = longdouble(elements[idx - 1].temperature)
                Td = longdouble(i.temperature)
                dT = conductiveHeatTransfer(thermal_conductivity, ds, i.inner_area, Tr - Td, case.specific_heat, i.mass) * dt
                i.temperature += dT
                elements[idx-1].temperature -= dT
            elif idx == 0 and not insulation: # innermost element receives heat from combustion gasses or the insulation, not a previous element
                thermal_conductivity = longdouble(convective_heat_transfer_coefficient)
                Tr = longdouble(Tc) # high temperature
                Td = longdouble(i.temperature) # low temperature
                i.temperature += conductiveHeatTransfer(thermal_conductivity, ds, i.inner_area, Tr - Td, case.specific_heat, i.mass) * dt
            elif idx == 0:
                insulation.erode(case.insulation_decomposition_temp, case.insulation_decomposition_rate, dt)
                if disp:
                    ins_s_curve.plot(t, insulation.thickness)
                if insulation.thickness > 0:
                    insulation.temperature += convectiveHeatTransfer(case.insulation_thermal_conductivity, insulation.inner_area, Tc - insulation.temperature, case.insulation_specific_heat, insulation.mass)
                    if disp:
                        insulation_temp_curve.plot(t, insulation.temperature)
                    thermal_conductivity = longdouble(case.insulation_thermal_conductivity)
                    dT = conductiveHeatTransfer(thermal_conductivity, insulation.thickness, case.inner_area, Tc - i.temperature, case.insulation_specific_heat, insulation.mass) * dt
                    insulation.temperature -= dT
                    i.temperature += dT
                else:
                    # If insulation has been eroded away, heat transfers from convection
                    i.temperature += convectiveHeatTransfer(convective_heat_transfer_coefficient, i.inner_area, Tc - i.temperature, case.specific_heat, i.mass) * dt
            else: # outermost element pushes heat to the outside as well as receiving heat from within
                thermal_conductivity = longdouble(case.thermal_conductivity)
                Tr = longdouble(elements[idx - 1].temperature)
                Td = longdouble(i.temperature)
                i.temperature += conductiveHeatTransfer(thermal_conductivity, ds, i.inner_area, elements[idx - 1].temperature - i.temperature, case.specific_heat, i.mass) * dt # heat gained from previous slice
                thermal_conductivity = longdouble(case.thermal_conductivity * 1e-3)
                Tr = longdouble(i.temperature)
                Td = longdouble(T)
                i.temperature -= conductiveHeatTransfer(thermal_conductivity, ds, i.outer_area, Tr - Td, case.specific_heat, i.mass) * dt # heat lost to surroundings
        if disp:
            max_temp_curve.plot(t, elements[0].temperature)
        print(f"t={t:.6f}, Tmax={elements[0].temperature:.3f}, realtime={time.perf_counter()-rt_0:.2f}")
        t+=dt

    for idx, i in enumerate(elements):
        final_temperature_curve.plot(i.thickness, i.temperature)
    # Termination of heat flow simulation

def main():
    heatFlowSim()

if __name__ == "__main__":
    main()
