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

    def __init__(self, sy: float, s: float, rm: float, lg: float, ta: float, rho: float, khr: float, cp: float):
         self.ta = ta # deg F, case allowable temperature
         self.sy = sy # ksi, yield strength
         self.s = s # in, case thickness
         self.rm = rm # in, case outer radius
         self.ri = rm - s # in, case inner radius
         self.lg = lg # pressure case length
         self.rho = rho # lb / in^3, case density
         self.mass = rho * pi * (self.rm**2 - self.ri**2) * self.lg # lb, case mass
         self.khr = khr # BTU / (hr * deg F * ft) thermal conductivity
         self.k = khr * 1 / 3600 * 1 / 12 # BTU / (s * deg F * in) thermal conductivity
         self.cp = cp # BTU / (lb * deg F) specific heat capacity
         self.ai = 2 * pi * self.ri * self.lg # in^2, internal area
         self.Py = 2 * sy * s / (rm*2) # ksi, yield pressure
         self.Pa = self.Py / PressureCase.CASE_PRESSURE_SAFE_FACTOR # ksi, case allowable pressure

    def insulate(self, in_s: float, in_td: float, in_dr: float, in_rho: float, in_k: float, in_cp: float):
         self.in_s = in_s # in, thickness of insulation
         self.in_td = in_td # degrees F, insulation decomposition temperature
         self.in_dr = in_dr # in / (s * deg F) insulation decomposition rate
         self.in_rho = in_rho # lb / in^3 insulation density
         self.in_k = in_k * 1 / 3600 * 1 / 12 # BTU / (s * deg F * in) insulation conductivity
         self.in_cp = in_cp # BTU / (lb * deg F) insulation specific heat capacity

class TempMatrixElement:
    def __init__(self, ri, rm, lg, rho, Te, s):
        self.mass = pi * (rm**2 - ri**2) * lg * rho
        self.Te = float(Te)
        self.rho = rho
        self.ri = ri
        self.rm = rm
        self.lg = lg
        self.ao = 2 * pi * rm * lg
        self.ai = 2 * pi * ri * lg
        self.s = s

    def erode(self, Td, dr, dt):
        if self.Te - Td > 0 and self.s > 0:
            self.s -= (self.Te - Td) * dr * dt 
            if self.s < 0:
                self.s = 0
            self.mass = pi * (self.rm**2 - (self.rm - self.s)**2) * self.lg * self.rho
            self.ai = 2 * pi * (self.rm - self.s) / 2 * self.lg

def tempFlow(ki, s, A, dTi, cp, mass):
    return longdouble(ki/s*A*dTi/(cp*mass))

case = PressureCase(36, 0.049, 1.5, 7.3, 900, 2.83e-1, 31.2, 1.17e-1)
case.insulate(0.010, 300, 7e-6, 0.0361, 0.5779e-7, 2.38e-4)

def heatFlowSim():
    t = 0
    dt = longdouble(1e-5)
    T = 90 # deg F, ambient temperature
    ds = longdouble(1e-3) # diameter step
    Tc = 2900 # deg F, chamber temperature

    Tmax = T # maximum temperature within the array

    k_flow = longdouble(3.397e-7) # BTU / (in * s * deg F) Convection heat transfer coefficient b/w fluid and casing/insulation

    elements = []
    insulation = None
    for idx, i in enumerate(arange(case.ri, case.rm, ds)):
        elements.append(TempMatrixElement(i, i+ds, case.lg, case.rho, T, idx*ds))
    if case.in_s:
        insulation = TempMatrixElement(case.ri - case.in_s * 2, case.ri, case.lg, case.in_rho, T, case.in_s)

    print(f"# of elements: {len(elements)}")

    tgraph = graph(title="Maximum temperature", xtitle="t", ytitle="T (deg F)", fast=False)
    max_curve = gcurve(label="Tmax", color=color.red, markers=True)

    sgraph = graph(title="Slice temperature", xtitle="inches from inside of casing", ytitle="T (deg F)", fast=False)
    scurve = gcurve(label="T at s", color=color.red)
    ic = gcurve(label="initial T at s", color=color.blue)
    fcurve = gcurve(label="final T at s", color=color.green)

    insgraph = graph(title="Insulation properties", xtitle="t", ytitle="T deg F, s in", fast=False)
    ins_s_curve = gcurve(label="Insulation thickness", color=color.blue)
    ins_T_curve = gcurve(label="Insulation temperature", color=color.red)

    max_idx = 0
    for idx, i in enumerate(elements):
        ic.plot(i.s, i.Te)
        max_idx = idx

    rt_0 = time.perf_counter()
    while elements[0].Te <= case.ta and t <= 2:
        disp = False
        if fmod(t,0.01) < 2*dt:
            disp = True
        for idx, i in enumerate(elements):
            if idx == 0 and not insulation: # innermost element receives heat from combustion gasses or the insulation, not a previous element
                k = longdouble(k_flow)
                Tr = longdouble(Tc) # high temperature
                Td = longdouble(i.Te) # low temperature
                i.Te += tempFlow(k, ds, i.ai, Tr - Td, case.cp, i.mass) * dt
            elif idx == 0 and insulation:
                insulation.erode(case.in_td, case.in_dr, dt)
                if disp:
                    ins_s_curve.plot(t, insulation.s)
                if insulation.s > 0:
                    insulation.Te += tempFlow(case.in_k, insulation.s, insulation.ai, Tc - insulation.Te, case.in_cp, insulation.mass)
                    if disp:
                        ins_T_curve.plot(t, insulation.Te)
                    k = longdouble(case.in_k)
                    Tr = longdouble(Tc)
                    Td = longdouble(i.Te)
                    dT = tempFlow(k, insulation.s, case.ai, Tr - Td, case.in_cp, insulation.mass) * dt
                    insulation.Te -= dT
                    i.Te += dT
                else:
                    k = longdouble(k_flow)
                    Tr = longdouble(Tc)  # high temperature
                    Td = longdouble(i.Te)  # low temperature
                    i.Te += tempFlow(k, ds, i.ai, Tr - Td, case.cp, i.mass) * dt
            elif idx == max_idx: # outermost element pushes heat to the outside as well as receiving heat from within
                k = longdouble(case.k)
                Tr = longdouble(elements[idx - 1].Te)
                Td = longdouble(i.Te)
                i.Te += tempFlow(k, ds, i.ai, Tr - Td, case.cp, i.mass) * dt # heat gained from previous slice
                k = longdouble(case.k * 1e-3)
                Tr = longdouble(i.Te)
                Td = longdouble(T)
                i.Te -= tempFlow(k, ds, i.ao, Tr - Td, case.cp, i.mass) * dt # heat lost to surroundings
            else:
                k = longdouble(case.k)
                Tr = longdouble(elements[idx - 1].Te)
                Td = longdouble(i.Te)
                dT = tempFlow(k, ds, i.ai, Tr - Td, case.cp, i.mass) * dt
                i.Te += dT
                elements[idx-1].Te -= dT
        if disp:
            max_curve.plot(t, elements[0].Te)
        print(f"t={t:.6f}, Tmax={elements[0].Te:.3f}, realtime={time.perf_counter()-rt_0:.2f}")
        t+=dt

    for idx, i in enumerate(elements):
        fcurve.plot(i.s, i.Te)
    # Termination of heat flow simulation

def main():
    heatFlowSim()

if __name__ == "__main__":
    main()
