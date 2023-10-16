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


     def __init__(self, sy: float, s: float, od: float, lg: float, ta: float, rho: float, khr: float, cp: float):
         self.ta = ta # deg F, case allowable temperature
         self.sy = sy # ksi, yield strength
         self.s = s # in, case thickness
         self.od = od # in, case OD
         self.id = od - 2 * s # in, case ID
         self.lg = lg # pressure case length
         self.rho = rho # lb / in^3, case density
         self.mass = rho * pi * ((self.od/2)**2 - (self.id/2)**2) * self.lg # lb, case mass
         self.khr = khr # BTU / (hr * deg F * ft) thermal conductivity
         self.k = khr * 1 / 3600 * 1 / 12 # BTU / (s * deg F * in) thermal conductivity
         self.cp = cp # BTU / (lb * deg F) specific heat capacity
         self.ai = 2 * pi * (self.id/2) * self.lg # in^2, internal area
         self.Py = 2 * sy * s / od # ksi, yield pressure
         self.Pa = self.Py / PressureCase.CASE_PRESSURE_SAFE_FACTOR # ksi, case allowable pressure

class TempMatrixElement:
    def __init__(self, di, od, lg, rho, Te, s):
        self.mass = pi * ((od/2) ** 2 - (di / 2) ** 2) * rho
        self.Te = float(Te)
        self.ao = 2 * pi * od / 2 * lg
        self.ai = 2 * pi * di / 2 * lg
        self.s = s

def tempFlow(ki, s, A, dTi, cp, mass):
    return longdouble(ki/s*A*dTi/(cp*mass))

t = 0
dt = longdouble(1e-6)
T = 90 # deg F, ambient temperature
ds = longdouble(1e-3) # diameter step
Tc = 4000 # deg F, chamber temperature

Tmax = T # maximum temperature within the array

case = PressureCase(45, 0.049, 1.5, 7.3, 900, 2.83e-1, 20.8, 1.17e-1)
k_flow = longdouble(3.397e-7) # Overall heat transfer coefficient b/w fluid and casing


elements = []
for idx, i in enumerate(arange(case.id, case.od, ds)):
    elements.append(TempMatrixElement(i, i+ds, case.lg, case.rho, T, idx*ds))

print(f"# of elements: {len(elements)}")

tgraph = graph(title="Maximum temperature", xtitle="t", ytitle="T (deg F)", fast=True)
max_curve = gcurve(label="Tmax", color=color.red)

sgraph = graph(title="Slice temperature", xtitle="s", ytitle="T (deg F)", fast=True)
scurve = gcurve(label="T at s", color=color.red)
ic = gcurve(label="initial T at s", color=color.blue)
fcurve = gcurve(label="final T at s", color=color.green)

max_idx = 0
for idx, i in enumerate(elements):
    ic.plot(i.s, i.Te)
    max_idx = idx

rt_0 = time.perf_counter()
while elements[0].Te <= case.ta and t <= 1:
    for idx, i in enumerate(elements):
        if idx == 0: # innermost element receives heat from combustion gasses, not a previous element
            k = longdouble(k_flow)
            Tr = longdouble(Tc) # high temperature
            Td = longdouble(i.Te) # low temperature
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
    #max_curve.plot(t, elements[0].Te)
    print(f"t={t}, Tmax={elements[0].Te}, realtime={time.perf_counter()-rt_0}")
    t+=dt

for idx, i in enumerate(elements):
    fcurve.plot(i.s, i.Te)