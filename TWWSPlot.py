import numpy as np
import matplotlib.pyplot as plt

# thermophysical properties at flight point
# Altitude: 4000 m
# Mach Number: 0.3
# Temperature: 262.15 K
# Pressure: 61640.79 Pa
# Freestream Velocity: 97.37 m/s or 350.54 km/h
# Air Density: 0.82 kg/m^3
# Dynamic Viscosity: 1.66e-05 kg/(m*s)
# Reynolds Number: 4.48e+06

# defining the S_wet function for the UCAV
def S_wet(W0, c = 0.8565, d = 0.5423):
    S_wet = 10**c * W0**d
    return S_wet

# defining a function for Swet based on the Regression from Gotten et. al
def Swet_UAV(MTOM, A = 2.5535, B = 0.6677):
    Swet = A * MTOM**B
    return Swet


# defining the Wing loading requirement for takeoff using Rocket Assisted Takeoff (RATO)
rho = 0.82 # air density in kg/m^3
# rho = 0.82 * 0.00194032 # air density in slugs/ft^3
rho = 0.82 * 0.062428 # air density in lb/ft^3
Clmax_takeoff = 2 # assumed maximum lift coefficient for takeoff
RATO_gvalue = 40 # gs in MKS for RATO
RATO_gvalue = 40 * 3.28084 # gs in ft/s^2 for RATO
t = 5 # time in seconds for RATO
WS_takeoff = 0.5 * rho * Clmax_takeoff * (RATO_gvalue)**2
print(f"\n\nThe Wing Loading Requirement for Takeoff using RATO is {WS_takeoff:.2f} lb/ft^2")

# assuming Cf = 0.0035 for the UCAV (Skin friction coefficient)
Cf = 0.011

# assuming Sref = 40 ft^2 for the UCAV (Reference wing area) (fig 4.4 of the Metabook)
Sref = 80 # ft^2

# calculating the CD0 for the UCAV
# swet = S_wet(900)
# CD0 = Cf * S_wet(900) / Sref
swet = Swet_UAV(900)
CD0 = Cf * Swet_UAV(900) / Sref
print(f"\n\nThe CD0 for the UCAV is {CD0:.4f}")

# defining the Oswald efficiency factors for the UCAV during clean and takeoff configurations 
# (as landing characteristics not necessary due to the final attack phase)
e_clean = 0.85
e_takeoff = 0.8

# defining the del_CD0 for the takeoff and clean configurations
delCD0_clean = 0
delCD0_takeoff = 0.01


# defining the Thrust to weight ratio requirement for climb segment
ks = 1.2 # ratio of flight speed to stall speed
AR = 10 # aspect ratio
CLmax_cruise = 1.2 # assumed maximum lift coefficient for cruise
G = 0.04 # gradient of the climb segment
TW_climb = ks**2 * (CD0 + delCD0_clean) / CLmax_cruise + (CLmax_cruise / (ks**2 * np.pi * AR * e_clean)) + G
print(f"\n\nThe Thrust to Weight Ratio Requirement for the Climb Segment is {TW_climb:.2f}")

# defining the Thrust to weight ratio function (wrt Wing Loading) for the ceiling
# q_ceiling = 0.5 * rho * 97.37**2
q_ceiling = 0.5 * rho * 319.48 ** 2 # lbf/ft^2
print(f"\n\nq_ceiling = {q_ceiling:.2f}")
G_ceiling = 0.001 # gradient of the ceiling segment
def TW_ceiling(WS):
    TW = q_ceiling * (CD0 + delCD0_clean) / (WS) + (WS / (np.pi * AR * e_clean * q_ceiling)) + G_ceiling
    return TW

# defining a Thrust to weight ratio function (wrt Wing Loading) for the cruise segment
def TW_cruise(WS):
    TW = q_ceiling * (CD0 + delCD0_clean) / (WS) + (WS / (np.pi * AR * e_clean * q_ceiling))
    return TW

# defining a Thrust to weight ratio function (wrt Wing Loading) for the maneuver segment
def TW_maneuver(WS, n = 5): # n is the load factor
    # debug print
    print("Debug print for maneuver segment:")
    partA = q_ceiling * (CD0 + delCD0_clean) / (WS)
    partB = (WS * n**2 / (np.pi * AR * e_clean * q_ceiling))

    TW = q_ceiling * (CD0 + delCD0_clean) / (WS) + (WS * n**2 / (np.pi * AR * e_clean * q_ceiling))
    return TW

# plotting the TW vs WS graph for the preliminary design of the UCAV
WS = np.linspace(1, 1000, 200)
TW_ceiling_values = TW_ceiling(WS)
TW_cruise_values = TW_cruise(WS)
TW_maneuver_values = TW_maneuver(WS)

fig, ax = plt.subplots(dpi = 300) #dpi = 300
ax.plot(WS, TW_ceiling_values, label = "Ceiling", color = "blue")
ax.plot(WS, TW_cruise_values, label = "Cruise", color = "orange")
ax.plot(WS, TW_maneuver_values, label = "Maneuver", color = "purple")

# plot a vertical line at the Wing Loading Requirement for Takeoff using RATO
ax.axvline(x = WS_takeoff, color = "red", label = "Takeoff")

# plot a horizontal line at the Thrust to Weight Ratio Requirement for the Climb Segment
ax.axhline(y = TW_climb, color = "green", label = "Climb")

# color the region to the left of the WS_takeoff and above the TW_maneuver
ax.fill_between(WS, TW_maneuver_values, 10, where = WS < WS_takeoff, color = "blue", alpha = 0.5)
# ax.set_xlim(0, 700)
ax.set_ylim(0, 1)
# mark a red dot at (1300, 0.7) for the design point
ax.plot(850, 0.45, "ro", label = "Design Point (850, 0.45)", markersize = 3)

ax.set_xlabel("Wing Loading (W/S)(lb/ft^2)")
ax.set_ylabel("Thrust to Weight Ratio (T/W)")
ax.set_title("T/W vs W/S Plot for the UCAV")
ax.legend()
ax.grid()
plt.savefig("TWvsWS_4.png")
plt.show()