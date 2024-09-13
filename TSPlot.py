import numpy as np
import matplotlib.pyplot as plt


def S_wet(W0, c = 0.8565, d = 0.5423):
    if W0 < 0:
        raise ValueError("W0 must be non-negative")    
    S_wet = 10**c * W0**d
    return S_wet

# defining a function for Swet based on the Regression from Gotten et. al
def Swet_UAV(MTOM, A = 2.5535, B = 0.6677):
    Swet = A * MTOM**B
    return Swet

# assuming Sref = 40 ft^2 for the UCAV (Reference wing area) (fig 4.4 of the Metabook)
# Sref = 40 # ft^2
rho = 0.82 # air density in kg/m^3
# rho = 0.82 * 0.00194032 # air density in slugs/ft^3
rho = 0.82 * 0.062428 # air density in lb/ft^3
# defining the Thrust to weight ratio function (wrt Wing Loading) for the ceiling
# q_ceiling = 0.5 * rho * 97.37**2
q_ceiling = 0.5 * rho * 319.48 ** 2 # lbf/ft^2
G_ceiling = 0.001 # gradient of the ceiling segment
delCD0_clean = 0
e_clean = 0.85
# calculating the CD0 for the UCAV
# CD0 = Cf * S_wet(793) / Sref

# defining a function to get the CD0 for a given S (the definition in metabook sounded dumb so let's just use this)
def getCD0(S, W0):
    if W0 < 0:
        raise ValueError("W0 must be non-negative")
    # assuming Cf = 0.011 for the UCAV (Skin friction coefficient)
    Cf = 0.011
    # Swet = S_wet(W0)
    Swet = Swet_UAV(W0)
    CD0 = Cf * Swet / S
    return CD0

# defining the TSFC at sea level
c_SL = 0.3 # tsfc in lb/(lbf.hr)
c_cruise = 0.5 # tsfc in lb/(lbf.hr)
c_loiter = 0.4 # tsfc in lb/(lbf.hr)

AR_global = 14

#defining the detailed fuel fraction function
def get_WfW0(S, T, W0):
    print("\033[1m\t\tWfW0 debug print:\033[0m")
    W1W0 = 1 - c_SL * (1/4) * 0.05 * T / W0 # fuel burned from running the engine for 15min at 5% max thrust
    print(f"\t\t\tW1W0 (Starting) = {W1W0:.2f}")
    W1 = W1W0 * W0
    W2W1 = 1 - c_SL * (1/60) * T / W1 # fuel burned from running the engine for 1 minutes at max thrust
    print(f"\t\t\tW2W1 (Takeoff) = {W2W1:.2f}")
    W3W2 = 0.985 # climb segment from historical data
    # computing the CD0 for the given S
    CD0 = getCD0(S, W0)
    print(f"\t\t\t\tCD0 = {CD0:.2f}")
    e = 0.85 # for the clean configuration
    AR = AR_global # assumed aspect ratio
    k = (1 / (np.pi * e * AR))
    CL = np.sqrt(CD0 / k) # assuming the aircraft is flying at the optimal CL for max L/D
    print(f"\t\t\t\tCL = {CL:.2f}")
    LoverD = 0.94 * CL / (CD0 + k * CL**2)
    print(f"\t\t\t\tLoverD = {LoverD:.2f}")
    V = 220 # cruise speed in mph
    R = 500 # range in miles
    W4W3 = np.exp(-R * c_cruise / (V * LoverD)) # cruise segment
    print(f"\t\t\tW4W3 (Cruise) = {W4W3:.2f}")
    E = 7 # endurance in hours
    W5W4 = np.exp(-E*c_loiter / LoverD) # loiter segment
    print(f"\t\t\tW5W4 (Loiter) = {W5W4:.2f}")
    W6W5 = 0.99 #descent
    W7W6 = 0.995 #landing
    W6W0 = W1W0 * W2W1 * W3W2 * W4W3 * W5W4 * W6W5 * W7W6
    WfW0 = 1 - W6W0
    return WfW0

def get_Wwarhead(W0, A = 0.9422, B = 0.2289): # from M. Voskulijil research paper
    if W0 < 0:
        raise ValueError("W0 must be non-negative")
    W_warhead = B * W0**A
    return W_warhead

# defining a function for the empty weight fraction regression
def emptyweightfraction(W0, A =1.53, C = -0.16):
    if W0 < 0:
        raise ValueError("W0 must be non-negative")
    W0 = W0 / 2.20462 # converting lb to kg
    WeW0 = (A) * (W0)**(C)
    return WeW0

# Wing density for the UCAV
wingdensity = 10 # lb/ft^2
enginedensity = 1.3

Wguess_global = 1500



# defining the iterating function to estimate TOGW as a function of T and S
def getTOGW(T, S, Wguess = Wguess_global, TWratio_designpoint = 0.45, WSratio_designpoint = 250):
    W0 = Wguess
    tol = 1e-6
    delta = 2 * tol
    i = 0
    T_designpoint = W0 * TWratio_designpoint
    S_designpoint = W0 / WSratio_designpoint
    print(f"\tDesign Point: T = {T_designpoint:.2f} lb, S = {S_designpoint:.2f} ft^2")
    print(f"\tCurrent Point: T = {T:.2f} lb, S = {S:.2f} ft^2")
    while delta > tol and i < 100:
        print(f"\t\tgetTOGW Iteration {i}: W0 = {W0:.2f} lb")
        WeW0 = emptyweightfraction(W0, 0.699, -0.051)
        print(f"\t\t\tWeW0 = {WeW0:.2f}")
        We = WeW0 * W0
        We = We + wingdensity * (S - S_designpoint)
        We = We + enginedensity * (T - T_designpoint)
        WfW0 = get_WfW0(S, T, W0)
        print(f"\t\t\tWfW0 = {WfW0:.2f}")
        # W_warhead = get_Wwarhead(W0)
        # print(f"\t\t\tW_warhead = {W_warhead:.2f}")
        W_warhead = 200 # lb
        W0_new = W_warhead / (1 - WeW0 - WfW0)
        delta = abs((W0_new - W0) / (W0_new))
        W0 = W0_new
        i += 1
        # print(f"\nIteration {i}: W0 = {W0:.2f} lb")
    return W0


def TW_ceiling(WS, S, W0, AR = AR_global, e_clean = 0.85, delCD0_clean = 0):
    # AR = 5 # assumed aspect ratio
    TW = q_ceiling * (getCD0(S, W0) + delCD0_clean) / (WS) + (WS / (np.pi * AR * e_clean * q_ceiling)) + G_ceiling
    return TW

# defining a Thrust to weight ratio function (wrt Wing Loading) for the cruise segment
def TW_cruise(WS, S, W0, AR = AR_global, e_clean = 0.85, delCD0_clean = 0):
    partA = q_ceiling * (getCD0(S, W0) + delCD0_clean) / (WS)
    partB = (WS / (np.pi * AR * e_clean * q_ceiling))
    TW = q_ceiling * (getCD0(S, W0) + delCD0_clean) / (WS) + (WS / (np.pi * AR * e_clean * q_ceiling))
    return TW

# defining a Thrust to weight ratio function (wrt Wing Loading) for the maneuver segment
def TW_maneuver(WS, S, W0, n = 5, AR = AR_global, e_clean = 0.85, delCD0_clean = 0):
    partA = q_ceiling * (getCD0(S, W0) + delCD0_clean) / (WS)
    partB = (WS / (np.pi * AR * e_clean * q_ceiling))
    TW = q_ceiling * (getCD0(S, W0) + delCD0_clean) / (WS) + (WS * n**2 / (np.pi * AR * e_clean * q_ceiling))
    return TW

# defining the TW requirement obtained for the climb segment
# TW_climb = 0.39


# # defining the Thrust to weight ratio requirement for climb segment
# ks = 1.2 # ratio of flight speed to stall speed
# AR = 10 # aspect ratio
# CLmax_cruise = 1.2 # assumed maximum lift coefficient for cruise
# G = 0.04 # gradient of the climb segment
# TW_climb = ks**2 * (CD0 + delCD0_clean) / CLmax_cruise + (CLmax_cruise / (ks**2 * np.pi * AR * e_clean)) + G

def TW_climb(WS, S, W0, AR = AR_global, e_clean = 0.85, delCD0_clean = 0):
    ks = 1.2 # ratio of flight speed to stall speed
    CLmax_cruise = 1.2 # assumed maximum lift coefficient for cruise
    G = 0.04 # gradient of the climb segment
    TW = ks**2 * (getCD0(S, W0) + delCD0_clean) / CLmax_cruise + (CLmax_cruise / (ks**2 * np.pi * AR * e_clean)) + G
    return TW

# defining the value obtained for WS_takeoff
WS_takeoff = 685.00 # N/m^2

# # trying out the getTOGW function
# T = 2000
# S = 100
# TOGW = getTOGW(T, S)
# print(f"\nThe estimated TOGW for T = {T} lbf and S = {S} ft^2 is {TOGW:.2f} lb, ie {TOGW/2.20462:.2f} kg")
# print(f"\n The wahread weight is {get_Wwarhead(TOGW):.2f} lb, ie {get_Wwarhead(TOGW)/2.20462:.2f} kg")

# defining a function for : iterating over values of S to obtain the constrint lines for the TW vs WS plot
def getTvals(Tguess = 1200):
    # Svals = np.linspace(5, 10, 11)
    Svals = np.linspace(45, 100, 11)
    # Svals = np.linspace(10, 50, 11)
    Tvals_ceiling = np.zeros_like(Svals)
    Tvals_cruise = np.zeros_like(Svals)
    Tvals_maneuver = np.zeros_like(Svals)
    Tvals_climb = np.zeros_like(Svals)
    Wguess = Wguess_global

    # loop for the ceiling constraint
    for i in range(len(Svals)):
        print(f"\nCeiling Constraint Loop: Iteration {i}")
        S = Svals[i]
        Tvals_ceiling[i] = Tguess
        tol = 1e-6
        converged = False
        W_last = Wguess_global
        while not converged:
            W = getTOGW(Tvals_ceiling[i], S, W_last)
            WS = W / S
            TW_ceiling_val = TW_ceiling(WS, S, W)
            print(f"\tTW_ceiling_val = {TW_ceiling_val:.2f}")
            Tnew = W * TW_ceiling_val
            W_last = W
            if abs(Tnew - Tvals_ceiling[i]) < tol:
                converged = True
            else:
                print(f"\t getTvals debug print: Value of T did not converge. \n\tTvals_ceiling[{i}] = {Tvals_ceiling[i]} \n\tTnew = {Tnew}")
            Tvals_ceiling[i] = Tnew
        print(f"Tvals_ceiling[{i}] = {Tvals_ceiling[i]}")

    # loop for the cruise constraint
    for i in range(len(Svals)):
        print(f"\nCruise Constraint Loop: Iteration {i}")
        S = Svals[i]
        Tvals_cruise[i] = Tguess
        tol = 1e-6
        converged = False
        while not converged:
            W = getTOGW(Tvals_cruise[i], S)
            WS = W / S
            TW_cruise_val = TW_cruise(WS, S, W)
            Tnew = W * TW_cruise_val
            if abs(Tnew - Tvals_cruise[i]) < tol:
                converged = True
            Tvals_cruise[i] = Tnew

    # loop for the maneuver constraint
    for i in range(len(Svals)):
        print(f"\nManeuver Constraint Loop: Iteration {i}")
        S = Svals[i]
        Tvals_maneuver[i] = Tguess
        tol = 1e-6
        converged = False
        while not converged:
            W = getTOGW(Tvals_maneuver[i], S)
            WS = W / S
            TW_maneuver_val = TW_maneuver(WS, S, W)
            Tnew = W * TW_maneuver_val
            if abs(Tnew - Tvals_maneuver[i]) < tol:
                converged = True
            Tvals_maneuver[i] = Tnew

    # loop for the climb constraint
    for i in range(len(Svals)):
        print(f"\nClimb Constraint Loop: Iteration {i}")
        S = Svals[i]
        Tvals_climb[i] = Tguess
        tol = 1e-6
        converged = False
        W_last = Wguess
        while not converged:
            W = getTOGW(Tvals_climb[i], S, W)
            WS = W / S
            TW_climb_val = TW_climb(WS, S, W)
            Tnew = W * TW_climb_val
            # W_last = W
            if abs(Tnew - Tvals_climb[i]) < tol:
                converged = True
            Tvals_climb[i] = Tnew

    return Svals, Tvals_ceiling, Tvals_cruise, Tvals_maneuver, Tvals_climb

# defining a function for : iterating over values of T to obtain the constrint lines for the TW vs WS plot
def getSvals(Sguess = 50):
    # Tvals = np.linspace(500, 2000, 10)
    Tvals = np.linspace(2000, 500, 10)
    Svals_takeoff = np.zeros_like(Tvals)

    # loop for the takeoff constraint
    for i in range(len(Tvals)):
        print(f"\nTakeoff Constraint Loop: Iteration {i}")
        T = Tvals[i]
        Svals_takeoff[i] = Sguess
        tol = 1e-6
        converged = False
        W_last = Wguess_global + 1000
        while not converged:
            W = getTOGW(T, Svals_takeoff[i], W_last)
            TWratio = T / W
            WS_takeoff_val = WS_takeoff
            Snew = W / WS_takeoff_val
            W_last = W
            if abs(Snew - Svals_takeoff[i]) < tol:
                converged = True
            else:
                print(f"\t getSvals debug print: Value of S did not converge. \n\tSvals_takeoff[{i}] = {Svals_takeoff[i]} \n\tSnew = {Snew}")
            Svals_takeoff[i] = Snew

    return Tvals, Svals_takeoff

# plotting the T vs S plot for the preliminary design of the UCAV
Svals, Tvals_ceiling, Tvals_cruise, Tvals_maneuver, Tvals_climb = getTvals(Tguess= 500)
# Tvals, Svals_takeoff = getSvals()

# write a function to get the contours of TOGW for the UCAV over the T vs S plot
# range of Svals = 30 to 100
# range of Tvals = 20 to 300
W_guess = 500
def get_TOGWcontours(Svals, Tvals, W_guess = 700):
    TOGWvals = np.zeros((len(Svals), len(Tvals)))
    for i in range(len(Svals)):
        for j in range(len(Tvals)):
            TOGWvals[i, j] = getTOGW(Tvals[j], Svals[i], W_guess)
    return TOGWvals

fig, ax = plt.subplots(dpi = 300)
ax.plot(Svals, Tvals_ceiling, label = "Ceiling", color = "blue")
ax.plot(Svals, Tvals_cruise, label = "Cruise", color = "orange")
ax.plot(Svals, Tvals_maneuver, label = "Maneuver", color = "purple")
ax.plot(Svals, Tvals_climb, label = "Climb", color = "green")
# ax.plot(60, 400, "ro", label = "Design Point (45, 200)", markersize = 3)
# ax.plot(Svals_takeoff, Tvals, label = "Takeoff", color = "red")

# plot the contours of TOGW
Svalscontours = np.linspace(40, 100, 100)
Tvalscontours = np.linspace(20, 10000, 100)
print("\n\n\n\n FINDING TOGW CONTOURS") 
TOGW_contours = get_TOGWcontours(Svalscontours, Tvalscontours)
CS = ax.contour(Svalscontours, Tvalscontours, TOGW_contours, colors = "black", linestyles = "dashed")
ax.clabel(CS, inline = True, fontsize = 8)

ax.set_xlabel("Wing Area (S)(ft^2)")
ax.set_ylabel("Thrust (T)(lbf)")
ax.set_title("T vs S Plot for the UCAV (Payload 200lbs)")
ax.legend()
ax.grid()
plt.savefig("TSPlot_5.png")
plt.show()

# find the TOWG for the UCAV at design point
# T_designpoint = 400
# S_designpoint = 60
# TOGW_designpoint = getTOGW(T_designpoint, S_designpoint)
# print(f"\nThe estimated TOGW for T = {T_designpoint} lbf and S = {S_designpoint} ft^2 is {TOGW_designpoint:.2f} lb, ie {TOGW_designpoint/2.20462:.2f} kg")