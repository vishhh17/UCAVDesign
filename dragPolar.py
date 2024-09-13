import numpy as np
import matplotlib.pyplot as plt

# defining the S_wet function for the UCAV
def S_wet(W0, c = 0.8565, d = 0.5423):
    S_wet = 10**c * W0**d
    return S_wet

# defining the drag polar function for the UCAV
def drag_polar(CL, CD0, e, AR = 5):
    CD = CD0 + (1 / (np.pi * e * AR)) * CL**2
    return CD

# assuming Cf = 0.0035 for the UCAV (Skin friction coefficient)
Cf = 0.0035

# assuming Sref = 40 ft^2 for the UCAV (Reference wing area) (fig 4.4 of the Metabook)
Sref = 40 # ft^2

# calculating the CD0 for the UCAV
CD0 = Cf * S_wet(793) / Sref
print(f"\nThe CD0 for the UCAV is {CD0:.4f}")

# defining the Oswald efficiency factors for the UCAV during clean and takeoff configurations 
# (as landing characteristics not necessary due to the final attack phase)
e_clean = 0.85
e_takeoff = 0.8

# defining the del_CD0 for the takeoff and clean configurations
delCD0_clean = 0
delCD0_takeoff = 0.01

# plotting the drag polars for the UCAV in clean and takeoff configurations
CL = np.linspace(-2, 2, 100)
CD_clean = drag_polar(CL, CD0 + delCD0_clean, e_clean)
CD_takeoff = drag_polar(CL, CD0 + delCD0_takeoff, e_takeoff)

# function to find the max CL/CD ratio for the UCAV, and return the max CL/CD with the corresponding CL and CD
def max_CLoverCD(CL, CD):
    CL_CD = CL / CD
    max_CL_CD = np.max(CL_CD)
    corres_CL = CL[np.argmax(CL_CD)]
    corres_CD = CD[np.argmax(CL_CD)]
    print(f"\nThe maximum CL/CD ratio is {max_CL_CD:.4f} at CL = {corres_CL:.4f} and CD = {corres_CD:.4f}")
    return max_CL_CD, corres_CL, corres_CD




# plotting the drag polars for the UCAV in clean and takeoff configurations with the max CL/CD marked for both polars and labelled
fig, ax = plt.subplots(dpi = 300)
ax.plot(CD_clean, CL, label = "Drag Polar : Clean Configuration")
ax.plot(CD_takeoff, CL, label = "Drag Polar : Takeoff Configuration")
# finding the max CL/CD for the clean configuration
cleanMaxCLoverCD, cleanCL, cleanCD = max_CLoverCD(CL, CD_clean)
# finding the max CL/CD for the takeoff configuration
takeoffMaxCLoverCD, takeoffCL, takeoffCD = max_CLoverCD(CL, CD_takeoff)
# marking the max CL/CD for the clean configuration
ax.plot(cleanCD, cleanCL, "ro", label = f"Max CL/CD : Clean Configuration\nCL/CD = {cleanMaxCLoverCD:.4f}")
# marking the max CL/CD for the takeoff configuration
ax.plot(takeoffCD, takeoffCL, "go", label = f"Max CL/CD : Takeoff Configuration\nCL/CD = {takeoffMaxCLoverCD:.4f}")
ax.set_xlabel("CD")
ax.set_ylabel("CL")
ax.set_title("Drag Polars for the UCAV")
ax.legend()
ax.grid()
plt.savefig("dragPolar.png")
plt.show()

