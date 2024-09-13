# writing a python script for calulcating the initial TOGW of the UCAV (using Ch2 in Aircraft Design Metabook)

import numpy as np

# defining a function for the empty weight fraction regression
def emptyweightfraction(W0, A = 0.97, C = -0.06):
    W0 = W0 / 2.20462 # converting lb to kg
    WeW0 = (A) * (W0)**(C)
    return WeW0

# defining the constants
c = 0.3 # tsfc in lb/(lbf.hr)
LoverD = 10 # L/D ratio
V = 220 # cruise speed in mph
R = 500 # range in miles
E = 10 # endurance in hours

# defining the known weight ratios
W1W0 = 0.97 #takeoff
W2W1 = 0.985 #climb
# W3W2 is the weight ratio for the cruise segment found using the Breguet range formula
W4W3 = 0.99 #descent
# W5W4 is the weight ratio for the loiter segment found using the Breguet range (Endurance) formula
W6W5 = 0.995 #landing

# calculating the weight ratios for the cruise and loiter segments
W3W2 = np.exp(-R*c/(V*LoverD))
print(f"\nThe weight ratio for the cruise segment is {W3W2:.4f}")
W5W4 = np.exp(-E*c/(LoverD))
print(f"\nThe weight ratio for the loiter segment is {W5W4:.4f}")

# calculating the fuel weight fraction with 6% fuel reserve
WfW0 = 1 - (W6W5 * W5W4 * W4W3 * W3W2 * W2W1 * W1W0) * 1.06
print(f"\nThe fuel weight fraction is {WfW0:.4f}")

Wpayload = 100 # payload weight in lb

# implementing a fixed point iteration to find the initial TOGW
# initial guess for the TOGW
W0 = 400
# tolerance
tol = 1e-6
delta = 2 * tol
i = 0
while delta > tol and i < 100:
    WeW0 = emptyweightfraction(W0, 1.53, -0.16)
    W0_new = Wpayload / (1 - WeW0 - WfW0)
    delta = abs(W0_new - W0) / abs(W0_new)
    W0 = W0_new
    i += 1
    print(f"\nIteration {i}: W0 = {W0:.2f} lb")

print(f"\n\nThe TOGW is {W0:.2f} lb, ie {W0/2.20462:.2f} kg")




