# Importing necessary libraries
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

# Data from the table
# MTOM is in kg, Swet is in m^2
data = {
    "MTOM_kg": [25, 50, 75, 130, 150, 170, 182, 205, 230, 630],
    "Swet_m2": [4.0, 4.8, 8.1, 9.8, 11.5, 8.4, 14.4, 13.3, 17.2, 33.9]
}

# Converting MTOM to lbs (1 kg = 2.20462 lbs) and Swet to ft^2 (1 m^2 = 10.7639 ft^2)
MTOM_lbs = np.array(data["MTOM_kg"]) * 2.20462
Swet_ft2 = np.array(data["Swet_m2"]) * 10.7639

# Taking logarithm of MTOM and Swet for linear regression in log-log space
log_MTOM = np.log(MTOM_lbs)
log_Swet = np.log(Swet_ft2)

# Fitting the linear regression model
log_model = LinearRegression()
log_model.fit(log_MTOM.reshape(-1, 1), log_Swet)

# Extracting the coefficients
B = log_model.coef_[0]  # Slope
A_log = log_model.intercept_  # Intercept in log space
A = np.exp(A_log)  # Convert back from log(A) to A

print(f"\n\nThe regression equation is: Swet = {A:.4f} * MTOM^{B:.4f}")

# Plotting the relationship on a log-log scale
plt.figure(figsize=(8, 6))
plt.scatter(MTOM_lbs, Swet_ft2, color='blue', label='Actual Swet')
plt.plot(MTOM_lbs, A * MTOM_lbs**B, color='red', label=f'Fitted Swet = {A:.2f} * MTOM^{B:.2f}')
plt.xscale('log')
plt.yscale('log')
plt.title('Log-Log Plot: Swet vs MTOM')
plt.xlabel('MTOM (lbs)')
plt.ylabel('Swet (ft²)')
plt.legend()
plt.grid(True, which="both", ls="--")
plt.show()

# defining a function for Swet based on the Regression from Gotten et. al
def Swet_UAV(MTOM, A = 2.5535, B = 0.6677):
    Swet = A * MTOM**B
    return Swet