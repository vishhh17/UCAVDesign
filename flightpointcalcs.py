def isa_atmosphere(altitude_m):
    # Constants
    R = 8.31447  # Universal gas constant [J/(mol*K)]
    g0 = 9.80665  # Gravitational acceleration [m/s^2]
    T0 = 288.15  # Sea level standard temperature [K]
    P0 = 101325  # Sea level standard pressure [Pa]
    L = 0.0065  # Temperature lapse rate [K/m]
    M = 0.0289644  # Molar mass of Earth's air [kg/mol]
    
    # Calculate temperature at altitude
    T = T0 - L * altitude_m
    
    # Calculate pressure at altitude
    P = P0 * (T / T0) ** (g0*M / (L * R))
    
    # Return temperature (K) and pressure (Pa)
    return T, P

def freestream_velocity_and_density(altitude_m, mach_number):
    R = 287.05  # Specific gas constant for air [J/(kg*K)]
    g0 = 9.80665  # Gravitational acceleration [m/s^2]
    T0 = 288.15  # Sea level standard temperature [K]
    P0 = 101325  # Sea level standard pressure [Pa]
    L = 0.0065  # Temperature lapse rate [K/m]

    # Speed of sound at sea level [m/s]
    a0 = 340.294
    
    # Calculate temperature and pressure at altitude
    T, P = isa_atmosphere(altitude_m)
    
    # Speed of sound at altitude [m/s]
    a = (a0 ** 2 * T / T0) ** 0.5
    
    # Freestream velocity [m/s]
    V = mach_number * a
    
    # Calculate air density at altitude [kg/m^3]
    rho = P / (R * T)
    
    return V, rho

def calculate_dynamic_viscosity(altitude_m):
    # Sutherland's constants for air
    C1 = 1.458e-6  # kg/(m·s·K^0.5)
    S = 110.4  # Sutherland's temperature in K
    T, _ = isa_atmosphere(altitude_m)
    
    # Sutherland's formula for dynamic viscosity
    mu = C1 * T**(3/2) / (T + S)
    
    return mu

def calculate_reynolds_number(altitude_m, mach_number, characteristic_length = 1.0):
    
    T, _ = isa_atmosphere(altitude_m)
    V, rho = freestream_velocity_and_density(altitude_m, mach_number)


    # Calculate dynamic viscosity
    mu = calculate_dynamic_viscosity(T)
    
    # Calculate Reynolds number
    Re = (rho * V * characteristic_length) / mu
    
    return Re

# write a function to calculate the coefficient of lift
def calculate_cl(weight, wing_area, flight_speed, air_density):
    force_due_to_gravity = weight * 9.81  # in Newtons
    coefficient_of_lift = (2 * force_due_to_gravity) / (air_density * flight_speed**2 * wing_area)
    return coefficient_of_lift

import math

def breguet_range(airspeed, specific_fuel_consumption, initial_weight, final_weight):
    """
    Calculate the range of an aircraft using the Breguet range formula for jet aircraft.

    Parameters:
    airspeed (float): True airspeed of the aircraft in m/s.
    specific_fuel_consumption (float): Specific fuel consumption of the engine in (kg/s)/N.
    initial_weight (float): Initial weight of the aircraft including fuel in kg.
    final_weight (float): Final weight of the aircraft after fuel burn in kg.

    Returns:
    float: Theoretical range of the aircraft in kilometers.
    """
    # Convert specific fuel consumption from (kg/s)/N to (kg/s)/kg by using gravity
    g = 9.81  # gravity in m/s^2
    sfc_kg_per_kg_s = specific_fuel_consumption * g  # Convert SFC to kg/(kg*s)

    # Calculate the natural logarithm of the weight ratio
    weight_ratio = initial_weight / final_weight
    ln_weight_ratio = math.log(weight_ratio)

    # Calculate the range
    range_km = (airspeed / sfc_kg_per_kg_s) * ln_weight_ratio / 1000  # Convert m to km

    return range_km

# write a function to find the final weight of the aircraft after fuel burn using the Breguet range formula
import numpy as np

def breguet_final_weight(airspeed, sfc_kg_per_kg_s, initial_weight, range_km):
    range_m = range_km * 1000  # Convert range from km to m
    L_over_D = 10
    # print("range_m: ", range_m)
    # print("(sfc_kg_per_kg_s/airspeed)*range_m: ", (sfc_kg_per_kg_s/L_over_D)*range_m)
    weight_ratio = np.exp((sfc_kg_per_kg_s/L_over_D)*range_m)
    # print("weight_ratio: ", weight_ratio)
    final_weight = initial_weight / weight_ratio
    return final_weight

# def breguet_final_weight(airspeed, sfc_kg_per_kg_s, initial_weight, range_km):
#     range_m = range_km * 1000  # Convert range from km to m
#     print("range_m", range_m)
#     print("range_m / (airspeed / sfc_kg_per_kg_s)", range_m / (airspeed / sfc_kg_per_kg_s))
#     weight_ratio = np.exp(-range_m / (airspeed / sfc_kg_per_kg_s))
#     final_weight = initial_weight * weight_ratio
#     return final_weight


# Example: Calculate freestream velocity and density at 10,000 meters altitude and Mach 0.85
altitude_m = 4000  # Altitude in meters
mach_number = 0.3  # Mach number

freestream_velocity, air_density = freestream_velocity_and_density(altitude_m, mach_number)
print(f"\n\nAltitude: {altitude_m} m")
print(f"Mach Number: {mach_number}")
print(f"Temperature: {isa_atmosphere(altitude_m)[0]:.2f} K")
print(f"Pressure: {isa_atmosphere(altitude_m)[1]:.2f} Pa")
print(f"Freestream Velocity: {freestream_velocity:.2f} m/s or {freestream_velocity * 18/5:.2f} km/h \n or {freestream_velocity * 2.237:.2f} mph or {freestream_velocity * 3.281:.2f} ft/s")
# air density in kg/m^3 and slugs/ft^3 and lb/ft^3
print(f"Air Density: {air_density:.2f} kg/m^3 or {air_density * 0.00194032:.6f} slugs/ft^3 or {air_density * 0.062428:.2f} lb/ft^3")
print(f"Dynamic Viscosity: {calculate_dynamic_viscosity(altitude_m):.2e} kg/(m*s)")
print(f"Reynolds Number: {calculate_reynolds_number(altitude_m, mach_number):.2e}")


# implement a nonlinear solver to find the altitude at which Re = 6.5e6 at Mach 0.734
# from scipy.optimize import fsolve

# def reynolds_number_residual(altitude_m, mach_number, re_target):
#     return calculate_reynolds_number(altitude_m, mach_number) - re_target

# # Target Reynolds number
# re_target = 6.5e6

# # Initial guess for altitude
# altitude_guess = 10000

# # Solve for altitude
# altitude_solution = fsolve(reynolds_number_residual, altitude_guess, args=(0.734, re_target))[0]
# print(f"\n\nAltitude for Re = {re_target:.2e} at Mach 0.734: {altitude_solution:.2f} m")

# # calculate Cl at 9000 m altitude and Mach 0.6
# weight = 225 #*(7/8) #kgs
# wing_area = 4.375 #m^2
# altitude_m = 100  # Altitude in meters
# mach_number = 0.55  # Mach number
# print(f"\n\nCalculating Coefficient of Lift at {altitude_m} m altitude and Mach {mach_number}...")
# print(f"Weight: {weight} kg")
# print(f"Wing Area: {wing_area} m^2")
# freestream_velocity, air_density = freestream_velocity_and_density(altitude_m, mach_number)
# print(f"Freestream Velocity: {freestream_velocity:.2f} m/s")
# print(f"Air Density: {air_density:.2f} kg/m^3")
# cl = calculate_cl(weight, wing_area, freestream_velocity, air_density)
# print(f"\nCoefficient of Lift at {altitude_m} m altitude and Mach {mach_number}: {cl:.2f}")




# # find the final weight of the aircraft after fuel burn using the Breguet range formula
# airspeed = 151.90  # m/s
# specific_fuel_consumption = 1.11*10**(-5)  # kg/(N*s)
# initial_weight = 350  # kg
# range_km = 400  # km
# final_weight = breguet_final_weight(airspeed, specific_fuel_consumption, initial_weight, range_km)
# print(f"\n\nFinding Final Weight of the Aircraft after Fuel Burn using Breguet Range Formula:")
# print(f"Initial Weight of the Aircraft: {initial_weight} kg")
# print(f"Range: {range_km} km")
# print(f"Final Weight of the Aircraft after Fuel Burn: {final_weight:.2f} kg")