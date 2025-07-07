# Differential-Equation-Integration-Comparison
Solves ODEs in MATLAB using Euler and RK4 integration, compares them with accurate solutions, and analyzes RMSEs for different step sizes.

---

# Overview
This repository contains MATLAB scripts for comparing **Euler** and **Runge-Kutta 4th Order (RK4)** methods for solving ordinary differential equations. The project checks both methods against an exact solution and high-accuracy reference data, calculates RMSEs, and visualizes the results.

---

# What’s Included
- **AccurateDataSP25.mat** — Reference solution data for Problem 1.
- **MATLAB Script** — Implements Euler and RK4 integration, compares numerical and exact solutions, computes RMSEs, and generates plots.
- **Functions** — Built-in Euler and RK4 integration functions.
- **answers structure** — Stores final computed values and RMSEs.

---

# Core Steps
**Euler Integration**
- Define ODE: dy/dx = -exp(-x²), y(0) = 1, from x = 0 to 2
- Integrate using Euler’s method with dx = 0.1
- Load accurate reference data
- Interpolate accurate data to match Euler grid
- Calculate RMSE vs reference
- Plot Euler vs Accurate

**Euler vs RK4 Integration**
- Define ODE: dy/dx = y * sin²(x), y(0) = π, from x = 0 to 3π
- Compute exact solution
- Integrate using Euler and RK4 for various step sizes (dx = π to π/16)
- Calculate RMSEs for each method and step size
- Plot Euler, RK4, and exact solution for each dx

---

# Save Results
- All final values and RMSEs are saved to an `answers` structure for easy reporting.

---

# Key Features
- Demonstrates step-by-step numerical integration
- Compares results with exact and reference solutions
- Calculates RMSEs for method accuracy
- Visualizes how integration accuracy improves with smaller steps
- Includes clear plots for method comparison

---

# How to Run
- Clone this repository.
- Place `AccurateDataSP25.mat` in the project folder.
- Open MATLAB and run the main script.
- Review generated plots and check the `answers` structure for RMSEs and final results.
