# Wind-Tunnel-Lift-Data-Analysis
Processes wind tunnel data, calculates lift with error propagation, and compares general method and Monte Carlo results for uncertainty analysis.

---

# Overview
This repository contains MATLAB scripts for analyzing wind tunnel lift data. The project processes raw lift measurements, removes startup and invalid data, computes lift using two approaches — a general analytical method with error propagation and a Monte Carlo simulation — and compares the results with clear visualizations.

---

# What’s Included
- **WT_Lift_Data.xlsx** — Raw wind tunnel lift measurements for three groups.
- **MATLAB Script** — Cleans data, computes weighted averages, calculates lift using both methods, runs the Monte Carlo simulation, and generates comparison plots.
- **Flowchart** — Visualizes the full workflow: from data cleaning to final result saving.
- **answers.mat** — Output file storing all final lift values and uncertainty estimates.

---

# Core Steps
**Load & Clean Data**
- Read `WT_Lift_Data.xlsx`
- Remove rows with NaNs and filter out startup regions below a lift threshold
- Plot raw vs. cleaned data to verify

**Compute Weighted Average**
- Calculate means and standard deviations for each group
- Compute the weighted average lift and its uncertainty
- Plot lift profiles with weighted mean and ±1σ bands

**General Method Calculation**
- Use air density, free-stream velocity, angle of attack, and wing area
- Apply lift equations with uncertainty propagation
- Plot contributions of each source of error

**Monte Carlo Simulation**
- Randomly sample input parameters within their uncertainty ranges
- Calculate lift for 10,000 samples
- Plot histogram of simulated lift results with mean ±1σ

**Compare Methods**
- Plot PDFs of Monte Carlo and General Method results side by side
- Create error bar plots to compare group means, weighted average, general method, and Monte Carlo lift estimates

---

# Save Results
- All final lift values, uncertainties, and figures are saved to `answers.mat` for easy reuse and reporting.

---

# Key Features
- Automatically cleans noisy wind tunnel data
- Calculates lift analytically and statistically
- Visualizes uncertainty clearly with multiple plots
- Compares multiple methods for validation
- Saves all results for reproducibility

---

# How to Run
- Clone this repository.
- Place `WT_Lift_Data.xlsx` in the project folder.
- Open MATLAB and run the main script.
- Review the generated plots and check `answers.mat` for all final lift results and uncertainty values.

