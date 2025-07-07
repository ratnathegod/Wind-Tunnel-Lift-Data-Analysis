# Wind-Tunnel-Lift-Data-Analysis
Processes wind tunnel data, calculates lift with error propagation, and compares general method and Monte Carlo results for uncertainty analysis.

# Overview
This repository contains MATLAB scripts for analyzing wind tunnel lift data. The project processes raw lift measurements, removes startup and invalid data, and computes lift using two approaches: a general analytical method with error propagation and a Monte Carlo simulation. It then compares the results, visualizes uncertainty, and saves final outputs for further analysis.

# What’s Included
- WT_Lift_Data.xlsx — Raw wind tunnel lift measurements for three groups.
- MATLAB Script — Cleans data, computes weighted averages, performs lift calculations, runs the Monte Carlo simulation, and generates comparison plots.
- Flowchart — Visualizes the full workflow: from data cleaning to final result saving.
- answers.mat — Output file storing final results for reuse.

# Core Steps
- Load & Clean Data
- Read WT_Lift_Data.xlsx
- Remove NaNs and startup region below a lift threshold
- Plot raw vs. cleaned data for verification

Compute Weighted Average
- Calculate group means and standard deviations
- Compute the weighted average lift and its uncertainty
- Plot lift profiles with weighted mean ± 1σ bands

General Method Calculation
- Use air density, free-stream velocity, angle of attack, and wing area
- Apply lift equations with error propagation
- Visualize uncertainty contributions

Monte Carlo Simulation
- Randomly sample input parameters within their uncertainties
- Compute lift for 10,000 samples
- Plot a histogram of lift distribution with mean ± 1σ

Compare Methods
- Plot PDFs of Monte Carlo and General Method
- Create error bar plots to compare group means, weighted average, general method, and Monte Carlo results side by side

# Save Results
- Save all outputs and statistics to answers.mat for reproducibility

# Key Features
- Cleans noisy wind tunnel data automatically
- Computes lift analytically and statistically
- Visualizes uncertainty with clear plots
- Compares multiple methods for validation
-  Saves outputs for easy reporting

# How to Run
- Clone this repository.
- Place WT_Lift_Data.xlsx in the project folder.
- Open MATLAB and run the main script.
- Check generated plots and answers.mat for final results.

