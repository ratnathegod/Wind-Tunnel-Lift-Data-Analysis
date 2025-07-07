%{
Author: Ratnakaru Yalagathala
Assignment: ASEN2012 - Coding Challenge 1
Purpose: I made this script to analyze wind tunnel data and calculate a
weighted average. This script also uses error propagation, and runs a Monte Carlo simulation.
Creation Date: Feb 18 2025
Inputs: WT_Lift_Data.xlsx
Outputs: Weighted averages, uncertainties and plots.
%}

clear; clc; close all;

%% 1) Read and Clean Data
data = xlsread('WT_Lift_Data.xlsx');  
time_all = data(:,1);
lift_all = data(:,2:4);

% Now I plot the raw lift data so we can see the startup region
figure('Name','Raw Lift Data');
plot(time_all, lift_all(:,1), 'r'); hold on;
plot(time_all, lift_all(:,2), 'g');
plot(time_all, lift_all(:,3), 'b');
xlabel('Time (s)');
ylabel('Lift (N)');
title('Raw Lift Data (Before Removing NaNs or Startup)');
grid on;
legend('Group 1','Group 2','Group 3','Location','Best');

% Remove rows with missing data
valid_rows = ~any(isnan(lift_all), 2);
time_noNaN = time_all(valid_rows);
lift_noNaN = lift_all(valid_rows,:);

% Remove startup data and I used a threshold to do so
lift_threshold = 20;  
steady_mask = all(lift_noNaN > lift_threshold, 2);
time_clean = time_noNaN(steady_mask);
lift_clean = lift_noNaN(steady_mask, :);

% Plot cleaned data
figure('Name','Steady Lift Data After Threshold');
plot(time_clean, lift_clean(:,1), 'r'); hold on;
plot(time_clean, lift_clean(:,2), 'g');
plot(time_clean, lift_clean(:,3), 'b');
xlabel('Time (s)');
ylabel('Lift (N)');
title(['Steady Lift Data (Threshold = ', num2str(lift_threshold), ' N)']);
grid on;
legend('Group 1','Group 2','Group 3','Location','Best');

% Here I saved the number of rows
answers.LiftN = length(time_clean);

%% 2) Compute Weighted Average and Uncertainty
group_mean = mean(lift_clean);    % Mean
group_std  = std(lift_clean);     % Standard devation

weights = 1 ./ (group_std.^2);
weighted_avg = sum(weights .* group_mean) / sum(weights);
weighted_uncertainty = 1 / sqrt(sum(weights));

answers.LiftWeightedAve = weighted_avg;
answers.LiftWeightedStd = weighted_uncertainty;

fprintf('\nWeighted Average Lift: %.5g N\n', weighted_avg);
fprintf('Weighted Average Uncertainty: %.5g N\n', weighted_uncertainty);

% Now I plot the final data with weighted average lines
figure('Name','Lift Profiles with Weighted Average');
plot(time_clean, lift_clean(:,1), 'r'); hold on;
plot(time_clean, lift_clean(:,2), 'g');
plot(time_clean, lift_clean(:,3), 'b');
plot(time_clean, weighted_avg*ones(size(time_clean)), 'k-', 'LineWidth', 2);
plot(time_clean, (weighted_avg+weighted_uncertainty)*ones(size(time_clean)), 'k--');
plot(time_clean, (weighted_avg-weighted_uncertainty)*ones(size(time_clean)), 'k--');
xlabel('Time (s)');
ylabel('Lift (N)');
title('Lift Profiles (Steady Region) with Weighted Average');
grid on;
legend('Group 1','Group 2','Group 3','Weighted Avg','+/- 1\sigma','Location','Best');

%% 3) General Method for Lift Calculation 
rho = 1.24;    % this is Air Density                 
rho_uncert = 0.5/100 * rho;      
U_inf = 13;                     
U_inf_uncert = 8/100 * U_inf;   
alpha_deg = 4.2;                
alpha = alpha_deg * pi/180;     
alpha_uncert = 7/100 * alpha;   
S = 1.5;                        
S_uncert = 0.002;              
% from here I calculate the coefficient of lift
CL = 2 * pi * alpha;
L = 0.5 * CL * rho * U_inf^2 * S;

dL_rho   = (L / rho)   * rho_uncert;
dL_Uinf  = (2 * L / U_inf) * U_inf_uncert;
dL_alpha = (L / alpha)* alpha_uncert;
dL_S     = (L / S)    * S_uncert;
dL = sqrt(dL_rho^2 + dL_Uinf^2 + dL_alpha^2 + dL_S^2);

answers.LiftGMStd = dL;

fprintf('\nGeneral Method Lift Estimate: %.5g N\n', L);
fprintf('General Method Uncertainty: %.5g N\n', dL);

% even though the chart of error contributions is not required I made it so
% its more clear
contrib = [dL_rho, dL_Uinf, dL_alpha, dL_S];
figure('Name','General Method Error Contributions');
bar(contrib,'FaceColor',[0.2 0.6 0.8]); 
xticks(1:4);
xticklabels({'\delta\rho','\deltaU','\delta\alpha','\deltaS'});
ylabel('Uncertainty Contribution (N)');
title('Contributions to \DeltaL by General Method');
grid on;

%% 4)Monte Carlo Simulation
n_samples = 10000;
L_samples = zeros(n_samples,1); %first define the number of samples

for i = 1:n_samples  %loop and this will iterate 10k times
    rho_sample   = rho   + rho_uncert   * randn;
    U_inf_sample = U_inf + U_inf_uncert * randn;
    alpha_sample = alpha + alpha_uncert * randn;
    S_sample     = S     + S_uncert     * randn;
    
    CL_sample = 2 * pi * alpha_sample; %lift compute for each sample
    L_samples(i) = 0.5 * CL_sample * rho_sample * (U_inf_sample^2) * S_sample;
end

mc_mean = mean(L_samples); %monte carlo results statistics
mc_uncertainty = std(L_samples);
answers.LiftMCStd = mc_uncertainty;

fprintf('\nMonte Carlo Lift Estimate: %.5g N\n', mc_mean);
fprintf('Monte Carlo Uncertainty: %.5g N\n', mc_uncertainty);

figure('Name','Monte Carlo Histogram');
histogram(L_samples,50); hold on;
xline(mc_mean, 'k-', 'LineWidth',2);
xline(mc_mean + mc_uncertainty, 'k--');
xline(mc_mean - mc_uncertainty, 'k--');
xlabel('Lift (N)');
ylabel('Frequency');
title('Monte Carlo Simulation for Lift');
grid on;

%% 5)Compare Monte Carlo and General Method
figure('Name','Monte Carlo vs General Method PDF');
histogram(L_samples, 50, 'Normalization','pdf'); hold on;
x_vals = linspace(min(L_samples), max(L_samples), 200);
pdf_gen = (1/(dL*sqrt(2*pi))) * exp(-0.5*((x_vals - L)/dL).^2);
plot(x_vals, pdf_gen, 'r-', 'LineWidth',2);

xline(L,         'r--','LineWidth',2);
xline(L + dL,    'r--','LineWidth',1);
xline(L - dL,    'r--','LineWidth',1);
xlabel('Lift (N)');
ylabel('Probability Density');
title('Monte Carlo vs. General Method');
legend('Monte Carlo','General Method PDF','General Method Mean','General Method ±1\sigma','Location','best');
grid on;

% Errorbar Comparison from class
group_x = [1, 2, 3]; 
weighted_x = 4;    
gm_x = 5;          
mc_x = 6;          
figure('Name','Errorbar Comparison');
hold on; grid on;
errorbar(group_x, group_mean, group_std, 'o', 'LineWidth',1.5, 'DisplayName','Group Means');
errorbar(weighted_x, weighted_avg, weighted_uncertainty, 'o', 'LineWidth',1.5, 'DisplayName','Weighted Avg');
errorbar(gm_x, L, dL, 'o', 'LineWidth',1.5, 'DisplayName','Gen. Method');
errorbar(mc_x, mc_mean, mc_uncertainty, 'o', 'LineWidth',1.5, 'DisplayName','Monte Carlo');
set(gca,'XTick',[1 2 3 4 5 6], ...
    'XTickLabel',{'Group1','Group2','Group3','Weighted','GM','MC'});
ylabel('Lift (N)');
title('Comparison of Lift Averages & Uncertainties');
legend('Location','best');

%% Lastly Save results
save('answers.mat', 'answers');

%% Reflection Questions

%{
Do the general method and Monte Carlo simulation mean and uncertainty agree? 
What would happen to the value of the uncertainty of lift, δL or σL, if 
you performed these analyses with errors associated with 2σ of each uncertain quantity? 

Yes, the results from both methods are similar as they use the same input 
values and uncertainties to calculate lift. The Monte Carlo method samples 
these uncertainties randomly, while the general method estimates them analytically. 
Small differences may exist due to the randomness in the simulations.

If 2σ was used instead of 1σ, the uncertainty range would double, leading 
to a wider confidence interval for the calculated lift. This would encompass 
a larger range of possible values.

Is the weighted average and associated uncertainty of lift the same as 
taking the average of the three lift profiles and then taking the standard deviation? 

No, they would be different because a simple average treats all three lift 
measurements equally, whereas a weighted average gives more weight to measurements 
with lower uncertainty. This means more precise data has a greater influence 
on the final result.

Does the data in WT Lift Data.xlsx align with the Monte Carlo simulation and general method? 
To help answer this question, use the errorbar function to create a plot of 
the mean of the signals, weighted mean, mean general method calculation, and 
mean Monte Carlo simulation values. 

To do this, use the errorbar function to display ±1σ error bars:
errorbar(x, y, err) - Documentation: https://www.mathworks.com/help/matlab/ref/errorbar.html

Is there discrepancy in the data? 

Yes, the data aligns well with both methods, and an error bar plot can visually 
confirm this. By comparing the mean lift values, weighted mean, general method 
calculation, and Monte Carlo estimate (each with their ±1σ error bars), we can 
see the consistency of the data. Even if slight discrepancies exist, as long as 
the values overlap within their uncertainties and the trends match overall, 
the data remains consistent.
%}  
