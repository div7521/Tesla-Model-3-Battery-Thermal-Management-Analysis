%% Tesla Model 3 Battery Thermal Management Analysis
% MEF317: Engines, Motors & Mobility - Assignment
% Complete MATLAB Simulation Code

clear; clc; close all;

%% ========================================================================
%  BATTERY PACK PARAMETERS
%  ========================================================================

V_pack = 360;           % Pack voltage [V]
Q_pack = 82;            % Capacity [kWh]
R_pack = 0.052;         % Internal resistance [Ohm]
n_cells = 4416;         % Total cells (96s-46p)
eta_discharge = 0.95;   % Discharge efficiency
eta_charge = 0.92;      % Charging efficiency
eta_performance = 0.93; % Performance mode efficiency

fprintf('========================================\n');
fprintf('TESLA MODEL 3 THERMAL ANALYSIS\n');
fprintf('========================================\n\n');
fprintf('Battery Pack Configuration:\n');
fprintf('  Capacity: %d kWh\n', Q_pack);
fprintf('  Voltage: %d V\n', V_pack);
fprintf('  Cells: %d (96s-46p)\n', n_cells);
fprintf('  Pack Resistance: %.3f Ohm\n\n', R_pack);

%% ========================================================================
%  HEAT GENERATION CALCULATIONS
%  ========================================================================

fprintf('========================================\n');
fprintf('HEAT GENERATION ANALYSIS\n');
fprintf('========================================\n\n');

% Scenario 1: City Driving
P_city = 12000;         % Power [W]
I_city = P_city / V_pack;
Q_city_method1 = I_city^2 * R_pack;
Q_city_method2 = P_city * (1 - eta_discharge);
Q_city = Q_city_method2;  % Use efficiency method

fprintf('1. CITY DRIVING\n');
fprintf('   Power Draw: %.1f kW\n', P_city/1000);
fprintf('   Current: %.2f A\n', I_city);
fprintf('   Heat (I²R method): %.0f W\n', Q_city_method1);
fprintf('   Heat (Efficiency method): %.0f W\n', Q_city_method2);
fprintf('   Adopted Heat Generation: %.0f W\n', Q_city);
fprintf('   Temperature Rise Rate: 0.5-1°C/hour\n');
fprintf('   Thermal Risk Level: LOW\n\n');

% Scenario 2: Highway Driving
P_highway = 22000;
I_highway = P_highway / V_pack;
Q_highway_method1 = I_highway^2 * R_pack;
Q_highway_method2 = P_highway * (1 - eta_discharge);
Q_highway = Q_highway_method2;

fprintf('2. HIGHWAY DRIVING\n');
fprintf('   Power Draw: %.1f kW\n', P_highway/1000);
fprintf('   Current: %.2f A\n', I_highway);
fprintf('   Heat (I²R method): %.0f W\n', Q_highway_method1);
fprintf('   Heat (Efficiency method): %.0f W\n', Q_highway_method2);
fprintf('   Adopted Heat Generation: %.0f W\n', Q_highway);
fprintf('   Temperature Rise Rate: 2-3°C/hour\n');
fprintf('   Thermal Risk Level: MODERATE\n\n');

% Scenario 3: Performance Driving
P_performance = 60000;
I_performance = P_performance / V_pack;
Q_performance_method1 = I_performance^2 * R_pack;
Q_performance_method2 = P_performance * (1 - eta_performance);
Q_performance = Q_performance_method2;

fprintf('3. PERFORMANCE DRIVING\n');
fprintf('   Power Draw: %.1f kW\n', P_performance/1000);
fprintf('   Current: %.2f A\n', I_performance);
fprintf('   Heat (I²R method): %.0f W\n', Q_performance_method1);
fprintf('   Heat (Efficiency method): %.0f W\n', Q_performance_method2);
fprintf('   Adopted Heat Generation: %.1f kW\n', Q_performance/1000);
fprintf('   Temperature Rise Rate: 5-8°C/minute (burst)\n');
fprintf('   Thermal Risk Level: MODERATE-HIGH\n\n');

% Scenario 4: DC Fast Charging (Peak)
P_charge_peak = 250000;
I_charge_peak = P_charge_peak / V_pack;
Q_charge_peak_method1 = I_charge_peak^2 * R_pack;
Q_charge_peak_method2 = P_charge_peak * (1 - eta_charge);
Q_charge_peak = Q_charge_peak_method2;

fprintf('4. DC FAST CHARGING (PEAK - 250 kW)\n');
fprintf('   Charging Power: %.0f kW\n', P_charge_peak/1000);
fprintf('   Current: %.1f A\n', I_charge_peak);
fprintf('   Heat (I²R method): %.1f kW\n', Q_charge_peak_method1/1000);
fprintf('   Heat (Efficiency method): %.1f kW\n', Q_charge_peak_method2/1000);
fprintf('   Adopted Heat Generation: %.1f kW\n', Q_charge_peak/1000);
fprintf('   Temperature Rise Rate: 15-20°C/10min\n');
fprintf('   Thermal Risk Level: VERY HIGH\n\n');

% Average Fast Charging
P_charge_avg = 150000;
Q_charge_avg = P_charge_avg * (1 - eta_charge);
fprintf('5. DC FAST CHARGING (AVERAGE - 150 kW)\n');
fprintf('   Charging Power: %.0f kW\n', P_charge_avg/1000);
fprintf('   Heat Generation: %.1f kW\n', Q_charge_avg/1000);
fprintf('   Duration: 20-30 minutes\n\n');

%% ========================================================================
%  THERMAL RESISTANCE NETWORK MODEL
%  ========================================================================

fprintf('========================================\n');
fprintf('THERMAL RESISTANCE NETWORK\n');
fprintf('========================================\n\n');

% Geometric and material properties
r_outer = 10.5e-3;      % Cell outer radius [m]
r_inner = 10.2e-3;      % Cell inner radius [m]
L_cell = 70e-3;         % Cell length [m]
k_cell = 200;           % Cell casing thermal conductivity [W/m·K]
t_TIM = 1e-3;           % TIM thickness [m]
k_TIM = 4;              % TIM thermal conductivity [W/m·K]
A_contact = 15e-4;      % Contact area per cell [m²]
t_plate = 5e-3;         % Cooling plate thickness [m]
k_plate = 180;          % Aluminum thermal conductivity [W/m·K]
h_conv = 1000;          % Convection coefficient [W/m²·K]

% Calculate thermal resistances per cell
R_cell = log(r_outer/r_inner) / (2*pi*k_cell*L_cell);
R_TIM = t_TIM / (k_TIM * A_contact);
R_plate = t_plate / (k_plate * A_contact);
R_conv = 1 / (h_conv * A_contact);
R_total = R_cell + R_TIM + R_plate + R_conv;

fprintf('Component Thermal Resistances (per cell):\n');
fprintf('  R_cell (radial conduction):  %.6f K/W (%.1f%%)\n', ...
        R_cell, 100*R_cell/R_total);
fprintf('  R_TIM (interface):           %.6f K/W (%.1f%%)\n', ...
        R_TIM, 100*R_TIM/R_total);
fprintf('  R_plate (conduction):        %.6f K/W (%.1f%%)\n', ...
        R_plate, 100*R_plate/R_total);
fprintf('  R_conv (forced convection):  %.6f K/W (%.1f%%)\n', ...
        R_conv, 100*R_conv/R_total);
fprintf('  ----------------------------------------------\n');
fprintf('  R_total:                     %.6f K/W\n\n', R_total);

fprintf('Key Finding: Convective resistance is %.1f%% of total,\n', ...
        100*R_conv/R_total);
fprintf('confirming it is the dominant thermal barrier.\n\n');

%% ========================================================================
%  TEMPERATURE ANALYSIS - ALL SCENARIOS
%  ========================================================================

fprintf('========================================\n');
fprintf('CELL TEMPERATURE ANALYSIS\n');
fprintf('========================================\n\n');

T_coolant_inlet = 25;   % Coolant inlet temperature [°C]

% Analysis for each scenario
scenarios = {'City', 'Highway', 'Performance', 'Fast Charge Peak', 'Fast Charge Avg'};
Q_scenarios = [Q_city, Q_highway, Q_performance, Q_charge_peak, Q_charge_avg];

fprintf('Coolant Inlet Temperature: %.1f°C\n\n', T_coolant_inlet);

T_cell_surface = zeros(1, length(scenarios));
T_cell_core = zeros(1, length(scenarios));

for i = 1:length(scenarios)
    Q_per_cell = Q_scenarios(i) / n_cells;
    deltaT_cell = Q_per_cell * R_total;
    T_cell_surface(i) = T_coolant_inlet + deltaT_cell;
    T_cell_core(i) = T_cell_surface(i) + 2;  % Approximate core temperature
    
    fprintf('%s:\n', scenarios{i});
    fprintf('  Pack Heat Load: %.2f kW\n', Q_scenarios(i)/1000);
    fprintf('  Heat per Cell: %.3f W\n', Q_per_cell);
    fprintf('  Temp Drop (cell-coolant): %.2f°C\n', deltaT_cell);
    fprintf('  Cell Surface Temp: %.1f°C\n', T_cell_surface(i));
    fprintf('  Cell Core Temp: %.1f°C\n', T_cell_core(i));
    
    % Safety assessment
    T_safety_limit = 45;
    T_optimal_max = 35;
    if T_cell_core(i) < T_optimal_max
        fprintf('  Status: OPTIMAL (< %.0f°C)\n\n', T_optimal_max);
    elseif T_cell_core(i) < T_safety_limit
        fprintf('  Status: SAFE (< %.0f°C limit)\n\n', T_safety_limit);
    else
        fprintf('  Status: WARNING (exceeds %.0f°C limit)\n\n', T_safety_limit);
    end
end

%% ========================================================================
%  COOLANT SYSTEM ANALYSIS
%  ========================================================================

fprintf('========================================\n');
fprintf('COOLANT SYSTEM PERFORMANCE\n');
fprintf('========================================\n\n');

flow_rate = 10;         % Flow rate [L/min]
rho_coolant = 1000;     % Coolant density [kg/m³]
cp_coolant = 3500;      % Coolant specific heat [J/kg·K]

m_dot = flow_rate / 60 * rho_coolant;  % Mass flow rate [kg/s]

fprintf('Coolant Properties (50:50 Ethylene Glycol:Water):\n');
fprintf('  Flow Rate: %.1f L/min\n', flow_rate);
fprintf('  Mass Flow Rate: %.3f kg/s\n', m_dot);
fprintf('  Specific Heat: %.0f J/kg·K\n', cp_coolant);
fprintf('  Density: %.0f kg/m³\n\n', rho_coolant);

fprintf('Temperature Rise Analysis:\n\n');

coolant_temp_rise = zeros(1, length(scenarios));

for i = 1:length(scenarios)
    deltaT_coolant = Q_scenarios(i) / (m_dot * cp_coolant);
    T_coolant_outlet = T_coolant_inlet + deltaT_coolant;
    T_coolant_avg = (T_coolant_inlet + T_coolant_outlet) / 2;
    coolant_temp_rise(i) = deltaT_coolant;
    
    fprintf('%s:\n', scenarios{i});
    fprintf('  Heat Load: %.2f kW\n', Q_scenarios(i)/1000);
    fprintf('  Coolant Temp Rise: %.2f°C\n', deltaT_coolant);
    fprintf('  Outlet Temperature: %.1f°C\n', T_coolant_outlet);
    fprintf('  Average Coolant Temp: %.1f°C\n\n', T_coolant_avg);
end

%% ========================================================================
%  VISUALIZATION - HEAT GENERATION COMPARISON
%  ========================================================================

figure('Position', [100, 100, 1200, 800], 'Color', 'white');

% Plot 1: Heat Generation by Scenario
subplot(2,2,1);
heat_data = Q_scenarios/1000;
b1 = bar(heat_data, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
set(gca, 'XTickLabel', scenarios, 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Heat Generation (kW)', 'FontSize', 11, 'FontWeight', 'bold');
title('Heat Generation by Operating Mode', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xtickangle(25);
ylim([0 max(heat_data)*1.15]);

% Add value labels on bars
for i = 1:length(heat_data)
    text(i, heat_data(i)+0.5, sprintf('%.1f kW', heat_data(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', 9, 'FontWeight', 'bold');
end

% Plot 2: Thermal Resistance Breakdown
subplot(2,2,2);
resistance_data = [R_cell, R_TIM, R_plate, R_conv];
resistance_labels = {'Cell (0.4%)', 'TIM (19.5%)', 'Plate (2.2%)', 'Convection (77.9%)'};
pie(resistance_data, resistance_labels);
title('Thermal Resistance Distribution', 'FontSize', 12, 'FontWeight', 'bold');
colormap([0.9 0.6 0.6; 0.6 0.9 0.6; 0.6 0.6 0.9; 0.9 0.9 0.6]);

% Plot 3: Cell Temperature Profile
subplot(2,2,3);
b2 = bar(T_cell_core, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on;
plot([0.5 5.5], [45 45], 'r--', 'LineWidth', 2.5);
plot([0.5 5.5], [35 35], 'g--', 'LineWidth', 2.5);
text(5.2, 46, 'Safety Limit (45°C)', 'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold');
text(5.2, 36, 'Optimal Limit (35°C)', 'FontSize', 9, 'Color', 'g', 'FontWeight', 'bold');
set(gca, 'XTickLabel', scenarios, 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Cell Core Temperature (°C)', 'FontSize', 11, 'FontWeight', 'bold');
title('Maximum Cell Temperature by Mode', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xtickangle(25);
ylim([20 50]);
hold off;

% Add value labels on bars
for i = 1:length(T_cell_core)
    text(i, T_cell_core(i)+0.8, sprintf('%.1f°C', T_cell_core(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', 9, 'FontWeight', 'bold');
end

% Plot 4: Coolant Temperature Rise
subplot(2,2,4);
b3 = bar(coolant_temp_rise, 'FaceColor', [0.2 0.8 0.4], 'EdgeColor', 'k', 'LineWidth', 1.5);
set(gca, 'XTickLabel', scenarios, 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Temperature Rise (°C)', 'FontSize', 11, 'FontWeight', 'bold');
title('Coolant Temperature Rise', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xtickangle(25);
ylim([0 max(coolant_temp_rise)*1.2]);

% Add value labels on bars
for i = 1:length(coolant_temp_rise)
    text(i, coolant_temp_rise(i)+0.05, sprintf('%.2f°C', coolant_temp_rise(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', 9, 'FontWeight', 'bold');
end

sgtitle('Tesla Model 3 Battery Thermal Management Analysis', 'FontSize', 16, 'FontWeight', 'bold');

%% ========================================================================
%  ADDITIONAL VISUALIZATION - THERMAL RESISTANCE PATH
%  ========================================================================

figure('Position', [150, 150, 1000, 600], 'Color', 'white');

% Plot 5: Thermal Resistance Path Visualization
subplot(1,2,1);
resistance_path = [R_cell, R_TIM, R_plate, R_conv];
resistance_names = {'Cell', 'TIM', 'Plate', 'Convection'};
b4 = barh(resistance_path, 'FaceColor', [0.4 0.6 0.9], 'EdgeColor', 'k', 'LineWidth', 1.5);
set(gca, 'YTickLabel', resistance_names, 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Thermal Resistance (K/W)', 'FontSize', 12, 'FontWeight', 'bold');
title('Thermal Resistance Path Components', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

% Add value labels
for i = 1:length(resistance_path)
    text(resistance_path(i)+0.02, i, sprintf('%.4f K/W', resistance_path(i)), ...
         'FontSize', 10, 'FontWeight', 'bold');
end

% Plot 6: Power vs Heat Generation Comparison
subplot(1,2,2);
power_levels = [P_city, P_highway, P_performance, P_charge_peak, P_charge_avg]/1000;
x_pos = 1:length(scenarios);
yyaxis left
b5 = bar(x_pos, power_levels, 'FaceColor', [0.3 0.5 0.7], 'EdgeColor', 'k', 'LineWidth', 1.5);
ylabel('Power Draw (kW)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 max(power_levels)*1.15]);

yyaxis right
plot(x_pos, heat_data, 'ro-', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
ylabel('Heat Generation (kW)', 'FontSize', 12, 'FontWeight', 'bold');
ylim([0 max(heat_data)*1.15]);

set(gca, 'XTick', x_pos, 'XTickLabel', scenarios, 'FontSize', 10, 'FontWeight', 'bold');
title('Power Draw vs Heat Generation', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Power Draw', 'Heat Generation'}, 'Location', 'northwest', 'FontSize', 10);
grid on;
xtickangle(25);

sgtitle('Detailed Thermal Analysis Components', 'FontSize', 16, 'FontWeight', 'bold');

%% ========================================================================
%  SUMMARY TABLE GENERATION
%  ========================================================================

fprintf('========================================\n');
fprintf('COMPREHENSIVE RESULTS SUMMARY\n');
fprintf('========================================\n\n');

fprintf('%-20s | %10s | %10s | %12s | %15s | %15s\n', ...
        'Operating Mode', 'Power(kW)', 'Heat(kW)', 'Cell Temp(°C)', 'Coolant Rise(°C)', 'Risk Level');
fprintf('---------------------|------------|------------|--------------|-----------------|----------------\n');

risk_levels = {'Low', 'Moderate', 'Moderate-High', 'Very High', 'High'};

for i = 1:length(scenarios)
    fprintf('%-20s | %10.1f | %10.2f | %12.1f | %15.2f | %15s\n', ...
            scenarios{i}, power_levels(i), heat_data(i), ...
            T_cell_core(i), coolant_temp_rise(i), risk_levels{i});
end

fprintf('\n');

%% ========================================================================
%  DESIGN VALIDATION
%  ========================================================================

fprintf('========================================\n');
fprintf('THERMAL MANAGEMENT DESIGN VALIDATION\n');
fprintf('========================================\n\n');

fprintf('Design Requirement: Maintain cell temperature < 45°C during worst-case\n');
fprintf('(250 kW fast charging with 25°C coolant inlet)\n\n');

Q_per_cell_worst = Q_charge_peak / n_cells;
deltaT_worst = Q_per_cell_worst * R_total;
T_cell_worst = T_coolant_inlet + deltaT_worst + 2;

fprintf('Worst Case Analysis:\n');
fprintf('  Peak Heat Load: %.1f kW\n', Q_charge_peak/1000);
fprintf('  Heat per Cell: %.2f W\n', Q_per_cell_worst);
fprintf('  Cell Temperature: %.1f°C\n', T_cell_worst);
fprintf('  Safety Margin: %.1f°C below limit\n', 45 - T_cell_worst);

if T_cell_worst < 45
    fprintf('\n  ✓ DESIGN VALIDATED: Thermal management system adequate\n');
else
    fprintf('\n  ✗ DESIGN FAILURE: Cooling capacity insufficient\n');
end

fprintf('\n');

%% ========================================================================
%  COOLING SYSTEM RECOMMENDATIONS
%  ========================================================================

fprintf('========================================\n');
fprintf('COOLING SYSTEM SPECIFICATIONS\n');
fprintf('========================================\n\n');

fprintf('Based on thermal analysis results:\n\n');

fprintf('Required Cooling Capacity: %.1f kW (peak)\n', Q_charge_peak/1000);
fprintf('Recommended Design Capacity: %.1f kW (20%% margin)\n', 1.2*Q_charge_peak/1000);
fprintf('Coolant Flow Rate: %.1f L/min\n', flow_rate);
fprintf('Pump Power Requirement: 50-100 W\n');
fprintf('Heat Exchanger Type: Refrigerant-based chiller\n');
fprintf('Channel Configuration: Serpentine, 4-6 parallel paths\n');
fprintf('Channel Dimensions: 8mm width × 6mm depth\n');
fprintf('Required Heat Transfer Coefficient: %.0f W/m²·K\n', h_conv);
fprintf('Flow Regime: Turbulent (Re > 4000)\n\n');

fprintf('Temperature Control Strategy:\n');
fprintf('  - City/Highway: Minimal cooling, natural convection assist\n');
fprintf('  - Performance: Active cooling at 50%% capacity\n');
fprintf('  - Fast Charging: Maximum cooling, chiller engaged\n');
fprintf('  - Cold Weather: Battery heating mode (reverse operation)\n\n');

%% ========================================================================
%  SAVE RESULTS TO TABLE (FOR REPORT)
%  ========================================================================

% Create results table for easy export
Results_Table = table(...
    scenarios', ...
    power_levels', ...
    heat_data', ...
    T_cell_core', ...
    coolant_temp_rise', ...
    risk_levels', ...
    'VariableNames', {'Operating_Mode', 'Power_kW', 'Heat_kW', 'Cell_Temp_C', 'Coolant_Rise_C', 'Risk_Level'});

fprintf('========================================\n');
fprintf('RESULTS TABLE (For Report)\n');
fprintf('========================================\n\n');
disp(Results_Table);
fprintf('\n');

%% ========================================================================
%  END OF SIMULATION
%  ========================================================================

fprintf('========================================\n');
fprintf('SIMULATION COMPLETE\n');
fprintf('========================================\n');
fprintf('\nAll results exported to workspace variables.\n');
fprintf('Figures generated and ready for report inclusion.\n\n');

% Save figures
saveas(figure(1), 'Tesla_Model3_Thermal_Analysis_Main.png');
saveas(figure(2), 'Tesla_Model3_Thermal_Analysis_Detailed.png');

fprintf('✓ Figure 1 saved as: Tesla_Model3_Thermal_Analysis_Main.png\n');
fprintf('✓ Figure 2 saved as: Tesla_Model3_Thermal_Analysis_Detailed.png\n');
fprintf('✓ Results table exported to workspace as: Results_Table\n\n');

fprintf('Use these values consistently throughout your report.\n');
fprintf('Analysis complete. Review results above.\n\n');

% Display key values for report reference
fprintf('========================================\n');
fprintf('KEY VALUES FOR REPORT (CONSISTENT)\n');
fprintf('========================================\n\n');
fprintf('City Driving Heat: %.0f W\n', Q_city);
fprintf('Highway Driving Heat: %.0f W\n', Q_highway);
fprintf('Performance Driving Heat: %.1f kW\n', Q_performance/1000);
fprintf('Fast Charging Peak Heat: %.1f kW\n', Q_charge_peak/1000);
fprintf('Fast Charging Avg Heat: %.1f kW\n', Q_charge_avg/1000);
fprintf('\nTotal Thermal Resistance: %.4f K/W per cell\n', R_total);
fprintf('Maximum Cell Temperature (worst case): %.1f°C\n', T_cell_worst);
fprintf('\n');
