%% UCA Design: Uniform-Like Steering + Multiple Nulls
clear; clc; close all;

%% 1. Parameters
N = 8;                  % Number of antenna elements
R_lambda = (8.4*1e-2)/(12.5*1e-2);        % Radius / Wavelength ratio
circumference = 2*pi*R_lambda; 
M = 5;                  % Number of Phase Modes
m_indices = (-(M-1)/2 : (M-1)/2)'; 

% --- Steering and Null Settings ---
phi_0_deg = -100;                 % Main lobe direction (Degrees)
phi_0 = deg2rad(phi_0_deg);       
null_angles_deg = [-25];     % Angles to be nulled (Degrees)
null_angles = deg2rad(null_angles_deg);
phi = linspace(-pi, pi, 4000);    % High resolution for precision

%% 2. Point Nulling (Constrained Least Squares)
num_nulls = length(null_angles);
C = zeros(num_nulls + 1, M);    
b = zeros(num_nulls + 1, 1);    

% Main lobe constraint
C(1, :) = exp(1j * m_indices * 0); 
b(1) = 1;

% Null constraints
for k = 1:num_nulls
    C(k+1, :) = exp(1j * m_indices * (null_angles(k) - phi_0));
    b(k+1) = 0;
end

% Calculating Phase Mode Weights
w_PM = C' * ((C * C') \ b);

%% 3. Physical Antenna Weights
w_physical = zeros(N, 1);
phi_n = 2 * pi * (0:N-1)' / N; 
for n = 1:N
    sum_val = 0;
    for i = 1:M
        m = m_indices(i);
        bessel_term = besselj(m, circumference);
        % Synthesis formula for UCA Phase Modes
        term = (w_PM(i) / (N * (1j^m) * bessel_term)) * exp(1j * m * phi_n(n));
        sum_val = sum_val + term;
    end
    w_physical(n) = sum_val;
end

% Display the Physical Weight Matrix
fprintf('--- Physical Antenna Weights (w_physical) ---\n');
disp(table((1:N)', abs(w_physical), rad2deg(angle(w_physical)), ...
    'VariableNames', {'Element_No', 'Amplitude', 'Phase_Degrees'}));

%% 4. Analysis and Beamwidth Calculations
pattern = zeros(size(phi));
for i = 1:M
    m = m_indices(i);
    pattern = pattern + (w_PM(i) * exp(1j * m * (phi - phi_0)));
end

pattern_linear = abs(pattern) / max(abs(pattern));
pattern_db = 20*log10(pattern_linear);

% --- HPBW and BWNN Calculation ---
[~, main_lob_idx] = max(pattern_db);

% HPBW (-3 dB points)
target_hp = -3;
idx_left = find(pattern_db(1:main_lob_idx) <= target_hp, 1, 'last');
idx_right = find(pattern_db(main_lob_idx:end) <= target_hp, 1, 'first') + main_lob_idx - 1;
hpbw = abs(rad2deg(phi(idx_right) - phi(idx_left)));

% BW Null-to-Null (First local minima around main lobe)
[~, local_mins] = findpeaks(-pattern_db);
null_left_idx = local_mins(find(local_mins < main_lob_idx, 1, 'last'));
null_right_idx = local_mins(find(local_mins > main_lob_idx, 1, 'first'));

if ~isempty(null_left_idx) && ~isempty(null_right_idx)
    bwnn = abs(rad2deg(phi(null_right_idx) - phi(null_left_idx)));
else
    bwnn = NaN; % If distinct nulls are not formed
end

fprintf('--- Beamwidth Analysis ---\n');
fprintf('HPBW: %.2f degrees\n', hpbw);
fprintf('BW Null-to-Null: %.2f degrees\n', bwnn);

%% 5. Visualization (dB Format)
figure('Name', 'UCA Analysis Results', 'Color', 'w', 'Position', [100, 100, 1200, 500]);

% Cartesian Plot
subplot(1,2,1);
plot(rad2deg(phi), pattern_db, 'LineWidth', 2, 'Color', [0 0.447 0.741]); 
grid on; hold on;
xline(phi_0_deg, '--g', 'Main Lobe', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
for k = 1:num_nulls
    xline(null_angles_deg(k), '--r', ['Null ' num2str(k)], 'LineWidth', 1.5);
end

% HPBW Visualization
yline(-3, 'k:', 'HPBW (-3dB)');
plot([rad2deg(phi(idx_left)) rad2deg(phi(idx_right))], [-3 -3], 'mo-', 'MarkerFaceColor', 'm');
ylim([-60 5]); xlim([-180 180]);
title(['Radiation Pattern (HPBW: ' num2str(hpbw, '%.1f') '^\circ)']);
xlabel('Angle (Degrees)'); ylabel('Gain (dB)');

% Polar Plot
subplot(1,2,2);
db_limit = -60; 
polar_data = pattern_db;
polar_data(polar_data < db_limit) = db_limit; 
polarplot(phi, polar_data, 'LineWidth', 2, 'Color', [0 0.447 0.741]);
hold on;
polarplot([phi_0 phi_0], [db_limit 0], '--g', 'LineWidth', 1.5); 
for k = 1:num_nulls
    polarplot([null_angles(k) null_angles(k)], [db_limit 0], '--r', 'LineWidth', 1.5); 
end
ax = gca; ax.RLim = [db_limit 0];
title('Polar Radiation Pattern (dB)');