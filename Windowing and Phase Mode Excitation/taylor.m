%% UCA Design: Taylor Windowing + Synthesis + Null Steering
clear; clc; close all;

%% 1. Parameters
N = 8;                           % Number of antenna elements (Increased for better null resolution)
R_lambda = (84e-3)/(0.125);                   % Radius / Wavelength ratio
circumference = 2*pi*R_lambda; 
M = 6;                            % Number of phase modes (Recommended M <= N)
m_indices = (-(M-1)/2 : (M-1)/2)'; 

% --- Taylor & Steering Settings ---
SLL_dB = 25;                      % Desired Side Lobe Level (dB)
nbar = 4;                         % Taylor parameter
phi_0_deg = 45;                   % Main lobe direction (Degrees)
phi_nulls_deg = [-70, 100];       % Directions to be NULLIFIED
phi_0 = deg2rad(phi_0_deg);       
phi_nulls = deg2rad(phi_nulls_deg);

phi = linspace(-pi, pi, 4000);    % Angular range for plotting

%% 2. Calculation of Taylor Phase Mode Weights
% Generates Taylor window weights for the phase modes
w_taylor = taylorwin(M, nbar, -SLL_dB);

% Basic mode weights steered towards the main lobe
w_PM_basic = w_taylor .* exp(-1j * m_indices * phi_0);

%% 3. Null Steering (Projection Matrix Method)
% Steering Matrix: Contains phase mode responses at null directions
A_null = zeros(length(phi_nulls), M);
for k = 1:length(phi_nulls)
    A_null(k, :) = exp(1j * m_indices * phi_nulls(k));
end

% Projection Matrix: C = I - A'(AA')^-1 A
% This forces the weight vector into the subspace orthogonal to the null directions
I = eye(M);
C = I - A_null' * ((A_null * A_null') \ A_null);

% New phase mode weights with nulls incorporated
w_PM_steered = C * w_PM_basic;

%% 4. Physical Antenna Weights (Synthesis)
w_physical = zeros(N, 1);
phi_n = 2 * pi * (0:N-1)' / N;    % Physical positions of elements
for n = 1:N
    sum_val = 0;
    for i = 1:M
        m = m_indices(i);
        bessel_term = besselj(m, circumference);
        
        % Stability check to avoid division by zero
        if abs(bessel_term) < 1e-6, bessel_term = sign(bessel_term)*1e-6; end
        
        % UCA Synthesis Formula
        term = (w_PM_steered(i) / (N * (1j^m) * bessel_term)) * exp(1j * m * phi_n(n));
        sum_val = sum_val + term;
    end
    w_physical(n) = sum_val;
end
w_physical = w_physical / max(abs(w_physical)); % Normalization

%% 5. Analysis and Pattern Calculation
pattern = zeros(size(phi));
for i = 1:M
    m = m_indices(i);
    pattern = pattern + (w_PM_steered(i) * exp(1j * m * phi));
end

% Convert to decibels (dB)
pattern_linear = abs(pattern) / max(abs(pattern));
pattern_db = 20*log10(pattern_linear);

% HPBW Calculation (Half Power Beam Width)
[~, main_lob_idx] = max(pattern_db);
idx_left = find(pattern_db(1:main_lob_idx) <= -3, 1, 'last');
idx_right = find(pattern_db(main_lob_idx:end) <= -3, 1, 'first') + main_lob_idx - 1;
hpbw = abs(rad2deg(phi(idx_right) - phi(idx_left)));

%% 6. Visualization
figure('Name', 'UCA Taylor + Null Steering', 'Color', 'w', 'Position', [100, 100, 1200, 500]);

% Cartesian Plot
subplot(1,2,1);
plot(rad2deg(phi), pattern_db, 'LineWidth', 2); 
grid on; hold on;
xline(phi_0_deg, '--g', 'Main Lobe', 'LineWidth', 1.5);
scatter(phi_nulls_deg, interp1(rad2deg(phi), pattern_db, phi_nulls_deg), 80, 'rx', 'LineWidth', 2, 'DisplayName', 'Null Points');
yline(-SLL_dB, 'r--', 'Target SLL');
ylim([-80 5]); xlim([-180 180]);
title(['Taylor + Null Steering (HPBW: ' num2str(hpbw, '%.1f') '^\circ)']);
xlabel('Angle (Degrees)'); ylabel('Gain (dB)');
legend('Radiation Pattern', 'Main Lobe', 'Null Points', 'Location', 'South');

% Polar Plot
subplot(1,2,2);
db_limit = -60; 
polar_data = pattern_db;
polar_data(polar_data < db_limit) = db_limit; % Clip for visibility
polarplot(phi, polar_data, 'LineWidth', 2);
ax = gca; ax.RLim = [db_limit 0];
title('UCA Polar Pattern (dB)');

% Print Results to Console
fprintf('--- Physical Element Weights ---\n');
disp(table((1:N)', abs(w_physical), rad2deg(angle(w_physical)), ...
    'VariableNames', {'Element_No', 'Magnitude', 'Phase_Degrees'}));