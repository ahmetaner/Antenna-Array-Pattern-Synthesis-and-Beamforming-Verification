%% UCA Design: Chebyshev Windowing + Multiple Nulls
clear; clc; close all;

%% 1. Parameters
N = 8;                  % Number of antenna elements
R_lambda = (8.4*1e-2)/(12.5*1e-2); % Radius / Wavelength ratio
circumference = 2*pi*R_lambda; 
M = 12;                  % Number of Phase Modes (Odd number recommended)
m_indices = (-(M-1)/2 : (M-1)/2)'; 

% --- Chebyshev and Null Settings ---
SLL_dB = 20;            % Desired Sidelobe Level (dB)
phi_0_deg = 0;          % Main lobe direction (Degrees)
phi_0 = deg2rad(phi_0_deg);       
null_angles_deg = [-120]; % Angles to be nulled (Degrees)
null_angles = deg2rad(null_angles_deg);
phi = linspace(-pi, pi, 4000); 

%% 2. Chebyshev Phase Mode Weights Calculation
% The chebwin function generates coefficients for M modes
w_cheb = chebwin(M, SLL_dB);

% Steering phase modes toward the main lobe direction (phi_0)
% Note: Chebyshev coefficients are real; we add phase for steering
w_PM = w_cheb .* exp(-1j * m_indices * phi_0);

%% 3. Null Steering (Optional: Adding Nulls over Chebyshev)
% If you wish to add nulls without disrupting the Chebyshev structure,
% Constrained Least Squares can use Chebyshev coefficients as the target.
% However, for pure Chebyshev performance, w_PM can be used directly.

%% 4. Physical Antenna Weights (Synthesis)
w_physical = zeros(N, 1);
phi_n = 2 * pi * (0:N-1)' / N; 

for n = 1:N
    sum_val = 0;
    for i = 1:M
        m = m_indices(i);
        bessel_term = besselj(m, circumference);
        
        % Check for stability in the mode cut-off region (if Bessel term is too small)
        if abs(bessel_term) < 1e-6, bessel_term = sign(bessel_term)*1e-6; end
        
        % UCA Synthesis Formula
        term = (w_PM(i) / (N * (1j^m) * bessel_term)) * exp(1j * m * phi_n(n));
        sum_val = sum_val + term;
    end
    w_physical(n) = sum_val;
end

% Normalize Weights (For physical limits)
w_physical = w_physical / max(abs(w_physical));

fprintf('--- Chebyshev Physical Weights ---\n');
disp(table((1:N)', abs(w_physical), rad2deg(angle(w_physical)), ...
    'VariableNames', {'Element_No', 'Amplitude', 'Phase_Degrees'}));

%% 5. Analysis
pattern = zeros(size(phi));
for i = 1:M
    m = m_indices(i);
    % Far-field pattern calculation
    pattern = pattern + (w_PM(i) * exp(1j * m * phi));
end

pattern_linear = abs(pattern) / max(abs(pattern));
pattern_db = 20*log10(pattern_linear);

% --- HPBW Calculation ---
[~, main_lob_idx] = max(pattern_db);
target_hp = -3;
idx_left = find(pattern_db(1:main_lob_idx) <= target_hp, 1, 'last');
idx_right = find(pattern_db(main_lob_idx:end) <= target_hp, 1, 'first') + main_lob_idx - 1;
hpbw = abs(rad2deg(phi(idx_right) - phi(idx_left)));

%% 6. Visualization
figure('Name', 'UCA Chebyshev Analysis', 'Color', 'w', 'Position', [100, 100, 1200, 500]);

% Cartesian Plot
subplot(1,2,1);
plot(rad2deg(phi), pattern_db, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); 
grid on; hold on;
xline(phi_0_deg, '--g', 'Main Lobe', 'LineWidth', 1.5);

for k = 1:length(null_angles_deg)
    xline(null_angles_deg(k), '--r', 'Desired Null', 'LineWidth', 1.2);
end

yline(-SLL_dB, 'k--', ['Target SLL: -' num2str(SLL_dB) 'dB']);
ylim([-SLL_dB-20 5]); xlim([-180 180]);
title(['Chebyshev Pattern (HPBW: ' num2str(hpbw, '%.1f') '^\circ)']);
xlabel('Angle (Degrees)'); ylabel('Gain (dB)');

% Polar Plot
subplot(1,2,2);
db_limit = -60; 
polar_data = pattern_db;
polar_data(polar_data < db_limit) = db_limit; 
polarplot(phi, polar_data, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
hold on;
ax = gca; ax.RLim = [db_limit 0];
title(['UCA Polar Plot (SLL: ' num2str(SLL_dB) ' dB)']);