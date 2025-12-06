clc;
clear;
close all;

% Parameters
c        = 3e8;              % Speed of light (m/s)
f        = 1e9;              % Frequency (1 GHz)
lambda   = c / f;            % Wavelength (m)
eta0     = 120 * pi;         % Intrinsic impedance of free space (Ohms)
dipoleL  = lambda / 2;       % Dipole length (λ/2)
k        = 2 * pi / lambda;  % Wave number (rad/m)

Z11      = 73.1 + 1i*42.5;   % Self-impedance of isolated λ/2 dipole
Z0       = 50.0 + 0i;        % Characteristic line impedance (50 Ohms)

% Distance Range between Elements (0 <= d <= λ)
dNormRange = linspace(0.001, 1.0, 400);   % Normalized spacing d/λ
dRange     = dNormRange * lambda;         % Physical spacing in meters

% Input Impedance Z_in(d)
% Mutual impedance function (vectorized)
ZmutFun = @(d_input) mutualImpedanceVec(d_input, k, dipoleL, eta0);

% b(d) = Z12(d),   c(d) = Z13(2d)
b = ZmutFun(dRange);          % Mutual impedance between adjacent elements
c = ZmutFun(2 * dRange);      % Mutual impedance between elements spaced 2d

a = Z11;

% Current ratio for parasitic elements: I1/I2 = I3/I2 = -b / (a + c)
currentRatio = -b ./ (a + c);

% Input impedance of driven element: Z_in = a + 2*b*(I1/I2)
Z_in = a + 2 .* b .* currentRatio;

% Reflection Coefficient Γ(d) at 50 Ohms 
gamma     = (Z_in - Z0) ./ (Z_in + Z0);
gammaMag  = abs(gamma);

% Plot |Γ| versus d/λ
figure;
plot(dNormRange, gammaMag, 'LineWidth', 1.2, 'DisplayName', '|Γ(d/λ)|');
hold on;

title('Reflection Coefficient Magnitude |Γ| vs Element Spacing d/λ (Z_0 = 50 Ω)');
xlabel('d / \lambda');
ylabel('Reflection Coefficient |Γ|');
xlim([0 1.0]);
ylim([0 1.0]);
grid on;

% Threshold line |Γ| = 0.3 (matching criterion)
gammaThresh = 0.3;
hThresh = plot(dNormRange, gammaThresh * ones(size(dNormRange)), ...
               'k', 'LineWidth', 2, 'DisplayName', '|Γ| = 0.3');
uistack(hThresh, 'bottom');

legend('Location', 'NorthEast');

% Μatching regions (|Γ| < 0.3)
fprintf('------------------- (c) Results ----------------\n');
fprintf('Reflection coefficient vs element spacing d (0 <= d <= λ)\n\n');

[gammaMin, idxMin] = min(gammaMag);
fprintf('Minimum |Γ| = %.3f at d/λ = %.3f\n\n', gammaMin, dNormRange(idxMin));

fprintf('Regions where |Γ| < 0.3 (good matching to 50 Ω):\n');

matchingMask = gammaMag < gammaThresh;
idx = find(matchingMask);
ranges = {};

% Find contiguous segments where |Γ| < 0.3
if ~isempty(idx)
    startIdx = idx(1);
    for i = 2:length(idx)
        if idx(i) ~= idx(i-1) + 1
            ranges{end+1} = [dNormRange(startIdx), dNormRange(idx(i-1))];
            startIdx = idx(i);
        end
    end
    ranges{end+1} = [dNormRange(startIdx), dNormRange(idx(end))];
end

if isempty(ranges)
    fprintf('  No spacing d/λ in [0, 1] satisfies |Γ| < 0.3.\n');
else
    for r = ranges
        dStart = r{1}(1);
        dEnd   = r{1}(2);
        fprintf('  d/λ from %.3f to %.3f\n', dStart, dEnd);
    end
end
fprintf('------------------------------------------------\n');

% Mutual Impedance Zmn
function Zm = mutualImpedanceVec(d_input, k, dipoleL, eta0)

    % Phase arguments U0, U1, U2
    u0 = k .* d_input;
    u1 = k .* (sqrt(d_input.^2 + dipoleL^2) + dipoleL);
    u2 = k .* (sqrt(d_input.^2 + dipoleL^2) - dipoleL);

    % Sine and cosine integrals
    Si_u0 = sinint(u0);
    Si_u1 = sinint(u1);
    Si_u2 = sinint(u2);

    Ci_u0 = cosint(u0);
    Ci_u1 = cosint(u1);
    Ci_u2 = cosint(u2);

    % Mutual impedance Zmn = Rmn + j Xmn
    Rmn = (eta0 / (4 * pi)) .* (2 .* Ci_u0 - Ci_u1 - Ci_u2);
    Xmn = -(eta0 / (4 * pi)) .* (2 .* Si_u0 - Si_u1 - Si_u2);

    Zm = Rmn + 1i * Xmn;
end
