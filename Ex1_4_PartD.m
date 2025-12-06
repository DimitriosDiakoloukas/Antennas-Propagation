clc;
clear;
close all;

% Parameters
c        = 3e8;               % Speed of light (m/s)
f        = 1e9;               % Frequency (1 GHz)
lambda   = c / f;             % Wavelength (m)
eta0     = 120 * pi;          % Intrinsic impedance of free space (Ohms)
dipoleL  = lambda / 2;        % Dipole length (λ/2)
k        = 2 * pi / lambda;   % Wave number (rad/m)
Z0       = 50.0 + 0i;         % Characteristic impedance of feed line (50 Ohms)

% 2D Parameter Space: element spacing d and height h
% d/λ in [0, 1], h/λ in [0, 1]
dNormRange = linspace(0.01, 1.0, 50);   % Normalized spacing d/λ
hNormRange = linspace(0.01, 1.0, 50);   % Normalized height h/λ

[D_norm, H_norm] = meshgrid(dNormRange, hNormRange);

dGrid = D_norm * lambda;               % Spacing in meters
hGrid = H_norm * lambda;               % Height in meters

% Self and Mutual Impedance Terms
% Self impedance of isolated half-wave dipole (thin-wire approximation)
Z11 = 73.1 + 1i * 42.5; 
Z22 = Z11;

% Mutual impedance function for parallel dipoles
ZmutFun = @(dist) mutualImpedanceParallel(dist, k, dipoleL, eta0);

% Mutual impedances in the structure
Z12 = ZmutFun(dGrid);
Z21 = Z12;

Z13 = ZmutFun(2 * dGrid);
Z23 = Z12;

Z14 = ZmutFun(sqrt((2 * dGrid).^2 + (2 * hGrid).^2));

Z15 = ZmutFun(sqrt(dGrid.^2 + (2 * hGrid).^2));
Z24 = Z15;
Z26 = Z15;

Z16 = ZmutFun(2 * hGrid);
Z25 = Z16;

% Induced current ratio on parasitic element (EMF relation)
currentRatio = (Z15 - Z12) ./ (Z11 + Z13 - Z14 - Z16);

% Input impedance Z_in(d, h)
Z_in = Z21 .* currentRatio + Z22 + Z23 .* currentRatio ...
       - Z24 .* currentRatio - Z25 - Z26 .* currentRatio; 

% Reflection Coefficient Γ(d, h)
gamma      = (Z_in - Z0) ./ (Z_in + Z0);
gammaMag   = abs(gamma);

% 2D Contour Plot of |Γ|
figure;
contourf(D_norm, H_norm, gammaMag, ...
         [0 0.1 0.2 0.3 0.4 0.5 0.6 0.8 1.0], 'LineWidth', 0.5);
colorbar;
hold on;

% Highlight contour |Γ| = 0.3
contour(D_norm, H_norm, gammaMag, [0.3 0.3], ...
        'k-', 'LineWidth', 2);

title('Contour of |Γ| vs Element Spacing d/\lambda and Height h/\lambda (Z_0 = 50 Ω)');
xlabel('d / \lambda');
ylabel('h / \lambda');
xlim([0 1.0]);
ylim([0 1.0]);
colormap jet;

% 3D Surface Plot of |Γ|
figure;
surf(D_norm, H_norm, gammaMag, gammaMag, 'EdgeColor', 'none');
view(3);
colorbar;

xlabel('d / \lambda');
ylabel('h / \lambda');
zlabel('Reflection Coefficient |Γ|');
title('3D Surface of |Γ| vs d/\lambda and h/\lambda (Z_0 = 50 Ω)');
zlim([0 1.0]);

% Find points where |Γ| < 0.3
[rows, cols] = find(gammaMag < 0.3);  

matchingPoints = [];

for i = 1:length(rows)
    d_val = D_norm(rows(i), cols(i));  % d/λ at matching point
    h_val = H_norm(rows(i), cols(i));  % h/λ at matching point
    matchingPoints = [matchingPoints; d_val, h_val]; %#ok<AGROW>
end

% Results
fprintf('------------------ (d) Results ----------------\n');
fprintf('Reflection coefficient |Γ| vs spacing d and height h\n\n');

minGamma = min(gammaMag(:));
[rowMin, colMin] = find(gammaMag == minGamma, 1);
fprintf('Minimum |Γ| = %.3f at d/λ = %.3f, h/λ = %.3f\n\n', ...
        minGamma, D_norm(rowMin, colMin), H_norm(rowMin, colMin));

fprintf('Matching points where |Γ| < 0.3:\n');
fprintf('  (d/λ,  h/λ)\n');
fprintf('  ----------------\n');

for i = 1:size(matchingPoints, 1)
    fprintf('  %.3f   %.3f\n', matchingPoints(i, 1), matchingPoints(i, 2));
end

fprintf('Total matching points: %d\n', size(matchingPoints, 1));
fprintf('------------------------------------------------\n');

% Mutual Impedance for Parallel λ/2 Dipoles
function Zm = mutualImpedanceParallel(dist, k, dipoleL, eta0)

    % Phase arguments
    u0 = k .* dist;
    u1 = k .* (sqrt(dist.^2 + dipoleL^2) + dipoleL);
    u2 = k .* (sqrt(dist.^2 + dipoleL^2) - dipoleL);

    % Sine and cosine integrals
    Si_u0 = sinint(u0);
    Ci_u0 = cosint(u0);

    Si_u1 = sinint(u1);
    Ci_u1 = cosint(u1);

    Si_u2 = sinint(u2);
    Ci_u2 = cosint(u2);

    % Mutual impedance Zmn = Rmn + j Xmn
    Rmn = (eta0 / (4 * pi)) .* (2 .* Ci_u0 - Ci_u1 - Ci_u2);
    Xmn = -(eta0 / (4 * pi)) .* (2 .* Si_u0 - Si_u1 - Si_u2);

    Zm = Rmn + 1i * Xmn;
end
