clc;
clear;
close all;

% Parameters
f = 1e9;                            % Operating frequency
c = physconst('light');             % Speed of light
lambda = c / f;                     % Wavelength
d = lambda / 2;                     % Element spacing (half-wavelength)
k = 2 * pi / lambda;                % Wavenumber

Nx = 24;                            % Elements along x-axis
Nz = 12;                            % Elements along z-axis

% 1D array-factor function (uniform linear array)
arrayFactor1D = @(N, psi) abs( ...
    (abs(psi) < 1e-9).*N + ...
    (abs(psi) >= 1e-9).* (sin(N*psi/2) ./ sin(psi/2)) );

% Angular grids for 3D patterns (and numerical integration)
thetaGrid = linspace(1e-3, pi-1e-3, 180);   % Elevation (avoid 0 and pi)
phiGrid   = linspace(0, 2*pi, 360);         % Azimuth
[THETA, PHI] = meshgrid(thetaGrid, phiGrid);

stepDeg = 1;                                % Effective step for integration (deg)
stepRad = deg2rad(stepDeg);                 % Effective step (rad)

% (E) --> ====================== ORDINARY END-FIRE ARRAY ======================

% Assume broadside along z-axis and end-fire along +x-axis
phaseZ_endfire = 0;
phaseX_endfire = -k * d;                    % Ordinary end-fire progressive phase

psiX = k * d * sin(THETA) .* cos(PHI) + phaseX_endfire;
psiZ = k * d * cos(THETA) + phaseZ_endfire;

AFx = arrayFactor1D(Nx, psiX);
AFz = arrayFactor1D(Nz, psiZ);

% Element factor: EF(θ) = cos(pi/2 * cosθ) / sinθ
elemFactor = abs( cos((pi/2) * cos(THETA)) ./ sin(THETA) );

E_field = AFz .* elemFactor .* AFx;
E_norm  = E_field / max(E_field(:));

% Convert to Cartesian coordinates for 3D plotting
X = E_norm .* sin(THETA) .* cos(PHI);
Y = E_norm .* sin(THETA) .* sin(PHI);
Z = E_norm .* cos(THETA);

figure;
surf(X, Y, Z, E_norm, 'EdgeColor', 'none');
axis equal;
axis vis3d;
xlabel('X'); ylabel('Y'); zlabel('Z');
colormap(jet);
colorbar;
title('Ordinary End-Fire Array – 3D Radiation Pattern');

% Directivity – analytical estimate (Kraus-type)
Theta_z = deg2rad(48.4 * lambda / (Nz * d));          % Approx. HPBW in z
Theta_x = deg2rad(105.4 * sqrt(lambda / (Nx * d)));   % Approx. HPBW in x

D_analytical_ef = pi^2 ./ (Theta_x * Theta_z);

% Directivity – numerical evaluation from definition
U_rad = E_field.^2;                           % Radiation intensity ∝ |E|^2
Umax_ef = max(U_rad(:));

powerDensityIntegrand = U_rad .* sin(THETA);  % dP = U sinθ dθ dφ
P_radiated_ef = sum(powerDensityIntegrand(:)) * stepRad * stepRad;

D_numerical_ef = 4 * pi * Umax_ef / P_radiated_ef;

D_analytical_ef_dBi = 10*log10(D_analytical_ef);
D_numerical_ef_dBi  = 10*log10(D_numerical_ef);

% (ST) --> ================= HANSEN–WOODYARD END-FIRE ARRAY ====================

% Same geometry, but Hansen–Woodyard additional phase for higher directivity
phaseZ_HW = 0;
phaseX_HW = -k * d - 2.92 / Nx;              % HW extra phase term

psiX = k * d * sin(THETA) .* cos(PHI) + phaseX_HW;
psiZ = k * d * cos(THETA) + phaseZ_HW;

AFx = arrayFactor1D(Nx, psiX);
AFz = arrayFactor1D(Nz, psiZ);

elemFactor = abs( cos((pi/2) * cos(THETA)) ./ sin(THETA) );

E_field = AFz .* elemFactor .* AFx;
E_norm  = E_field / max(E_field(:));

X = E_norm .* sin(THETA) .* cos(PHI);
Y = E_norm .* sin(THETA) .* sin(PHI);
Z = E_norm .* cos(THETA);

figure;
surf(X, Y, Z, E_norm, 'EdgeColor', 'none');
axis equal;
axis vis3d;
xlabel('X'); ylabel('Y'); zlabel('Z');
colormap(jet);
colorbar;
title('Hansen–Woodyard End-Fire Array – 3D Radiation Pattern');

% Directivity – analytical estimate for HW end-fire
Theta_z = deg2rad(48.4 * lambda / (Nz * d));          % Same z HPBW
Theta_x = deg2rad( 2 * acosd(1 - 0.1398 * lambda / (Nx * d)) );

D_analytical_hw = pi^2 ./ (Theta_x * Theta_z);

% Directivity – numerical evaluation from definition
U_rad = E_field.^2;
Umax_hw = max(U_rad(:));

powerDensityIntegrand = U_rad .* sin(THETA);
P_radiated_hw = sum(powerDensityIntegrand(:)) * stepRad * stepRad;

D_numerical_hw = 4 * pi * Umax_hw / P_radiated_hw;

D_analytical_hw_dBi = 10*log10(D_analytical_hw);
D_numerical_hw_dBi  = 10*log10(D_numerical_hw);

% ====================== SUMMARY PRINTOUT ==============================

fprintf('\n================================= PART E - ST SUMMARY =============================================\n');
fprintf('| Array Type               |  D_analytical [dBi] |  D_numerical [dBi] |  Difference Comparison [dB] |\n');
fprintf('-----------------------------------------------------------------------------------------------------\n');
fprintf('| Ordinary end-fire        | %19.2f | %18.2f | %27.2f |\n', ...
        D_analytical_ef_dBi, D_numerical_ef_dBi, ...
        abs(D_numerical_ef_dBi - D_analytical_ef_dBi));
fprintf('| Hansen–Woodyard end-fire | %19.2f | %18.2f | %27.2f |\n', ...
        D_analytical_hw_dBi, D_numerical_hw_dBi, ...
        abs(D_numerical_hw_dBi - D_analytical_hw_dBi));
fprintf('-----------------------------------------------------------------------------------------------------\n');
