clc;
clear;
close all;

% Parameters
f = 1e9;                        % Operating frequency
c = physconst('light');         % Speed of light
lambda = c / f;                 % Wavelength
d = lambda / 2;                 % Element spacing (half-wavelength)
k = 2 * pi / lambda;            % Wavenumber

Nx = 24;                        % Number of elements along x-axis
Nz = 12;                        % Number of elements along z-axis

% Steering angles (same as in pattern calculations)
theta_m_angles = deg2rad([90, 60]); 
phi_m_angles   = deg2rad([90, 60, 30]);

% Array factor function for a uniform linear array
get_AF = @(N, psi) abs( (abs(psi)<1e-9).*N + (abs(psi)>=1e-9).*(sin(N*psi/2)./sin(psi/2)) );

% Angular grid for numerical integration (Riemann sum over the sphere)
d_deg = 1;                         % Angular step in degrees
d_rad = deg2rad(d_deg);            % Angular step in radians

theta_vec = deg2rad(0.5 : d_deg : 179.5);   % Avoid exactly 0° and 180° (singularities)
phi_vec   = deg2rad(0 : d_deg : 360);       % Full azimuth range

[THETA, PHI] = meshgrid(theta_vec, phi_vec);

% Directivity – Methods (c) and (d)
tol = 1e-2;   % Tolerance for comparing floating-point angles

fprintf('\n===================== DIRECTIVITY SUMMARY (c) & (d) =====================\n');
fprintf('| theta_m | phi_m |   (1) Scaling   |    (2) HPBW    |  (3) Integral  |Difference Comparison (3) - (2)|\n');
fprintf('|  [deg]  | [deg] |      [dBi]      |      [dBi]     |      [dBi]     |                [dB]           |\n');
fprintf('-------------------------------------------------------------------------------------------------------\n');

for theta_m = theta_m_angles
    for phi_m = phi_m_angles
        
        % Coordinate transformation due to axis relabeling x'->z, y'->x, z'->y
        [theta_dot, phi_dot] = get_transformed_angles(theta_m, phi_m);
        
        % HPBW constants in radians (taken from analytical / plot data)
        % for x_dot-axis (corresponds to z-axis)
        HPBW90_x_dot   = deg2rad(9.0);    % broadside
        HPBW60_x_dot   = deg2rad(10.0);
        HPBW40_9_x_dot = deg2rad(15.0);
        HPBW56_3_x_dot = deg2rad(11.5);
        
        % for y_dot-axis (corresponds to x-axis)
        HPBW90_y_dot   = deg2rad(5.0);    % broadside
        HPBW60_y_dot   = deg2rad(6.0);
        HPBW30_y_dot   = deg2rad(10.0);
        HPBW48_6_y_dot = deg2rad(6.4);
        HPBW25_7_y_dot = deg2rad(13.0);

        % (1) Directivity from effective aperture scaling

        if abs(theta_dot - 0) < tol && abs(phi_dot - pi/2) < tol
            % Broadside on x'_z' plane
            Dx_dot = 2 * Nz * d / lambda;
            Dy_dot = 2 * Nx * d / lambda;
        elseif  abs(theta_dot - deg2rad(30)) < tol && abs(phi_dot - pi/2) < tol
            Dx_dot = 2 * Nz * d / lambda;
            Dy_dot = 2 * Nx * d / lambda * (HPBW90_y_dot / HPBW60_y_dot);
        elseif  abs(theta_dot - deg2rad(60)) < tol && abs(phi_dot - pi/2) < tol
            Dx_dot = 2 * Nz * d / lambda;
            Dy_dot = 2 * Nx * d / lambda * (HPBW90_y_dot / HPBW30_y_dot);
        elseif  abs(theta_dot - deg2rad(30)) < tol && abs(phi_dot - 0) < tol    
            Dx_dot = 2 * Nz * d / lambda * (HPBW90_x_dot / HPBW60_x_dot);
            Dy_dot = 2 * Nx * d / lambda;
        elseif  abs(theta_dot - deg2rad(41.4)) < tol && abs(phi_dot - deg2rad(40.9)) < tol 
            Dx_dot = 2 * Nz * d / lambda * (HPBW90_x_dot / HPBW40_9_x_dot);
            Dy_dot = 2 * Nx * d / lambda * (HPBW90_y_dot / HPBW48_6_y_dot);
        elseif  abs(theta_dot - deg2rad(64.3)) < tol && abs(phi_dot - deg2rad(56.3)) < tol 
            Dx_dot = 2 * Nz * d / lambda * (HPBW90_x_dot / HPBW56_3_x_dot);
            Dy_dot = 2 * Nx * d / lambda * (HPBW90_y_dot / HPBW25_7_y_dot);
        end

        D1 = pi * cos(theta_dot) * Dx_dot * Dy_dot;
        D1_dBi = 10*log10(D1);

        % (2) Directivity via HPBW approximation (Kraus-type formula)

        if abs(theta_dot - 0) < tol && abs(phi_dot - pi/2) < tol
            Theta_x_dot = deg2rad(48.4 * lambda / (Nz * d));
            Theta_y_dot = deg2rad(48.4 * lambda / (Nx * d));
        elseif  abs(theta_dot - deg2rad(30)) < tol && abs(phi_dot - pi/2) < tol
            Theta_x_dot = deg2rad(48.4 * lambda / (Nz * d));
            Theta_y_dot = HPBW60_y_dot;
        elseif  abs(theta_dot - deg2rad(60)) < tol && abs(phi_dot - pi/2) < tol
            Theta_x_dot = deg2rad(48.4 * lambda / (Nz * d));
            Theta_y_dot = HPBW30_y_dot;
        elseif  abs(theta_dot - deg2rad(30)) < tol && abs(phi_dot - 0) < tol    
            Theta_x_dot = HPBW60_x_dot;
            Theta_y_dot = deg2rad(48.4 * lambda / (Nx * d));
        elseif  abs(theta_dot - deg2rad(41.4)) < tol && abs(phi_dot - deg2rad(40.9)) < tol 
            Theta_x_dot = HPBW40_9_x_dot;
            Theta_y_dot = HPBW48_6_y_dot;
        elseif  abs(theta_dot - deg2rad(64.3)) < tol && abs(phi_dot - deg2rad(56.3)) < tol 
            Theta_x_dot = HPBW56_3_x_dot;
            Theta_y_dot = HPBW25_7_y_dot;
        end

        Theta_h_dot = 1 / ( cos(theta_dot) * sqrt( (cos(phi_dot)^2 / Theta_x_dot^2) + (sin(phi_dot)^2 / Theta_y_dot^2) ) );
        Psi_h_dot   = 1 / ( sqrt( (sin(phi_dot)^2 / Theta_x_dot^2) + (cos(phi_dot)^2 / Theta_y_dot^2) ) );

        D2 = pi^2 / (Theta_h_dot * Psi_h_dot);
        D2_dBi = 10*log10(D2);

        % (3) Directivity from definition (numerical double integral)

        % Phase shifts for steering towards (theta_m, phi_m)
        delta_x = -k * d * sin(theta_m) * cos(phi_m);
        delta_z = -k * d * cos(theta_m);

        Psi_x = k * d * sin(THETA) .* cos(PHI) + delta_x;
        Psi_z = k * d * cos(THETA) + delta_z;

        AF_x = get_AF(Nx, Psi_x);
        AF_z = get_AF(Nz, Psi_z);

        % Element factor: EF = cos(pi/2 * cos(theta)) / sin(theta)
        EF = abs( cos((pi/2) * cos(THETA)) ./ sin(THETA) );

        E_total = AF_z .* EF .* AF_x;

        % Radiation intensity proportional to |E|^2
        U = E_total.^2;
        U_max = max(max(U));

        % Numerical integration over the sphere (Riemann sum)
        integrand = U .* sin(THETA);
        P_rad = sum(sum(integrand)) * d_rad * d_rad;

        % Directivity from definition 
        D3 = 4 * pi * U_max / P_rad;
        D3_dBi = 10*log10(D3);
        
        diff_val = abs(D3_dBi - D2_dBi);
        fprintf('| %7.0f | %5.0f | %15.2f | %14.2f | %14.2f |                        %6.2f |\n', ...
            rad2deg(theta_m), rad2deg(phi_m), D1_dBi, D2_dBi, D3_dBi, diff_val);

    end
end
fprintf('-------------------------------------------------------------------------------------------------------\n');

% Helper function: coordinate transformation
function [theta_dot, phi_dot] = get_transformed_angles(theta_m, phi_m)

    % Original direction cosines (old coordinate system)
    x_old = sin(theta_m) .* cos(phi_m);
    y_old = sin(theta_m) .* sin(phi_m);
    z_old = cos(theta_m);

    tol = 1e-10;
    if abs(x_old) < tol, x_old = 0; end
    if abs(y_old) < tol, y_old = 0; end
    if abs(z_old) < tol, z_old = 0; end

    % New axes mapping: x'->z, y'->x, z'->y
    x_new = z_old;
    y_new = x_old;
    z_new = y_old;

    % Convert back to spherical coordinates (theta_dot, phi_dot)
    theta_dot = acos(z_new);
    if x_new == 0
        phi_dot = pi / 2;
    else
        phi_dot = atan2(y_new, x_new);
        phi_dot = mod(phi_dot, 2*pi);
    end
end
