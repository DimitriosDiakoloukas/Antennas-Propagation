clc;
clear;
close all;

% Parameters
f = 1e9;                      % Operating frequency
c = physconst('light');       % Speed of light
lambda = c / f;               % Wavelength
d = lambda / 2;               % Element spacing (half-wavelength)
k = 2 * pi / lambda;          % Wavenumber

% Array size
Nx = 24;                      % Number of elements along x-axis
Nz = 12;                      % Number of elements along z-axis

% Steering directions (in radians)
thetaSteerVec = deg2rad([90, 60]);       % Elevation steering angles
phiSteerVec   = deg2rad([90, 60, 30]);   % Azimuth steering angles

% Angular grids for 3D pattern
thetaGrid = linspace(1e-3, pi - 1e-3, 180);   % Avoid 0 and pi to prevent singularities
phiGrid   = linspace(0, 2*pi, 360);

[THETA, PHI] = meshgrid(thetaGrid, phiGrid);

% Array factor function for a linear array
get_AF = @(N, psi) abs( (abs(psi) < 1e-9).*N + (abs(psi) >= 1e-9).* (sin(N*psi/2) ./ sin(psi/2)) );

% Loop over all steering directions and plot 3D patterns
for thetaSteer = thetaSteerVec
    for phiSteer = phiSteerVec
        
        % Phase shifts for beam steering (depend on the target direction)
        delta_x = -k * d * sin(thetaSteer) * cos(phiSteer);
        delta_z = -k * d * cos(thetaSteer);
        
        % Progressive phase in x and z for the 3D pattern
        psi_x = k * d * sin(THETA) .* cos(PHI) + delta_x;
        psi_z = k * d * cos(THETA) + delta_z;
        
        AF_x = get_AF(Nx, psi_x);
        AF_z = get_AF(Nz, psi_z);
        
        % Element factor: EF = cos(pi/2 * cos(theta)) / sin(theta)
        EF = abs( cos((pi/2) * cos(THETA)) ./ sin(THETA) );
        EF(abs(sin(THETA)) < 1e-5) = 0;    % Avoid numerical issues near theta = 0, pi
        
        % Total field (magnitude)
        E_total = AF_z .* EF .* AF_x;
        E_norm  = E_total / max(E_total(:));   % Normalize to 1
        
        % Convert to Cartesian coordinates for 3D plotting
        X = E_norm .* sin(THETA) .* cos(PHI);
        Y = E_norm .* sin(THETA) .* sin(PHI);
        Z = E_norm .* cos(THETA);
        
        % 3D surface plot
        figure;
        surf(X, Y, Z, E_norm, 'EdgeColor', 'none'); 
        
        axis equal;
        axis vis3d;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        colormap(jet);
        colorbar;
        grid on;
        % view(45, 30);   % Add this to switch angle view if you want
              
        title(sprintf(['3D Radiation Pattern (Array Factor)\n', ...
                       'Steering: \\theta_m = %.0f^\\circ, \\phi_m = %.0f^\\circ'], ...
                       rad2deg(thetaSteer), rad2deg(phiSteer)));
    end
end
