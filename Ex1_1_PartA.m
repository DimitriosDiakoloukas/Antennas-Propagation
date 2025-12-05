clc;
clear;
close all;

% Parameters
f = 1e9;
c = physconst('light');
lambda = c / f;
d = lambda / 2;
k = 2 * pi / lambda;

Nx = 24;
Nz = 12;

% Dipole element definition
elemDip = dipole('Length', 0.48*lambda, 'Width', lambda/100);
elemDip.Tilt = 90; 
elemDip.TiltAxis = [1 0 0]; 

% Rectangular array (in XZ-plane)
antArray = rectangularArray('Element', elemDip, ...
                            'Size', [Nz, Nx], ...       
                            'RowSpacing', d, ...        
                            'ColumnSpacing', d);        

% Rotate the whole array so that it lies in the XZ plane
antArray.Tilt = 90; 
antArray.TiltAxis = [1 0 0];

show(antArray)

% Steering target angles (in radians)
thetaSteerVec = deg2rad([90, 60]); 
phiSteerVec   = deg2rad([90, 60, 30]);

% Angle grids
phi_scan   = deg2rad(0:0.1:360);   % For horizontal cuts
theta_scan = deg2rad(0:0.5:180);   % For vertical cuts

theta_cutoff = pi/2;               % Horizontal cut at theta = 90 deg

% Anonymous function for array factor
get_AF = @(N, psi) abs( (abs(psi) < 1e-9).*N + (abs(psi) >= 1e-9).* (sin(N*psi/2) ./ sin(psi/2)) );

% (a) Horizontal radiation patterns - single figure with 6 subplots
figure;
plotIdx = 1;

for theta_steer = thetaSteerVec
    for phi_steer = phiSteerVec
        
        % Phase shifts for steering (depend on the TARGET direction)
        delta_x = -k * d * sin(theta_steer) * cos(phi_steer);
        delta_z = -k * d * cos(theta_steer);
        
        % Progressive phase in x and z for the horizontal cut
        psi_x = k * d * sin(theta_cutoff) * cos(phi_scan) + delta_x;
        psi_z = round(k * d * cos(theta_cutoff)) + delta_z;  % same expression as original code
                
        AF_x = get_AF(Nx, psi_x);
        AF_z = get_AF(Nz, psi_z);
        
        % Element pattern at theta_cutoff
        E0 = abs(cos((pi/2) * cos(theta_cutoff)) ./ sin(theta_cutoff)); 
        
        % Total field for horizontal cut
        E_horizontal = AF_z .* E0 .* AF_x;
        E_horizontal = E_horizontal / max(E_horizontal);
        
        % Subplot for each (theta_steer, phi_steer) pair
        subplot(2, 3, plotIdx);
        polarplot(phi_scan, E_horizontal)

        title(sprintf(['Horizontal Pattern (\\theta = 90^\\circ)\n', ...
                       'Steering: \\theta_m = %.0f^\\circ, \\phi_m = %.0f^\\circ'], ...
                       rad2deg(theta_steer), rad2deg(phi_steer)));

        thetalim([0 360]);
        grid on;
        
        plotIdx = plotIdx + 1;
    end
end

% (a) Vertical radiation patterns - single figure with 6 subplots
figure;
plotIdx = 1;

for theta_steer = thetaSteerVec
    for phi_steer = phiSteerVec

        % Phase shifts for steering (depend on the TARGET direction)
        delta_x = -k * d * sin(theta_steer) * cos(phi_steer);
        delta_z = -k * d * cos(theta_steer);

        % Front side (phi = phi_steer)
        psi_x_front = k * d * sin(theta_scan) * cos(phi_steer) + delta_x;
        psi_z_front = k * d * cos(theta_scan) + delta_z;
        
        AF_x_front = get_AF(Nx, psi_x_front);
        AF_z_front = get_AF(Nz, psi_z_front);
        
        % Element pattern for the vertical cut
        num = cos((pi/2) * cos(theta_scan));
        den = sin(theta_scan);
        E0 = abs(num ./ den);
        E0(abs(den) < 1e-5) = 0;      % Avoid singularity at theta = 0
        
        E_front = AF_z_front .* E0 .* AF_x_front;

        % Back side (phi = phi_steer + 180 deg)
        phi_back = phi_steer + pi;    % +180 degrees
        
        psi_x_back = k * d * sin(theta_scan) * cos(phi_back) + delta_x;
        psi_z_back = k * d * cos(theta_scan) + delta_z;
        
        AF_x_back = get_AF(Nx, psi_x_back);
        AF_z_back = get_AF(Nz, psi_z_back);
        
        E_back = AF_z_back .* E0 .* AF_x_back;

        % Combine front and back to get full 0–360° vertical pattern
        E_total = [E_front, fliplr(E_back)];
        E_norm  = E_total / max(E_total);
           
        theta_plot = [theta_scan, theta_scan + pi];
        
        % Subplot for each (theta_steer, phi_steer) pair
        subplot(2, 3, plotIdx);
        polarplot(theta_plot, E_norm)
        
        title(sprintf(['Vertical Pattern (\\phi-plane = %.0f^\\circ / %.0f^\\circ)\n', ...
                       'Steering: \\theta_m = %.0f^\\circ'], ...
                       rad2deg(phi_steer), rad2deg(phi_steer)+180, rad2deg(theta_steer)));
        
        thetalim([0 360]);
        rlim([0 1]);
        grid on;
        
        ax = gca;
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        
        plotIdx = plotIdx + 1;
    end
end
