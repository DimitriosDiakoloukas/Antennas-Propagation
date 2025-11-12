clc;
clear;
close all;

% 24x12 planar array of vertical half-wavelength dipoles
% End-fire configuration (maximum radiation along +x axis)

% Antenna parameters
Mx = 24;                 % elements along x-axis
Mz = 12;                 % elements along z-axis
lambda = 1;              % normalization
dx = lambda/2;           
dz = lambda/2;
k = 2*pi/lambda;         % wavenumber

% Steering direction: end-fire along +x
theta_m = 90;            % degrees
phi_m   = 0;             % degrees

% Element positions (centered)
mx = (0:Mx-1) - (Mx-1)/2;
mz = (0:Mz-1) - (Mz-1)/2;
xpos = mx * dx;
zpos = mz * dz;

% Angle grids (1 degree steps for simple Riemann integration)
theta_deg = 1:1:179;
phi_deg   = 0:1:359;

% Compute array factor and element pattern
AF = planar_AF(deg2rad(theta_deg), deg2rad(phi_deg), xpos, zpos, k, theta_m, phi_m);
Fe = elem_halfwave_z(deg2rad(theta_deg));
Fe = Fe(:) * ones(1, numel(phi_deg));

% Total field pattern (normalized)
E = abs(Fe .* AF);
E = E / max(E, [], "all");

% 3D radiation pattern 
[TH, PH] = ndgrid(deg2rad(theta_deg), deg2rad(phi_deg));
R = E;
X = R .* sin(TH).*cos(PH);
Y = R .* sin(TH).*sin(PH);
Z = R .* cos(TH);

figure('Color','w');
surf(X, Y, Z, E, 'EdgeColor','none');
axis equal off; view(45,25);
colormap parula; camlight headlight; lighting gouraud;
title('3D Radiation Pattern (End-fire along +x)');
colorbar;

% Directivity (analytical & computational)
% Analytical
% Αpproximate D = N * D_elem
N = Mx * Mz;
D_elem = 1.64;                  % λ/2 vertical dipole
D_analytical = N * D_elem;

% Computational
% Riemann sum over sphere
U = E.^2;
th = deg2rad(theta_deg);
ph = deg2rad(phi_deg);
dth = deg2rad(1);
dph = deg2rad(1);
Ptot = sum(U .* (sin(th(:)) * ones(1, numel(ph))), 'all') * dth * dph;
D_computed = 4*pi / Ptot;

fprintf('End-fire along +x:\n');
fprintf('Analytical directivity:    D = %.3f  (%.2f dBi)\n', D_analytical, 10*log10(D_analytical));
fprintf('Computed directivity (Riemann): D = %.3f  (%.2f dBi)\n', D_computed, 10*log10(D_computed));

% Helper functions 
function AF = planar_AF(theta, phi, xpos, zpos, k, theta_m_deg, phi_m_deg)
    thm = deg2rad(theta_m_deg); phm = deg2rad(phi_m_deg);
    theta = theta(:).'; phi = phi(:).';
    [TH, PH] = ndgrid(theta, phi);
    [X, Z] = ndgrid(xpos, zpos);
    psi = -k*(X*sin(thm)*cos(phm) + Z*cos(thm));
    AF = zeros(size(TH));
    for ix = 1:numel(xpos)
        for iz = 1:numel(zpos)
            phs = k*(xpos(ix)*sin(TH).*cos(PH) + zpos(iz)*cos(TH)) + psi(ix,iz);
            AF = AF + exp(1j*phs);
        end
    end
end

function Fe = elem_halfwave_z(theta)
    eps = 1e-9;
    theta = max(min(theta, pi-eps), eps);
    Fe = abs(cos((pi/2)*cos(theta))./sin(theta));
    Fe = Fe / max(Fe);
end
