clc;
clear;
close all;

% Parameters
c      = 3e8;              % Speed of light (m/s)
f       = 1e9;              % Frequency (1 GHz)
lambda  = c / f;           % Wavelength (m)
eta0    = 120 * pi;         % Intrinsic impedance of free space (Ohms)
dipoleL = lambda / 2;       % Half-wave dipole length (m)
k       = 2 * pi / lambda;  % Wave number (rad/m)
r_obs   = 50;               % Observation distance (m)

% Array Configuration (3 parallel λ/2 dipoles)
elemSpacing  = 3 * lambda / 4;    % Element spacing (d = 0.75λ)
elemSpacing13 = 2 * elemSpacing;  % Distance between element 1 and 3 (2d = 1.5λ)

% Self-impedance of isolated half-wave dipole
Z11_self = 73.1 + 1i * 42.5;      

% Mutual impedances Zmn (Z12, Z13) using mutualImpedance
ZmutFun = @(d_input) mutualImpedance(d_input, k, dipoleL, eta0);

% Z12 = Z(d = 3λ/4)
Z12_mut = ZmutFun(elemSpacing);   

% Z13 = Z(d = 3λ/2)
Z13_mut = ZmutFun(elemSpacing13); 

% Notation a, b, c as in theory
a = Z11_self;  
b = Z12_mut;   
c = Z13_mut;   

% Calculate Current Ratios (I1/I2 = I3/I2)
% Passive element equation (V1 = 0): a*I1 + b*I2 + c*I3 = 0
% Symmetry: I1 = I3  =>  (a + c)*I1 = -b*I2
% Current ratio: I1/I2 = -b / (a + c)
currentRatio = -b / (a + c); 

% Set current of the driven element
I2 = 1.0; 
I1 = currentRatio * I2; 
I3 = I1;

% Input Impedance Z_in of Driven Element
% Z_in = V2 / I2
% V2 = b*I1 + a*I2 + b*I3
% => Z_in = a + 2*b*(I1/I2)
Z_in = a + 2 * b * currentRatio; 

fprintf('------------------- (b) EMF Method Results ----------------\n');
fprintf('Self and mutual impedances (Ω):\n');
fprintf('---->  a = Z11 = %7.2f %+.2fj\n', real(a), imag(a));
fprintf('---->  b = Z12 = %7.2f %+.2fj\n', real(b), imag(b));
fprintf('---->  c = Z13 = %7.2f %+.2fj\n', real(c), imag(c));
fprintf('\n');

fprintf('Current ratio (I1/I2 = I3/I2):\n');
fprintf('  |I1/I2| = %.3f,  angle(I1/I2) = %.1f degrees\n', ...
    abs(currentRatio), rad2deg(angle(currentRatio)));
fprintf('\n');

fprintf('Input impedance of driven dipole:\n');
fprintf('  Z_in = %7.2f %+.2fj  Ω\n', real(Z_in), imag(Z_in));
fprintf('-----------------------------------------------------------\n\n');

% Radiation Pattern in the Horizontal Plane (θ = 90°)
theta_h = pi/2; 
phi_h   = linspace(0, 2*pi, 360);

% Element factor magnitude (shape only; absolute scale not critical after normalization)
E0_mag = abs(eta0 * I2 * k * dipoleL * sin(theta_h) / (4*pi*r_obs));
        
% Array factor (3 elements at positions -d, 0, +d along x-axis)
AF_h = I2 ...
     + I1 .* exp(-1i * k * elemSpacing * cos(phi_h) * sin(theta_h)) ...
     + I3 .* exp( 1i * k * elemSpacing * cos(phi_h) * sin(theta_h));

AF_mag = abs(AF_h);

% Normalize total field
E_h     = AF_mag .* E0_mag;
E_h_norm = E_h / max(E_h);

figure;
polarplot(phi_h, E_h_norm, 'LineWidth', 1.2);
title('Normalized Horizontal Radiation Pattern (d = 0.75\lambda)');
ax = gca;
ax.ThetaDir = 'counterclockwise';
ax.ThetaZeroLocation = 'right'; 

% Mutual Impedance Zmn
function Zm = mutualImpedance(d_input, k, dipoleL, eta0)
    % Phase arguments U0, U1, U2
    u0 = k * d_input;
    u1 = k * (sqrt(d_input^2 + dipoleL^2) + dipoleL);
    u2 = k * (sqrt(d_input^2 + dipoleL^2) - dipoleL);

    % Sine and cosine integrals
    Si_u0 = sinint(u0);
    Si_u1 = sinint(u1);
    Si_u2 = sinint(u2);
    Ci_u0 = cosint(u0);
    Ci_u1 = cosint(u1);
    Ci_u2 = cosint(u2);

    % Mutual impedance Zmn = Rmn + j Xmn
    Rmn = (eta0 / (4 * pi)) * (2 * Ci_u0 - Ci_u1 - Ci_u2);
    Xmn = -(eta0 / (4 * pi)) * (2 * Si_u0 - Si_u1 - Si_u2);

    Zm = Rmn + 1i * Xmn;
end
