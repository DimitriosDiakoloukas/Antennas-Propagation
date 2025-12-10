clc; 
clear;
close all;   

% Parameters
c = 3e8;                       % Speed of light (m/s)
freq = 1e9;                    % Operating frequency (1 GHz)
lambda = c / freq;             % Wavelength (m)
eta = 120 * pi;                % Intrinsic impedance of free space (Ohms)
dipoleLength = lambda / 2;     % Half-wave dipole length (m)
k = 2 * pi / lambda;           % Wave number (rad/m)

% Distance Range Between Dipoles (0 to 3λ)
distNorm = linspace(1e-5, 3, 500);  % Normalized distance d/λ
distMeters = distNorm * lambda;    % Convert to meters

% Compute Auxiliary Phase Arguments u0, u1, u2
u0 = k .* distMeters;
u1 = k .* (sqrt(distMeters.^2 + dipoleLength^2) + dipoleLength);
u2 = k .* (sqrt(distMeters.^2 + dipoleLength^2) - dipoleLength);

% Compute Sine and Cosine Integral Functions Si(x), Ci(x)
Si_u0 = sinint(u0);
Si_u1 = sinint(u1);
Si_u2 = sinint(u2);

Ci_u0 = cosint(u0);
Ci_u1 = cosint(u1);
Ci_u2 = cosint(u2);

% Compute Mutual Impedance Zm = Rm + j Xm for λ/2 Dipoles
prefactor = eta / (4 * pi);   % Simplified factor for λ/2 dipoles

realZm =  prefactor .* (2*Ci_u0 - Ci_u1 - Ci_u2);      % Mutual resistance
imagZm = -prefactor .* (2*Si_u0 - Si_u1 - Si_u2);      % Mutual reactance

Zm = realZm + 1i * imagZm;

% Plot Mutual Impedance Components
figure;

plot(distNorm, realZm, 'LineWidth', 1.2, ...
     'DisplayName', 'Real Part R_{m}(d)');
hold on;

plot(distNorm, imagZm, 'LineWidth', 1.2, ...
     'DisplayName', 'Imaginary Part X_{m}(d)');

grid on;

title('Mutual Impedance of Two Parallel Half-Wave Dipoles');
xlabel('d / \lambda');
ylabel('Mutual Impedance  Z_{m}(d)  [Ω]');

legend('Location', 'NorthEast');

xlim([0 3]);
ylim([-60 100]);
