clc;
clear;
close all;

% SLL patterns (theta and psi)
figure;

% -20 dB
dBTarget = -20;
p20 = [3.135379927369486; 4.665416378742279; 9.686730524803112; 9.060534176008748];
%weights from second optimization
%p20 = [5.7998; 5.4897; 7.4711; 8.1992];

subplot(3, 2, 1);
SLL_plot(p20, dBTarget);

subplot(3, 2, 2);
plot_AF_psi(p20, dBTarget);

% -30 dB
dBTarget = -30;
p30 = [2.195133541389668; 3.543791802708922; 7.060181489069492; 9.736232397045448];
%weights from second optimization
%p30 = [3.9001; 5.7543; 8.4134; 9.2211];

subplot(3, 2, 3);
SLL_plot(p30, dBTarget);

subplot(3, 2, 4);
plot_AF_psi(p30, dBTarget);

% -40 dB
dBTarget = -40;
p40 = [2.016247175487722; 3.245847814847494; 5.818461419707745; 8.682570442685023];
%weights from second optimization
%p40 = [2.4398; 5.1124; 7.9882; 9.9542];
subplot(3, 2, 5);
SLL_plot(p40, dBTarget);

subplot(3, 2, 6);
plot_AF_psi(p40, dBTarget);


% Side-lobe behavior for different current sets
figure;

subplot(3, 2, 1);
plot_AF_psi_currents([1 2 3 4]);

subplot(3, 2, 3);
plot_AF_psi_currents([1 3 5 7]);

subplot(3, 2, 5);
plot_AF_psi_currents([2 2 3 4]);

subplot(3, 2, 2);
plot_AF_psi_currents([1 1 3 5]);

subplot(3, 2, 4);
plot_AF_psi_currents([2 4 6 8]);

subplot(3, 2, 6);
plot_AF_psi_currents([1 2 1 3]);


% |AF| (in dB) vs psi for a given SLL target
function plot_AF_psi(p, dBTarget)
    % Parameters
    N = 10;
    lambda = 1;
    d = lambda/2;
    k = 2*pi/lambda;

    % Symmetric currents
    I = [1, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1];

    theta_deg = linspace(-180, -90, 901);
    psi = k*d*cosd(theta_deg);

    AF = zeros(1, numel(theta_deg));
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi);
    end

    AF_dB = 20*log10(abs(AF) / max(abs(AF)));

    plot(psi, AF_dB);
    grid on;

    xlabel('Phase variable \psi (rad)');
    ylabel('Array Factor (dB)');
    title(sprintf('Normalized |AF| vs \\psi (Target SLL = %.1f dB)', dBTarget));

    ylim([-60 0]);
end


% |AF| (in dB) vs psi for a given current vector p
function plot_AF_psi_currents(p)
    % Parameters
    N = 10;
    lambda = 1;
    d = lambda/2;
    k = 2*pi/lambda;

    % Symmetric currents
    I = [1, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1];

    theta_deg = linspace(-180, -90, 901);
    psi = k*d*cosd(theta_deg);

    AF = zeros(1, numel(theta_deg));
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi);
    end

    AF_dB = 20*log10(abs(AF) / max(abs(AF)));

    plot(psi, AF_dB);
    grid on;

    xlabel('Phase variable \psi (rad)');
    ylabel('Array Factor (dB)');
    title(sprintf('Normalized |AF| vs \\psi, I = [%g %g %g %g %g %g %g %g %g %g]', I));

    ylim([-60 0]);
end
