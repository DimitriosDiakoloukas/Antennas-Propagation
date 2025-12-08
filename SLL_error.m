function [error_value] = SLL_error(p, dBTarget)
    % Parameters
    N = 10;
    lambda = 1;
    d = lambda/2;
    k = 2*pi/lambda;
    delta = 0;

    % Symmetric current excitation
    I = [1, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1];

    % Target SLL (linear scale)
    SLL_target_lin = 10^(dBTarget/20);

    % Compute array factor
    theta_deg = linspace(0, 90, 91);
    psi = k*d*cosd(theta_deg) + delta;

    AF = zeros(1, length(theta_deg));
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi);
    end

    % Normalize magnitude
    AF_norm = abs(AF) / max(abs(AF));

    % Find all lobes
    [peaks, ~] = findpeaks(AF_norm);

    % Remove main lobe (global maximum)
    [max_peak, max_idx] = max(peaks);
    if max_peak > 0.98
        peaks(max_idx) = [];
    end

    % Compute mean squared error vs target SLL
    if isempty(peaks)
        error_value = 1e6;                     % penalize if no side lobes found
    else
        squared_error = (peaks - SLL_target_lin).^2;
        error_value = mean(squared_error);
    end

end
