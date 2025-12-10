function SLL_plot(p, dBTarget)
    % Parameters
    N = 10;                 
    lambda = 1;             
    d = lambda/2;           
    k = 2*pi/lambda;        
    delta = 0;              
    
    % Symmetric current excitation
    I = [1, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1];

    % Compute array factor
    theta_deg = linspace(0, 90, 901);
    psi = k*d*cosd(theta_deg) + delta;

    AF = zeros(1, length(theta_deg));
    for n = 1:N
        AF = AF + I(n) * exp(1j*(n-1)*psi);
    end
    
    % Normalize AF and convert to dB
    AF_norm = abs(AF) / max(abs(AF));
    AF_dB = 20*log10(AF_norm);

    % Plot array factor
    plot(theta_deg, AF_dB); hold on;

    % Plot target SLL
    yline(dBTarget, '--g');

    % Identify and display side lobes
    sll_threshold = 10^(dBTarget/20);
    [peaks, locs] = findpeaks(AF_norm, 'MinPeakHeight', sll_threshold*1e-8);

    for i = 1:length(peaks)
        peak_angle = theta_deg(locs(i));
        peak_dB = 20*log10(peaks(i));
        if peak_dB < -1   % ignore main lobe
            scatter(peak_angle, peak_dB, 70, 'filled', 'MarkerFaceColor', 'k');
        end
    end

    % Formatting
    grid on;
    xlabel('Angle \theta (degrees)');
    ylabel('Array Factor (dB)');
    title(sprintf('Normalized Radiation Pattern (Target SLL = %.1f dB)', dBTarget));
    legend('Array Factor', 'Target SLL');
    ylim([-60 0]);
    xlim([0 90]);
    hold off;

end
