clc;
clear;
close all;

% 24x12 planar array of vertical half-wavelength dipoles
% Part 1.1(a): Horizontal (vs φ) and Vertical (vs θ) cuts
% Part 1.1(b): Full 3D polar plot of the radiation pattern
% Part 1.1(c): Directivity via full-sphere integration and HPBW approximation
% Part 1.1(d): Directivity via simple 1 degree Riemann sum on a sphere and comparison with (c)
% All patterns normalized to maximum (in dB for 2D plots)


% Antenna Parameters
Mx = 24;                 % number of elements along x-axis (horizontal)
Mz = 12;                 % number of elements along z-axis (vertical)
lambda = 1;              % set λ = 1 for normalization
dx = lambda/2;           % spacing along x
dz = lambda/2;           % spacing along z
k = 2*pi/lambda;         % wavenumber

% Angle grids for (a)
phi_deg   = 0:0.25:360;  % for horizontal pattern (at fixed θ = θm)
theta_deg = 0:0.25:180;  % for vertical pattern (at fixed φ = φm)

% Angle grids for (b) 3D (avoid singularities exactly at 0, π)
theta3D_deg = linspace(0.25, 179.75, 361);
phi3D_deg   = 0:1:360;
dB_floor    = -30;       % clip floor for 3D color/size

% Steering direction combinations (θm, φm)
theta_m_list = [90, 60];       % in degrees
phi_m_list   = [90, 60, 30];   % in degrees

% Element positions (centered around zero)
mx = (0:Mx-1) - (Mx-1)/2;      % index along x
mz = (0:Mz-1) - (Mz-1)/2;      % index along z
xpos = mx * dx;                % physical x positions
zpos = mz * dz;                % physical z positions

fprintf(' (theta_m, phi_m) |   D_integral   |   D_HPBW   |  HPBW_theta  HPBW_phi\n');

% Loop over all requested steering directions
for thm = theta_m_list
    for phm = phi_m_list

        % (a) 2D CUTS
        % Horizontal pattern: θ = θm, varying φ
        th = deg2rad(thm);
        ph = deg2rad(phi_deg);

        % Array factor AF(θm, φ)
        AF_h = planar_AF(th, ph, xpos, zpos, k, thm, phm);

        % Element pattern of a vertical λ/2 dipole (depends only on θ)
        % For horizontal cut (θ = θm), it is constant and replicated
        Fe_h = elem_halfwave_z(th) * ones(size(ph));

        % Total pattern magnitude
        Patt_h = abs(Fe_h .* AF_h);
        Patt_h = Patt_h / max(Patt_h);  % maximum normalization
        Patt_h_dB = 20*log10(max(Patt_h, 1e-6));

        % Vertical pattern: φ = φm, varying θ
        ph = deg2rad(phm);
        th_vec = deg2rad(theta_deg);

        AF_v = planar_AF(th_vec, ph, xpos, zpos, k, thm, phm);
        Fe_v = elem_halfwave_z(th_vec);

        Patt_v = abs(Fe_v .* AF_v);
        Patt_v = Patt_v / max(Patt_v);
        Patt_v_dB = 20*log10(max(Patt_v, 1e-6));

        % Plot (a)
        figure('Color','w','Name',sprintf('2D cuts: theta_m=%d°, phi_m=%d°',thm,phm));
        tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

        % Horizontal (vs φ)
        nexttile;
        plot(phi_deg, Patt_h_dB,'LineWidth',1.4); grid on;
        xlim([0 360]); ylim([-60 0]);
        xlabel('\phi (deg)'); ylabel('Normalized |E| (dB)');
        title(sprintf('Horizontal pattern at \\theta = %d^\\circ', thm));

        % Vertical (vs θ) for φ = φm
        nexttile;
        plot(theta_deg, Patt_v_dB,'LineWidth',1.4); grid on;
        xlim([0 180]); ylim([-60 0]);
        xlabel('\theta (deg)'); ylabel('Normalized |E| (dB)');
        title(sprintf('Vertical pattern at \\phi = %d^\\circ', phm));

        sgtitle(sprintf('24\\times12 \\lambda/2 dipoles, d = \\lambda/2, steer to (\\theta_m, \\phi_m) = (%d^\\circ, %d^\\circ)',...
                        thm, phm));

        % (b) 3D POLAR PATTERN
        % Array factor on full (θ,φ) grid
        AF_3D = planar_AF(deg2rad(theta3D_deg), deg2rad(phi3D_deg), xpos, zpos, k, thm, phm); % Nθ×Nφ

        % Element pattern over θ (replicate over φ)
        Fe_3D = elem_halfwave_z(deg2rad(theta3D_deg)); % 1×Nθ
        Fe_3D = Fe_3D(:) * ones(1, numel(phi3D_deg)); % Nθ×Nφ

        % Total pattern, normalize and convert to dB
        Patt3D = abs(Fe_3D .* AF_3D);
        Patt3D = Patt3D / max(Patt3D, [], "all");
        Patt3D_dB = 20*log10(max(Patt3D, 1e-6));
        Patt3D_dB = max(Patt3D_dB, dB_floor);            % clip floor
        R = 10.^(Patt3D_dB/20);                          % radius for 3D polar

        % Spherical to Cartesian
        [TH, PH] = ndgrid(deg2rad(theta3D_deg), deg2rad(phi3D_deg));
        X = R .* sin(TH).*cos(PH);
        Y = R .* sin(TH).*sin(PH);
        Z = R .* cos(TH);

        % 3D surface
        figure('Color','w','Name',sprintf('3D: theta_m=%d°, phi_m=%d°',thm,phm));
        s = surf(X, Y, Z, Patt3D_dB, 'EdgeColor','none');
        axis equal off; view(45,25); camlight headlight; lighting gouraud;
        colormap parula;
        cb = colorbar; cb.Label.String = 'Normalized |E| (dB)';
        title(sprintf('3D pattern, steer to (\\theta_m,\\phi_m)=(%d^\\circ,%d^\\circ)', thm, phm));

        % (c) DIRECTIVITY (two methods)
        th_rad = deg2rad(theta3D_deg);
        ph_rad = deg2rad(phi3D_deg);
        dth = (th_rad(end) - th_rad(1)) / (numel(th_rad)-1);
        dph = (ph_rad(end) - ph_rad(1)) / (numel(ph_rad)-1);
        U = Patt3D.^2;                                
        S = sin(th_rad(:)) * ones(1, numel(ph_rad));  % sinθ weights
        Ptot = sum(U(:) .* S(:)) * dth * dph;         % total radiated (normalized)
        D_integral = 4*pi / Ptot;

        HPBW_phi   = hpbw_1d(phi_deg,   Patt_h_dB);    % horizontal cut
        HPBW_theta = hpbw_1d(theta_deg, Patt_v_dB);    % vertical cut
        D_HPBW = 41253 / max(HPBW_theta,1e-6) / max(HPBW_phi,1e-6);

        fprintf('     (%3d°, %3d°) |   %10.3f   |   %8.3f |   %8.2f deg   %8.2f deg\n',...
                 thm, phm, D_integral, D_HPBW, HPBW_theta, HPBW_phi);

        % (d) DIRECTIVITY via simple 1° Riemann sum --------
        thetaR_deg = 1:1:179;
        phiR_deg   = 0:1:359;

        AF_R = planar_AF(deg2rad(thetaR_deg), deg2rad(phiR_deg), ...
                         xpos, zpos, k, thm, phm);                     % Nθ×Nφ
        Fe_R = elem_halfwave_z(deg2rad(thetaR_deg));                    % 1×Nθ
        Fe_R = Fe_R(:) * ones(1, numel(phiR_deg));                      % Nθ×Nφ

        E_R = Fe_R .* AF_R;
        U_R = abs(E_R).^2;                                             
        U_R = U_R / max(U_R, [], 'all');    

        thR = deg2rad(thetaR_deg);
        phR = deg2rad(phiR_deg);
        dthR = deg2rad(1);
        dphR = deg2rad(1);

        % Riemann sum over sphere
        Ptot_R = sum( U_R .* (sin(thR(:)) * ones(1, numel(phR))) , 'all') * dthR * dphR;
        D_Riemann1deg = 4*pi / Ptot_R;

        % Comparisons to (c)
        err_int_vs_R  = 100 * abs(D_integral - D_Riemann1deg) / D_integral;
        err_hpbw_vs_R = 100 * abs(D_HPBW     - D_Riemann1deg) / max(D_Riemann1deg,1e-12);

        fprintf('--> Riemann 1 deg |   %10.3f   |  err vs (c-int): %6.3f%%,  vs (c-HPBW): %6.3f%%\n', ...
                D_Riemann1deg, err_int_vs_R, err_hpbw_vs_R);

    end
end


% Helper Functions 
function AF = planar_AF(theta, phi, xpos, zpos, k, theta_m_deg, phi_m_deg)
    % theta, phi may be scalars or vectors (radians).
    % Returns AF with size [numel(theta) × numel(phi)] when both are vectors.

    thm = deg2rad(theta_m_deg);
    phm = deg2rad(phi_m_deg);

    theta = theta(:).';            % 1×Nθ
    phi   = phi(:).';              % 1×Nφ
    [TH, PH] = ndgrid(theta, phi); % Nθ×Nφ

    [X, Z] = ndgrid(xpos, zpos);   % Mx×Mz

    % Steering phase for maximum at (θm, φm)
    psi = -k*( X*sin(thm)*cos(phm) + Z*cos(thm) ); % Mx×Mz

    % Sum over all elements
    AF = zeros(size(TH));
    for ix = 1:numel(xpos)
        for iz = 1:numel(zpos)
            phs = k*( xpos(ix)*sin(TH).*cos(PH) + zpos(iz)*cos(TH) ) + psi(ix,iz);
            AF = AF + exp(1j*phs);
        end
    end

    % If one dimension is singleton, return a row vector
    if size(AF,1)==1 || size(AF,2)==1
        AF = AF(:).';
    end
end

function Fe = elem_halfwave_z(theta)
    % Vertical λ/2 dipole element pattern (z-axis)
    % |F(θ)| = |cos((π/2)cosθ) / sinθ|, 0<θ<π, normalized
    eps = 1e-9;
    theta = max(min(theta, pi-eps), eps);
    Fe = abs( cos( (pi/2)*cos(theta) ) ./ sin(theta) );
    Fe = Fe / max(Fe);
end

function r = deg2rad(d)
    r = d*pi/180;
end

function bw = hpbw_1d(ang_deg, patt_dB)
    % −3 dB beamwidth (degrees) around the main lobe in a 1D cut.
    pmax = max(patt_dB);
    target = pmax - 3;
    [~, imax] = max(patt_dB);
    
    % left crossing
    i1 = imax; 
    while i1>1 && patt_dB(i1) > target
        i1 = i1-1;
    end
    if i1==1
        aL = ang_deg(1);
    else
        aL = lin_interp_at(patt_dB(i1), ang_deg(i1), patt_dB(i1+1), ang_deg(i1+1), target);
    end
    
    % right crossing
    i2 = imax;
    while i2<numel(patt_dB) && patt_dB(i2) > target
        i2 = i2+1;
    end
    if i2==numel(patt_dB)
        aR = ang_deg(end);
    else
        aR = lin_interp_at(patt_dB(i2-1), ang_deg(i2-1), patt_dB(i2), ang_deg(i2), target);
    end
    
    bw = mod(aR - aL, 360);
    if bw > 180
        bw = 360 - bw;
    end
end

function x = lin_interp_at(y1, x1, y2, x2, y)
    % Linear interpolation to find x where y(x)=target between two samples.
    if abs(y2-y1) < 1e-12
        x = (x1+x2)/2;
    else
        t = (y - y1)/(y2 - y1);
        x = x1 + t*(x2 - x1);
    end
end
