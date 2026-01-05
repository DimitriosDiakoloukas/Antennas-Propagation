clc; clear; close all;

% ===== Design frequency (Δ-Λ case) =====
c0 = 3e8;
f0 = 400e6;
lambda = c0/f0;                 % 0.75 m

% ===== Helix (axial mode) =====
N  = 10;
C  = lambda;
D  = C/pi;
R  = D/2;
S  = lambda/4;

% ===== Wire diameter (assignment ~lambda/100) =====
d = lambda/100;
a = d/2;                        % NEC uses radius

% ===== Feed gap =====
gap = lambda/40;

% ===== Ground disk radius (λ/2) =====
Rg  = lambda/2;

% ===== Spider-web ground (8 radials + 4 rings) =====
Nrad = 8;
Ncir = 4;
phis  = (0:Nrad-1) * (2*pi/Nrad);
ringR = (1:Ncir) * (Rg/Ncir);   % Rg/4, Rg/2, 3Rg/4, Rg

% ===== Helix segmentation: 17 seg/turn =====
seg_per_turn = 17;
seg_helix = seg_per_turn * N;
t = linspace(0, 2*pi*N, seg_helix+1);

hx = R*cos(t);
hy = R*sin(t);
hz = gap + (S/(2*pi))*t;

% exact start point
hx(1) = R; hy(1) = 0; hz(1) = gap;

% ===== Feed wire tag =====
feed_tag = 100;
seg_feed = 1;       % 1 segment feed wire
feed_seg = 1;       % excite segment 1

% ===== Frequencies for question 2γ =====
freqs_MHz = [0.2 0.7 1.0 1.3 2.0 3.0] * (f0/1e6);   % [80, 280, 400, 520, 800, 1200]
labels    = ["080MHz","280MHz","400MHz","520MHz","800MHz","1200MHz"];

% ===== Pattern sampling (simple & stable) =====
% Vertical cut: phi=0, theta 0..180 step 5 deg => 37 points
Ntheta_v = 37; theta0_v = 0; dtheta_v = 5;
Nphi_v   = 1;  phi0_v   = 0; dphi_v   = 0;

% 3D: theta 0..180 step 5 , phi 0..360 step 5 
Ntheta_3d = 37; theta0_3d = 0; dtheta_3d = 5;
Nphi_3d   = 73; phi0_3d   = 0; dphi_3d   = 5;

% ===== Produce 6 NEC files =====
for i = 1:numel(freqs_MHz)
    fn = sprintf('helix_pat_%s.nec', labels(i));
    fid = fopen(fn,'w');

    fprintf(fid,'CE Helix patterns for f=%.1f MHz (Question 2g)\n', freqs_MHz(i));

    % --- Geometry (same every time) ---
    write_helix_and_ground( ...
        fid, hx,hy,hz, seg_helix, ...
        feed_tag,seg_feed, R,gap, ...
        Nrad,Ncir, phis, ringR, a);

    % --- Frequency (single) ---
    fprintf(fid,'FR 0 1 0 0 %.3f 0\n', freqs_MHz(i));

    % --- Excitation on feed wire ---
    fprintf(fid,'EX 0 %d %d 0 1 0\n', feed_tag, feed_seg);

    % --- Radiation patterns ---
    % Vertical plane containing z-axis: phi=0, theta sweep
    fprintf(fid,'RP 0 %d %d 1000 %g %g %g %g 0 0\n', ...
        Ntheta_v, Nphi_v, theta0_v, phi0_v, dtheta_v, dphi_v);

    % 3D pattern
    fprintf(fid,'RP 0 %d %d 1000 %g %g %g %g 0 0\n', ...
        Ntheta_3d, Nphi_3d, theta0_3d, phi0_3d, dtheta_3d, dphi_3d);

    fprintf(fid,'EN\n');
    fclose(fid);

    fprintf('OK: created %s\n', fn);
end

disp('Done. Open each helix_pat_XXX.nec in 4nec2 and view Far-field (vertical + 3D).');

% =====================================================================
% Helper: writes geometry with node-correct spider ground (no crossings)
% =====================================================================
function write_helix_and_ground(fid, hx,hy,hz, seg_helix, feed_tag,seg_feed, R,gap, Nrad,Ncir, phis, ringR, a)

    tag = 1;

    % --- HELIX polyline (skip tag 100) ---
    for k = 1:seg_helix
        if tag == feed_tag
            tag = tag + 1;
        end
        fprintf(fid,'GW %d 1 %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
            tag, hx(k), hy(k), hz(k), hx(k+1), hy(k+1), hz(k+1), a);
        tag = tag + 1;
    end

    % --- FEED wire: (R,0,0) -> (R,0,gap) ---
    fprintf(fid,'GW %d %d %.6f %.6f 0 %.6f %.6f %.6f %.6f\n', ...
        feed_tag, seg_feed, R, 0,  R, 0, gap, a);

    % --- GROUND radials segmented at ring radii (and at r=R on phi=0) ---
    for i = 1:Nrad
        phi = phis(i);

        rr = [0, ringR];               % must include ring nodes
        if i == 1 && ~any(abs(rr - R) < 1e-12)
            rr = sort([rr, R]);        % ensure a node exactly at feed radius
        end

        for k = 1:(numel(rr)-1)
            r1 = rr(k);   r2 = rr(k+1);
            x1 = r1*cos(phi); y1 = r1*sin(phi);
            x2 = r2*cos(phi); y2 = r2*sin(phi);

            fprintf(fid,'GW %d 1 %.6f %.6f 0 %.6f %.6f 0 %.6f\n', ...
                tag, x1, y1, x2, y2, a);
            tag = tag + 1;
        end
    end

    % --- Rings as Nrad-sided polygons at each ring radius ---
    for kc = 1:Ncir
        Rc = ringR(kc);
        for i = 1:Nrad
            phi1 = phis(i);
            phi2 = phis(mod(i,Nrad)+1);

            x1 = Rc*cos(phi1); y1 = Rc*sin(phi1);
            x2 = Rc*cos(phi2); y2 = Rc*sin(phi2);

            fprintf(fid,'GW %d 1 %.6f %.6f 0 %.6f %.6f 0 %.6f\n', ...
                tag, x1, y1, x2, y2, a);
            tag = tag + 1;
        end
    end

    fprintf(fid,'GE 0\n');   % free space
end
