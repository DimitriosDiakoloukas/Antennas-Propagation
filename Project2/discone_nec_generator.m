clc; clear; close all;

% ====== Δ-Λ case από εκφώνηση/αναφορά ======
lambda = 0.6;                 % m  -> f0 = 500 MHz
Rd     = 0.34*lambda;         % disk radius
Lc     = 0.5*lambda;          % cone wire length
theta  = 30*pi/180;           % half-angle (2θ0=60°)
gap    = lambda/20;           % λ/20 κάτω το cone junction
d      = lambda/200;          % wire diameter
a      = d/2;                 % NEC uses radius

% Cone geometry
h     = Lc*cos(theta);
Rcone = Lc*sin(theta);

% ====== Segmentation (όπως αιτιολογεί: ΔL≈1.58 cm στον κώνο στα 2 GHz) ======
% Lc = 0.30 m -> seg_c = 19 => ΔL = 0.30/19 = 0.015789 m ≈ 1.58 cm
seg_c = 19;
seg_d = 14;

% ====== Frequency sweep για (α): 0.5f0 έως 4f0, step 10 MHz ======
f_start = 250;  % MHz
f_stop  = 2000; % MHz
f_step  = 10;   % MHz
n_steps = round((f_stop - f_start)/f_step) + 1;

% ====== Junctions όπως λέει η εκφώνηση ======
% Disk junction at z=0
z_disk_junction = 0;

% Cone junction at z = -gap (λ/20 κάτω από το disk junction)
z_cone_junction = -gap;

% Cone skirt plane is further down by cone height h
z_cone_skirt = z_cone_junction - h;

% ====== Feed wire: μεταξύ των δύο junctions, ΜΕ 1 segment ======
seg_feed = 1;
feed_tag = 17;     % μετά τα 16 σύρματα (8 disk + 8 cone)
feed_seg = 1;      % αφού έχει 1 segment

% ----------------------------
% 1) discone_geometry.nec @ 500 MHz (για F3 screenshot)
% ----------------------------
fid = fopen('discone_geometry.nec','w');
fprintf(fid,'CE Discone (exact per assignment: cone junction at -lambda/20, 1-seg feed)\n');

tag = 1;

% Disk: 8 radials from junction (0,0,0) to (Rd*cosφ, Rd*sinφ, 0)
for i = 0:7
    phi = i*2*pi/8;
    x2  = Rd*cos(phi);
    y2  = Rd*sin(phi);

    fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        tag, seg_d, ...
        0, 0, z_disk_junction, ...
        x2, y2, z_disk_junction, ...
        a);
    tag = tag + 1;
end

% Cone: 8 wires start at common cone junction (0,0,-gap) and end on skirt circle
for i = 0:7
    phi = i*2*pi/8;
    x2  = Rcone*cos(phi);
    y2  = Rcone*sin(phi);

    fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        tag, seg_c, ...
        0, 0, z_cone_junction, ...
        x2, y2, z_cone_skirt, ...
        a);
    tag = tag + 1;
end

% Feed wire: one segment from disk junction to cone junction
fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
    feed_tag, seg_feed, ...
    0, 0, z_disk_junction, ...
    0, 0, z_cone_junction, ...
    a);

fprintf(fid,'GE 0\n');
fprintf(fid,'FR 0 1 0 0 500 0\n');               % για να γράφει 500 MHz στο F3
fprintf(fid,'EX 0 %d %d 0 1 0\n', feed_tag, feed_seg); % EX στο 1-seg feed wire
fprintf(fid,'EN\n');
fclose(fid);

% ----------------------------
% 2) discone_sweep.nec (για Ερώτημα α)
% ----------------------------
fid = fopen('discone_sweep.nec','w');
fprintf(fid,'CE Discone SWEEP (250-2000 MHz step 10 MHz) exact per assignment\n');

tag = 1;

% Disk
for i = 0:7
    phi = i*2*pi/8;
    x2  = Rd*cos(phi);
    y2  = Rd*sin(phi);

    fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        tag, seg_d, ...
        0, 0, z_disk_junction, ...
        x2, y2, z_disk_junction, ...
        a);
    tag = tag + 1;
end

% Cone
for i = 0:7
    phi = i*2*pi/8;
    x2  = Rcone*cos(phi);
    y2  = Rcone*sin(phi);

    fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        tag, seg_c, ...
        0, 0, z_cone_junction, ...
        x2, y2, z_cone_skirt, ...
        a);
    tag = tag + 1;
end

% Feed (1 segment)
fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
    feed_tag, seg_feed, ...
    0, 0, z_disk_junction, ...
    0, 0, z_cone_junction, ...
    a);

fprintf(fid,'GE 0\n');

% FR sweep: FR 0 N 0 0 f_start step
fprintf(fid,'FR 0 %d 0 0 %.3f %.3f\n', n_steps, f_start, f_step);

% Excitation on feed wire
fprintf(fid,'EX 0 %d %d 0 1 0\n', feed_tag, feed_seg);

fprintf(fid,'EN\n');
fclose(fid);

disp('OK: Created discone_geometry.nec (500 MHz) and discone_sweep.nec (250-2000 MHz step 10 MHz)');
disp(['Cone ΔL = ', num2str(Lc/seg_c*100), ' cm (target ~1.58 cm)']);
