clc; clear; close all;

% ====== Δ-Λ case ======
lambda = 0.6;                 % m -> f0 = 500 MHz
Rd     = 0.34*lambda;         % 0.204
Lc     = 0.5*lambda;          % 0.3
theta  = 30*pi/180;           % half-angle (2θ0=60°)
gap    = lambda/20;           % 0.03
d      = lambda/200;          % 0.003
a      = d/2;                 % NEC radius

% Cone geometry
h     = Lc*cos(theta);
Rcone = Lc*sin(theta);

% Segmentation (όπως τεκμηριώνεται)
seg_d = 14;
seg_c = 19;

% Junctions
z_disk = 0;
z_cone = -gap;
z_skirt = z_cone - h;

% Feed wire (1 segment) μεταξύ των junctions
feed_tag = 17;
feed_seg = 1;

freqs = [500 1000 1500 2000]; % MHz

for fMHz = freqs
    fname = sprintf('discone_ff_%d.nec', fMHz);
    fid = fopen(fname,'w');

    fprintf(fid,'CE Discone FarField @ %d MHz (same geometry, single-frequency)\n', fMHz);

    tag = 1;

    % ---- Disk: 8 radials ----
    for i = 0:7
        phi = i*2*pi/8;
        x2  = Rd*cos(phi);
        y2  = Rd*sin(phi);

        fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
            tag, seg_d, ...
            0,0,z_disk, ...
            x2,y2,z_disk, ...
            a);
        tag = tag + 1;
    end

    % ---- Cone: 8 wires ----
    for i = 0:7
        phi = i*2*pi/8;
        x2  = Rcone*cos(phi);
        y2  = Rcone*sin(phi);

        fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
            tag, seg_c, ...
            0,0,z_cone, ...
            x2,y2,z_skirt, ...
            a);
        tag = tag + 1;
    end

    % ---- Feed wire: 1 segment ----
    fprintf(fid,'GW %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        feed_tag, 1, ...
        0,0,z_disk, ...
        0,0,z_cone, ...
        a);

    fprintf(fid,'GE 0\n');

    % Single-frequency
    fprintf(fid,'FR 0 1 0 0 %.3f 0\n', fMHz);

    % Excitation on feed wire (segment 1)
    fprintf(fid,'EX 0 %d %d 0 1 0\n', feed_tag, feed_seg);

    fprintf(fid,'EN\n');
    fclose(fid);

    fprintf('Wrote %s\n', fname);
end
