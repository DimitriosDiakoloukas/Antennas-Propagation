clc; clear; close all;

% Parameters Init
f0 = 12;
lambda = 300 / f0;
L_total = 5 * lambda;
h = lambda / 4;

a_wire = (lambda / 200) / 2;
a_feed = a_wire / 4;

segs_hor = 101;
segs_vert = 5;

RL = 317.54;

% Sweep: 0.5f0 .. 1.5f0
f_start = 0.5 * f0;     % MHz
f_stop  = 1.5 * f0;     % MHz
f_step  = 0.1;          % MHz
n_steps = round((f_stop - f_start)/f_step) + 1;

% Two cases: perfect & good
cases = {
    struct('name','perfect','GN','GN  1  0  0  0  0  0  0  0')
    struct('name','good',   'GN','GN  2  0  0  0  17  0.015  0  0')
};

for k = 1:numel(cases)
    filename = sprintf('traveling_wave_antenna_geometry_%s.nec', cases{k}.name);
    fid = fopen(filename, 'w');

    fprintf(fid, 'CM Traveling Wave Horizontal Antenna (Beverage Type)\n');
    fprintf(fid, 'CM Source at Top (junction), Load at Bottom (ground connection)\n');
    fprintf(fid, 'CM Sweep %.2f..%.2f MHz step %.2f MHz\n', f_start, f_stop, f_step);
    fprintf(fid, 'CE\n');

    % Geometry
    fprintf(fid, 'GW  1  %d  %.4f  0.0000  %.4f  %.4f  0.0000  %.4f  %.4f\n', ...
        segs_hor, -L_total/2, h, L_total/2, h, a_wire);

    fprintf(fid, 'GW  2  %d  %.4f  0.0000  %.4f  %.4f  0.0000  0.0000  %.4f\n', ...
        segs_vert, L_total/2, h, L_total/2, a_feed);

    fprintf(fid, 'GW  3  %d  %.4f  0.0000  %.4f  %.4f  0.0000  0.0000  %.4f\n', ...
        segs_vert, -L_total/2, h, -L_total/2, a_feed);

    fprintf(fid, 'GE  1\n');

    % Ground
    fprintf(fid, '%s\n', cases{k}.GN);

    % Load
    fprintf(fid, 'LD  0  2  5  5  %.1f  0  0\n', RL);

    % Excitation
    fprintf(fid, 'EX  0  3  1  0  1.0  0.0\n');

    % Frequency sweep
    fprintf(fid, 'FR  0  %d  0  0  %.2f  %.2f\n', n_steps, f_start, f_step);

    fprintf(fid, 'EN\n');
    fclose(fid);

    disp(['File ', filename, ' is created']);
end
