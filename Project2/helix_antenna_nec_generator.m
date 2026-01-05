clc; clear; close all;

% ===== Constants =====
c0 = 3e8;
f0 = 400e6;
lambda = c0/f0;

% ===== Helix parameters =====
N  = 10;
C  = lambda;
D  = C/pi;
R  = D/2;
S  = lambda/4;

% ===== Wire =====
d = lambda/100;
a = d/2;

% ===== Feed & ground =====
gap = lambda/40;
Rg  = lambda/2;

% ===== Spider-web ground =====
Nrad = 8;
Ncir = 4;

phis  = (0:Nrad-1) * (2*pi/Nrad);
ringR = (1:Ncir) * (Rg/Ncir);

% ===== Helix segmentation =====
seg_per_turn = 17;
seg_helix = seg_per_turn * N;
t = linspace(0, 2*pi*N, seg_helix+1);

hx = R*cos(t);
hy = R*sin(t);
hz = gap + (S/(2*pi))*t;

hx(1) = R; hy(1) = 0; hz(1) = gap;

% ===== Feed =====
feed_tag = 100;
feed_seg = 1;
seg_feed = 1;

% ===== Ground nodes =====
nodeX = zeros(Nrad, Ncir);
nodeY = zeros(Nrad, Ncir);

for i = 1:Nrad
    for k = 1:Ncir
        nodeX(i,k) = ringR(k)*cos(phis(i));
        nodeY(i,k) = ringR(k)*sin(phis(i));
    end
end

feed_gnd_x = R;
feed_gnd_y = 0;

% ================= FILE 1 =================
fid = fopen('helix_geometry.nec','w');
fprintf(fid,'CE Helix + spider-web ground\n');

tag = 1;

for k = 1:seg_helix
    if tag == feed_tag, tag = tag + 1; end
    fprintf(fid,'GW %d 1 %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        tag, hx(k), hy(k), hz(k), hx(k+1), hy(k+1), hz(k+1), a);
    tag = tag + 1;
end

fprintf(fid,'GW %d %d %.6f %.6f 0 %.6f %.6f %.6f %.6f\n', ...
    feed_tag, seg_feed, feed_gnd_x, feed_gnd_y, feed_gnd_x, feed_gnd_y, gap, a);

for i = 1:Nrad
    rr = [0, ringR];
    if i == 1 && ~any(abs(rr - R) < 1e-12)
        rr = sort([rr, R]);
    end

    for k = 1:numel(rr)-1
        r1 = rr(k); r2 = rr(k+1);

        if r1 == 0
            x1 = 0; y1 = 0;
        elseif any(abs(r1-ringR)<1e-12)
            idx = find(abs(r1-ringR)<1e-12,1);
            x1 = nodeX(i,idx); y1 = nodeY(i,idx);
        else
            x1 = feed_gnd_x; y1 = feed_gnd_y;
        end

        if any(abs(r2-ringR)<1e-12)
            idx = find(abs(r2-ringR)<1e-12,1);
            x2 = nodeX(i,idx); y2 = nodeY(i,idx);
        else
            x2 = feed_gnd_x; y2 = feed_gnd_y;
        end

        fprintf(fid,'GW %d 1 %.6f %.6f 0 %.6f %.6f 0 %.6f\n', tag, x1, y1, x2, y2, a);
        tag = tag + 1;
    end
end

for k = 1:Ncir
    for i = 1:Nrad
        ip = mod(i,Nrad)+1;
        fprintf(fid,'GW %d 1 %.6f %.6f 0 %.6f %.6f 0 %.6f\n', ...
            tag, nodeX(i,k), nodeY(i,k), nodeX(ip,k), nodeY(ip,k), a);
        tag = tag + 1;
    end
end

fprintf(fid,'GE 0\n');
fprintf(fid,'FR 0 1 0 0 400 0\n');
fprintf(fid,'EX 0 %d %d 0 1 0\n', feed_tag, feed_seg);
fprintf(fid,'EN\n');
fclose(fid);

disp('OK: helix_geometry.nec created');

% ================= FILE 2 =================
f_start = 120;
f_stop  = 800;
f_step  = 5;
n_steps = round((f_stop - f_start)/f_step) + 1;

fid = fopen('helix_sweep.nec','w');
fprintf(fid,'CE Helix + spider-web ground sweep\n');

tag = 1;

for k = 1:seg_helix
    if tag == feed_tag, tag = tag + 1; end
    fprintf(fid,'GW %d 1 %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', ...
        tag, hx(k), hy(k), hz(k), hx(k+1), hy(k+1), hz(k+1), a);
    tag = tag + 1;
end

fprintf(fid,'GW %d %d %.6f %.6f 0 %.6f %.6f %.6f %.6f\n', ...
    feed_tag, seg_feed, feed_gnd_x, feed_gnd_y, feed_gnd_x, feed_gnd_y, gap, a);

for i = 1:Nrad
    rr = [0, ringR];
    if i == 1 && ~any(abs(rr-R)<1e-12)
        rr = sort([rr,R]);
    end

    for k = 1:numel(rr)-1
        r1 = rr(k); r2 = rr(k+1);

        if r1 == 0
            x1 = 0; y1 = 0;
        elseif any(abs(r1-ringR)<1e-12)
            idx = find(abs(r1-ringR)<1e-12,1);
            x1 = nodeX(i,idx); y1 = nodeY(i,idx);
        else
            x1 = feed_gnd_x; y1 = feed_gnd_y;
        end

        if any(abs(r2-ringR)<1e-12)
            idx = find(abs(r2-ringR)<1e-12,1);
            x2 = nodeX(i,idx); y2 = nodeY(i,idx);
        else
            x2 = feed_gnd_x; y2 = feed_gnd_y;
        end

        fprintf(fid,'GW %d 1 %.6f %.6f 0 %.6f %.6f 0 %.6f\n', tag, x1, y1, x2, y2, a);
        tag = tag + 1;
    end
end

for k = 1:Ncir
    for i = 1:Nrad
        ip = mod(i,Nrad)+1;
        fprintf(fid,'GW %d 1 %.6f %.6f 0 %.6f %.6f 0 %.6f\n', ...
            tag, nodeX(i,k), nodeY(i,k), nodeX(ip,k), nodeY(ip,k), a);
        tag = tag + 1;
    end
end

fprintf(fid,'GE 0\n');
fprintf(fid,'FR 0 %d 0 0 %.3f %.3f\n', n_steps, f_start, f_step);
fprintf(fid,'EX 0 %d %d 0 1 0\n', feed_tag, feed_seg);
fprintf(fid,'EN\n');
fclose(fid);

fprintf('Sweep: %d–%d MHz, step %d MHz (%d points)\n', f_start, f_stop, f_step, n_steps);
