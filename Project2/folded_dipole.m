clc; clear; close all;

% ========= ΕΠΙΛΟΓΗ ΟΜΑΔΑΣ =========
% Δ-Λ -> f0 = 60 MHz (από εκφώνηση Θέμα 3α)
f0 = 60e6;

c0 = 3e8;
lambda = c0/f0;

% ========= ΓΕΩΜΕΤΡΙΑ =========
L = lambda/2;                 % folded dipole συνολικό μήκος
d = lambda/200;               % διάμετρος σύρματος (ζητούμενο)
a = d/2;                      % NEC θέλει ακτίνα

% 3 περιπτώσεις απόστασης (s)
s_list = [lambda/100, lambda/20, lambda/4];
name_list = {'s_lam100', 's_lam20', 's_lam4'};

% ========= SWEEP 0.5f0 έως 1.5f0 =========
% Για f0=60 MHz => 30..90 MHz (όπως φαίνεται και στο report σου). :contentReference[oaicite:1]{index=1}
f_start_MHz = 0.5*(f0/1e6);
f_stop_MHz  = 1.5*(f0/1e6);
f_step_MHz  = 0.5;                % μπορείς 0.5 ή 1 MHz. (όσο μικρότερο τόσο πιο "ομαλό")
n_steps = round((f_stop_MHz - f_start_MHz)/f_step_MHz) + 1;

% ========= ΚΑΤΑΤΜΗΣΗ =========
% Θέλουμε ΟΠΩΣΔΗΠΟΤΕ περιττό πλήθος segments ώστε να υπάρχει "κεντρικό" segment
seg_long = 41;                    % καλό default (περιττό)
feed_seg = (seg_long+1)/2;        % κεντρικό segment

% ========= ΠΑΡΑΓΩΓΗ 3 NEC =========
for idx = 1:numel(s_list)
    s = s_list(idx);
    out_nec = sprintf('folded_dipole_%s.nec', name_list{idx});

    fid = fopen(out_nec,'w');
    fprintf(fid,'CE Folded dipole (lambda/2), d=lambda/200, spacing s=%.6f m\n', s);
    fprintf(fid,'CE Sweep %.1f..%.1f MHz step %.1f MHz\n', f_start_MHz, f_stop_MHz, f_step_MHz);

    % Συντεταγμένες:
    % δύο παράλληλα σύρματα στον άξονα z, σε x = ±s/2
    x1 = -s/2;  x2 = +s/2;
    zA = -L/2;  zB = +L/2;

    tag1 = 1;   % driven leg
    tag2 = 2;   % folded leg
    tag3 = 3;   % top short
    tag4 = 4;   % bottom short

    % 1) driven leg (x=-s/2)
    fprintf(fid,'GW %d %d %.6f 0 %.6f  %.6f 0 %.6f  %.6f\n', ...
        tag1, seg_long, x1, zA,  x1, zB, a);

    % 2) folded leg (x=+s/2)
    fprintf(fid,'GW %d %d %.6f 0 %.6f  %.6f 0 %.6f  %.6f\n', ...
        tag2, seg_long, x2, zA,  x2, zB, a);

    % 3) short at top (z=+L/2) connecting the two legs
    fprintf(fid,'GW %d 1 %.6f 0 %.6f  %.6f 0 %.6f  %.6f\n', ...
        tag3, x1, zB,  x2, zB, a);

    % 4) short at bottom (z=-L/2)
    fprintf(fid,'GW %d 1 %.6f 0 %.6f  %.6f 0 %.6f  %.6f\n', ...
        tag4, x1, zA,  x2, zA, a);

    % End geometry
    fprintf(fid,'GE 0\n');

    % Frequency sweep
    fprintf(fid,'FR 0 %d 0 0 %.3f %.3f\n', n_steps, f_start_MHz, f_step_MHz);

    % Excitation: voltage source στο κέντρο του driven leg (tag1, feed_seg)
    fprintf(fid,'EX 0 %d %d 0 1 0\n', tag1, feed_seg);

    fprintf(fid,'EN\n');
    fclose(fid);

    fprintf('OK: wrote %s  (s=%.4f m)\n', out_nec, s);
end

fprintf('\nΆνοιξε κάθε .nec στο 4nec2 και πάτα Calculate -> (F5) Impedance για Zin.\n');
fprintf('Για Reflection/SWR άλλαξε το reference Z0 στο plot (π.χ. 50Ω/200Ω/300Ω) και διάλεξε αυτό που ελαχιστοποιεί |Γ| κοντά σε Xin=0.\n');
