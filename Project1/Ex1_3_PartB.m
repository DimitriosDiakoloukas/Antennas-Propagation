clc;
clear;
close all;

% Basic parameters
N  = 10;
lam = 1;
d   = lam/2;
k   = 2*pi/lam;
phi = 0;

% --------------------  SLL=-20 dB  --------------------
p20 = [3.135379927369486; 4.665416378742279; 9.686730524803112; 9.060534176008748];
%weights from second optimization
%p20 = [5.7998; 5.4897; 7.4711; 8.1992];
D20      = dir_full(p20);
D20_dBi  = 10*log10(D20);

fprintf("\n--- SLL=-20 dB ---\n");
fprintf("Non-uniform:  %.3f dBi\n", D20_dBi);

pU = [1 1 1 1];
D20_uni     = dir_full(pU);
D20_uni_dBi = 10*log10(D20_uni);
fprintf("Uniform (1.2): %.3f dBi\n", D20_uni_dBi);

D_uni_cf     = dir_uniform(k, N, d, phi);
D_uni_cf_dBi = 10*log10(D_uni_cf);
fprintf("Uniform (1.1): %.3f dBi\n", D_uni_cf_dBi);


% --------------------  SLL=-30 dB  --------------------
p30 = [2.195133541389668; 3.543791802708922; 7.060181489069492; 9.736232397045448];
%weights from second optimization
%p30 = [3.9001; 5.7543; 8.4134; 9.2211];
D30      = dir_full(p30);
D30_dBi  = 10*log10(D30);

fprintf("\n--- SLL=-30 dB ---\n");
fprintf("Non-uniform:  %.3f dBi\n", D30_dBi);

D30_uni     = dir_full(pU);
D30_uni_dBi = 10*log10(D30_uni);
fprintf("Uniform (1.2): %.3f dBi\n", D30_uni_dBi);

D_uni_cf_dBi = 10*log10(dir_uniform(k, N, d, phi));
fprintf("Uniform (1.1): %.3f dBi\n", D_uni_cf_dBi);


% --------------------  SLL=-40 dB  --------------------
p40 = [2.016247175487722; 3.245847814847494; 5.818461419707745; 8.682570442685023];
%weights from second optimization
%p40 = [2.4398; 5.1124; 7.9882; 9.9542];
D40      = dir_full(p40);
D40_dBi  = 10*log10(D40);

fprintf("\n--- SLL=-40 dB ---\n");
fprintf("Non-uniform:  %.3f dBi\n", D40_dBi);

D40_uni     = dir_full(pU);
D40_uni_dBi = 10*log10(D40_uni);
fprintf("Uniform (1.2): %.3f dBi\n", D40_uni_dBi);

D_uni_cf_dBi = 10*log10(dir_uniform(k, N, d, phi));
fprintf("Uniform (1.1): %.3f dBi\n", D_uni_cf_dBi);


% full directivity for non-uniform 4-term symmetric taper
function D = dir_full(p)

    N  = 10;
    lam = 1;
    d   = lam/2;
    k   = 2*pi/lam;
    phi = 0;

    I0 = 1;
    I  = [I0 p(1) p(2) p(3) p(4) p(4) p(3) p(2) p(1) I0];

    S = sum(I);
    num = k*d*(S^2);

    den = 0;
    for n = 0:N-1
        for m = 0:N-1
            Inm = I(n+1)*I(m+1)*exp(1j*(n-m)*phi);
            if n==m
                s = k*d;
            else
                s = sin((n-m)*k*d)/(n-m);
            end
            den = den + Inm*s;
        end
    end

    D = num/real(den);
end


% closed-form directivity (uniform array)
function D = dir_uniform(k, N, d, phi)

    Si = @(x) integral(@(t) sin(t)./t, 0, x);

    a = N*(-k*d + phi)/2;
    b = N*( k*d + phi)/2;

    if abs(b) < 1e-5
        D = (N*k*d)/((sin(a)^2)/a - Si(2*a));
    else
        D = (N*k*d)/((sin(a)^2)/a - (sin(b)^2)/b + Si(2*b) - Si(2*a));
    end
end
