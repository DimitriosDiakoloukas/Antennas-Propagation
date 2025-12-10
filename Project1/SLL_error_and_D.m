function f = SLL_error_and_D(p, dBTargetlv)
    % basic array params
    Ne  = 10;
    lam = 1.0;
    s   = lam/2;
    k0  = 2*pi/lam;
    phi = 0;

    % symmetric currents
    Ivec = [1.0, p(1), p(2), p(3), p(4), p(4), p(3), p(2), p(1), 1.0];

    % objective 1 --> SLL error

    SLL_t = 10^(dBTargetlv/20);

    th   = linspace(0,90,91);
    psi  = k0*s*cosd(th) + phi;
    AF   = zeros(1,numel(th));

    for n = 1:Ne
        AF = AF + Ivec(n)*exp(1j*(n-1)*psi);
    end

    AFmag = abs(AF)/max(abs(AF));
    [pk,~] = findpeaks(AFmag);

    [pk_max, pk_idx] = max(pk);
    if pk_max > 0.98
        pk(pk_idx) = [];
    end

    if isempty(pk)
        err_val = 1e6;
    else
        err_val = mean((pk - SLL_t).^2);
    end

    % objective 2 --> directivity
    den = 0;
    S_I = sum(Ivec);
    num = k0*s*(S_I^2);
    
    for n = 0:Ne-1
        for m = 0:Ne-1
            Iterm = Ivec(n+1)*Ivec(m+1)*exp(1j*(n-m)*phi);
            if n == m
                s_term = k0*s;
            else
                dm     = n-m;
                dm_kd  = dm*k0*s;
                s_term = sin(dm_kd)/dm;
            end
            den = den + Iterm*s_term;
        end
    end

    D = num/real(den);
    f(1) = err_val;
    f(2) = -D;
end
