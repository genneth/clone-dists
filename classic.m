function classic()

N = 49; % max number of cells
lambda = 1.1;
r = 0.08;
rho = 0.22;
%gamma = lambda*rho/(1-rho);
gamma = 0.0;

pmn = zeros(N+1,N+1);
pmn(1+1,0+1) = 1; % initial conditions

pi_s = zeros(2*N, 1); % clone sizes (non-extinct)
pm = zeros(1,N+1); % P_m

% some utilities to convert to and from a 2D matrix
sz = size(pmn);
    function i = s2i(m,n)
        i = sub2ind(sz, m+1, n+1);
    end
    function [m,n] = i2s(i)
        [m1,n1] = ind2sub(sz, i)
        m = m1-1; n = n1-1;
    end

T = sparse((N+1)^2,(N+1)^2);

% form the transition matrix
for m=0:N
    for n=0:N
        if (m-1) >= 0
            T(s2i(m,n), s2i(m-1,n)) = lambda*r*(m-1);
        end
        if (m+1) <= N && (n-2) >= 0
            T(s2i(m,n), s2i(m+1,n-2)) = lambda*r*(m+1);
        end
        if (n-1) >= 0
            T(s2i(m,n), s2i(m,n-1)) = lambda*(1-2*r)*m;
        end
        if (n+1) <= N
            T(s2i(m,n), s2i(m,n+1)) = gamma*(n+1);
        end
        T(s2i(m,n),s2i(m,n)) = -lambda*m - gamma*n;
    end
end

%surf(pmn); pause(0.01);
while 1
    pmn1 = reshape(pmn, (N+1)^2, 1);
    pmn1 = expv(1.0, T, pmn1);
    pmn = reshape(pmn1, sz);
    pm = sum(pmn, 2);
    for s = [1:2*N]
        pi_s(s) = 0;
        for m = [0:N]
            n = s - m;
            if n >= 0 && n <= N
                pi_s(s) = pi_s(s) + pmn(m+1,n+1);
            end
        end
    end
    %surf(pmn);
    stairs([1:2*N], pi_s);
    pause(0.1);
end

end