function decaying()

N = 49; % max number of cells
G = 39; % max number of generations
lambda = 1.1;
r = 0.08;

pmg = zeros(N+1,G+1);
pmg(1+1,0+1) = 1; % initial conditions

% some utilities to convert to and from a 2D matrix
sz = size(pmg);
    function i = s2i(m,g)
        i = sub2ind(sz, m+1, g+1);
    end
    function [m,n] = i2s(i)
        [m1,g1] = ind2sub(sz, i)
        m = m1-1; n = g1-1;
    end

T = sparse((N+1)*(G+1),(N+1)*(G+1));

% form the transition matrix
for m=0:N
    for g=0:G
        if (m-1) >= 0 && (g-1) >= 0
            T(s2i(m,g), s2i(m-1,g-1)) = lambda*r*(m-1);
        end
        if (m+1) <= N && (g-1) >= 0
            T(s2i(m,g), s2i(m+1,g-1)) = lambda*r*(m+1);
        end
        if (g-1) >= 0
            T(s2i(m,g), s2i(m,g-1)) = lambda*(1-2*r)*m;
        end
        T(s2i(m,g),s2i(m,g)) = -lambda*m;
    end
end

surf(pmg); pause(0.01);
while 1
    pmg1 = reshape(pmg, (N+1)*(G+1), 1);
    pmg1 = expv(1.0, T, pmg1);
    pmg = reshape(pmg1, sz);
    surf(pmg);
    pause(0.1);
end

end