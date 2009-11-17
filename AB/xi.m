function xi = xi(lambda, r, gamma, t, z)

% we seek xi((1-z) exp(-g tau)) where:
g = gamma/lambda;
tau = t*lambda;
    % g xi' = xi - r/u * (xi+u-1)^2
    % and d(xi)/dt = xi' * (-g u)
    function dxi = dxi(t0, y0)
        xi0 = complex(y0(1), y0(2));
        u = (1-z) * exp(-g*t0);
        xi1 = - u * xi0 + r * (xi0 + u-1)^2;
        dxi = [real(xi1); imag(xi1)];
    end

options = odeset('RelTol',1e-7,'AbsTol',1e-14);
[ts, xis] = ode45(@dxi, [0 tau], [real(z); imag(z)], options);

%plot(ts, xis(:,1), ts, xis(:,2));

xi = complex(xis(end,1), xis(end,2));

end