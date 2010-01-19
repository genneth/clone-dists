function [xis, ts] = xi(lambda, r, gamma, t, z)

% we seek xi((1-z) exp(-g tau)) where:
g = gamma/lambda;
tau = t*lambda;

if tau == 0
    xi = z;
    return;
end

    % g xi' = xi - r/u * (xi+u-1)^2
    % and d(xi)/dt = xi' * (-g u)
    function dxi = dxi(t0, y0)
        xi0 = complex(y0(1), y0(2));
        u = (1-z) * exp(-g*t0);
        xi1 = - u * xi0 + r * (xi0 + u-1)^2;
        dxi = [real(xi1); imag(xi1)];
    end

options = odeset('RelTol',1e-7,'AbsTol',1e-14);
if nargout <= 1
    [ts, xis] = ode45(@dxi, [0 tau/2 tau], [real(z); imag(z)], options);
    xis = xis(end,:);
else
    [ts, xis] = ode45(@dxi, [0 tau], [real(z); imag(z)], options);
end

ts = ts ./ lambda;
xis = complex(xis(:,1), xis(:,2));
xi = xis(end);

end