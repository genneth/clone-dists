function [xis, ts] = xi(lambda, r, gamma, t, z)

% we seek xi((1-z) exp(-g tau)) where:
g = gamma/lambda;
tau = t*lambda;
if numel(t) == 1
    taus = [0 tau/2 tau];
else
    if t(1) > 0
        taus = [0 tau];
    else
        taus = tau;
    end
end

if numel(t) == 1 && t == 0
    ts = 0;
    xis = z;
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
if numel(taus) ~= 2
    [ts, xis] = ode45(@dxi, taus, [real(z); imag(z)], options);
else
    [ts, xis] = ode45(@dxi, [taus(1) (taus(1)+taus(2))/2 taus(2)], [real(z); imag(z)], options);
    ts = [ts(1); ts(3)];
    xis = [xis(1,:); xis(3,:)];
end
if numel(t) == 1 && nargout <= 1
    xis = xis(end,:); % ts does not matter
else
    if t(1) > 0
        ts = ts(2:end);
        xis = xis(2:end,:);
    end
end

ts = ts ./ lambda;
xis = complex(xis(:,1), xis(:,2));

end
