function f = generating_function(r, gamma, mu, ts, x0, y0, z0)

% solve \dot{f} = [r(x-y)^2 + x(y-1)] f_x + \gamma (z-y) f_y + \mu (1-z) f_z
% by method of characteristics:
%     \dot{x} = x(1-y) - r(x-y)^2
%     \dot{y} = \gamma (y-z)
%     \dot{z} = \mu (z-1)
% starting at t = t0, x = x0, y = y0, z=z0
% ending at t = 0; initial condition: f(0, x, y, z) = x

% two modifications:
%  1. t -> -t; the equations are t-translation invariant, and we would like to
%     be able to return the answer for all ts
%  2. use exact solutions to y and z

assert(ts(1) >= 0, 'must start at a positive time');
assert(numel(ts) > 1, 'must give a range of times to compute');

if ts(1) > 0
    times = [0 ts];
else
    times = ts;
end

    function z = z(t)
        z = exp(-mu * t) * (z0-1) + 1;
    end

    function y = y(t)
        y = ( (z0-1)*gamma * exp(-mu*t) ...
           + ((y0-z0)*gamma + (1-y0)*mu) * exp(-gamma*t)) ...
          / (gamma - mu) ...
          + 1;
    end

    function dX = deriv(t, X)
        dX = r * (X - y(t))^2 + X * (y(t) - 1);
    end

X0 = x0;
[T,X] = ode45(@deriv, times, X0, ...
    odeset('RelTol', 1e-10, 'AbsTol', 1e-7));

if ts(1) > 0
    f = X(2:end);
else
    f = X;
end

end
