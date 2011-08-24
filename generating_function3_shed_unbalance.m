function f = generating_function3_shed_unbalance(r, delta, gamma, mu, ts, x0, y0, z0)

% backward Kolmogrov:
%   \dot{x} = (r+delta)x^2 + (1-2r)xy + (r-delta)y^2 - x
%   \dot{y} = \gamma (z-y)
%   \dot{z} = \mu (1-z)
% starting at t = 0, x = x0, y = y0, z = z0
% return x at end

assert(ts(1) >= 0, 'must start at a positive time');
assert(~(numel(ts) == 2 && ts(1) == 0), 'not allowed to give ts=[0 t]');

[rows,~] = size(ts);
if rows ~= 1
    ts = ts';
end

if numel(ts) == 1
    times = [0 ts/2 ts];
elseif ts(1) > 0
    times = [0 ts];
else
    times = ts;
end

[times,~,ind] = unique(times);

    function y = y(t)
        y = ( (z0-1)*gamma * exp(-mu*t) ...
           + ((y0-z0)*gamma + (1-y0)*mu) * exp(-gamma*t)) ...
          / (gamma - mu) ...
          + 1;
    end

    function dX = deriv(t, X)
        dX = (r+delta) * X^2 + (1-2*r)*X*y(t) + (r-delta) * y(t)^2 - X;
    end

X0 = x0;
[~,X] = ode45(@deriv, times, X0, ...
    odeset('RelTol', 3e-14, 'AbsTol', 1e-14));

X = X(ind);

if numel(ts) == 1
    f = X(end);
elseif ts(1) > 0
    f = X(2:end);
else
    f = X;
end

end
