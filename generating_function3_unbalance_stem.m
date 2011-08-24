function f = generating_function3_unbalance_stem(r, delta, gamma, ts, x0, y0, z0)

% backward Kolmogrov:
%   \dot{x} = (r+delta)x^2 + (1-2r)xy + (r-delta)y^2 - x
%   \dot{y} = \gamma (1-y)
%   \dot{z} = -2*delta (xz - z)
% starting at t = 0, x = x0, y = y0, z = z0
% return z at end

assert(ts(1) >= 0, 'must start at a positive time');
assert(~(numel(ts) == 2 && ts(1) == 0), 'not allowed to give ts=[0 t]');
assert(delta < 0, 'must have subcritical CP');

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

    function yt = y(t)
        yt = exp(-gamma * t) * (y0-1) + 1;
    end

    function dX = deriv(t, X)
        x = X(1); z = X(2);
        dX = [(r+delta) * x^2 + (1-2*r)*x*y(t) + (r-delta) * y(t)^2 - x;
            -2*delta*(z*x - z)];
    end

X0 = [x0; z0];
[~,X] = ode45(@deriv, times, X0, ...
    odeset('RelTol', 3e-14, 'AbsTol', 1e-14));

X = X(ind,2);
if numel(ts) == 1
    f = X(end);
elseif ts(1) > 0
    f = X(2:end);
else
    f = X;
end

end
