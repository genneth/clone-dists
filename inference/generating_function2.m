function f = generating_function2(r, gamma, ts, x0, y0)

% solve \dot{f} = [r(x-y)^2 + x(y-1)] f_x + \gamma (1-y) f_y
% by method of characteristics:
%     \dot{x} = x(1-y) - r(x-y)^2
%     \dot{y} = \gamma (y-1)
% starting at t = t0, x = x0, y = y0
% ending at t = 0; initial condition: f(0, x, y) = x

% two modifications:
%  1. t -> -t; the equations are t-translation invariant, and we would like to
%     be able to return the answer for all ts
%  2. use exact solution to y

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

    function yt = y(t)
        yt = exp(-gamma * t) * (y0-1) + 1;
    end

    function dX = deriv(t, X)
        dX = r * (X - y(t))^2 + X * (y(t) - 1);
    end

X0 = x0;
[~,X] = ode45(@deriv, times, X0, ...
    odeset('RelTol', 1e-10, 'AbsTol', 1e-7));

if numel(ts) == 1
    f = X(end,1);
elseif ts(1) > 0
    f = X(2:end,1);
else
    f = X(:,1);
end

end
