function f = generating_function(r1, r2, ts, x0, y0)

% solve \dot{f} = (r1 x^2 + (1-r1-r2) x y + r2 y^2 - x) f_x
% by method of characteristics:
%     \dot{x} = -(r1 x^2 + (1-r1-r2) x y + r2 y^2 - x)
%     \dot{y} = 0
% starting at t = t0, x = x0, y = y0
% ending at t = 0; initial condition: f(0, x, y) = x

% modification:
%  1. t -> -t; the equations are t-translation invariant, and we would like to
%     be able to return the answer for all ts

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

y=y0;
    
    function dX = deriv(~, X)
        dX = r1 * X^2 + (1-r1-r2) * X * y + r2 * y^2 - X;
    end

[~,X] = ode45(@deriv, times, x0, ...
    odeset('RelTol', 1e-10, 'AbsTol', 1e-7));

if numel(ts) == 1
    f = X(end,1);
elseif ts(1) > 0
    f = X(2:end,1);
else
    f = X(:,1);
end

end
