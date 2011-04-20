function f = generating_function_dilution(r1, r2, gamma, ts, xs0, y0)

% Generic 2-type branching process, with dividing A and differentiated B
% cells. The complication is that we can only track cells for so long
% before the (unreplenished) plasmid/dye is diluted such that we can no
% longer see it. We implement this as a Markovian process with n+1 dividing
% cell types (for a process where we can track n-fold halvings of dye).

assert(ts(1) >= 0, 'must start at a positive time');
assert(~(numel(ts) == 2 && ts(1) == 0), 'not allowed to give ts=[0 t]');

[rows,~] = size(ts);
if rows ~= 1
    ts = ts';
end

% The vector xs0 should be of length n+1
n = numel(xs0) - 1;

if numel(ts) == 1
    times = [0 ts/2 ts];
elseif ts(1) > 0
    times = [0 ts];
else
    times = ts;
end

    function r = g(x1,y,x0)
        r = r1*x1*x1 + (1-r1-r2)*x1*y + r2*y*y - x0;
    end

    function yt = y(t)
        yt = exp(-gamma * t) * (y0-1) + 1;
    end

    function dX = deriv(t, X)
        dX = zeros(size(X));
        for i=0:(n-1)
            dX(i+1) = g(X(i+2),y(t),X(i+1));
        end
        dX(n+1) = g(1,1, X(n+1));
    end

X0 = xs0;
[~,X] = ode45(@deriv, times, X0, ...
    odeset('RelTol', 3e-14, 'AbsTol', 1e-14));

if numel(ts) == 1
    f = X(end,1);
elseif ts(1) > 0
    f = X(2:end,1);
else
    f = X(:,1);
end

end
