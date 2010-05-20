function f = generating_function(gamma, q, t0, x0, y0)

% solve \dot{f} = x (y-1) f_x + \gamma [qx + (1-q) - y] f_y
% by method of characteristics:
%     \dot{x} = x (1-y)
%     \dot{y} = \gamma [y - qx + (q-1)]
% starting at t = t0, x = x0, y = y0
% ending at t = 0; initial condition: f(0, x, y) = x

    function dY = deriv(t, Y)
        x = Y(1); y = Y(2);
        dx = x * (1-y);
        dy = gamma * (y - q*x - (1-q));
        dY = [dx; dy];
    end

X0 = [x0; y0];
[T,X] = ode45(@deriv, [t0 t0/2 0], X0, ...
    odeset('RelTol', 1e-10, 'AbsTol', 1e-7));

f = X(end, 1);

end