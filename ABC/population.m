function P = population(Tl, Tr, Tg, lambda, r, gamma, t, P0)

% dP/dt = T P ==> P = P0 * exp(T t)

% P = expv(t, lambda*(Tl + r * Tr) + gamma*Tg, P0, 1e-10);

% dP/dt = T P + u 

u = zeros(size(P0));
u(1) = -lambda * 0.01;
u(2) = lambda * 0.01;
P = phiv(t, lambda*(Tl + r * Tr) + gamma*Tg, u, P0, 1e-10);

end
