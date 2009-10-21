function P = population(Tl, Tr, Tg, lambda, r, gamma, t, P0)

% dP/dt = T P ==> P = P0 * exp(T t)

P = expv(t, lambda*(Tl + r * Tr) + gamma*Tg, P0);

end
