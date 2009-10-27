function P = population(Tl, Tr, Tg, lambda, r, gamma, t, P0)

% dP/dt = T P ==> P = P0 * exp(T t)
sz = size(P0);
P1 = reshape(P0, prod(sz), 1);
P = expv(t, lambda*(Tl + r * Tr) + gamma*Tg, P1, 1e-10);
if sum(P) < 0.9 
    warning('losing probability: %f; lambda t = %f, gamma t = %f', sum(P), lambda*t, gamma*t);
end
P = reshape(P, sz);

end
