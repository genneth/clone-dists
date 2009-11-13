function P = population(Tl, Tr, Tg, lambda, r, gamma, t, P0)

% dP/dt = T P ==> P = P0 * exp(T t)
sz = size(P0);
P1 = reshape(P0, prod(sz), 1);
P = expv(t, lambda*(Tl + r * Tr) + gamma*Tg, P1, 1e-10);
if sum(P) < 0.99 
    warning('losing probability: %f; lambda t = %f, gamma t = %f', 1 - sum(P), lambda*t, gamma*t);
end
P = max(P, 0); % no negative probabilities
P = reshape(P, sz);

end
