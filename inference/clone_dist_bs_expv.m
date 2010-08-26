function pbs = clone_dist_bs_expv(Tl, Tr, Tg, P0, r, gamma, ts)

% ts is measured in mean division times
% gamma is the dimensionless stratification rate

% dP/dt = T P ==> P = P0 * exp(T t)

p = path;
path(strcat(pwd, '/clone_dist_bs_expv_private'), p);

P = zeros(numel(P0), numel(ts));
t = 0;
k = untetra(length(P0)) - 1;
pbs = zeros(k+1,k+1,numel(ts));
for i = 1:numel(ts)
    P(:,i) = expv(ts(i) - t, (Tl + r * Tr) + gamma*Tg, P0, 1e-10);
    pbs(:,:,i) = Pbs(P(:,i));
    t = ts(i); P0 = P(:,i);
end

path(p);

end
