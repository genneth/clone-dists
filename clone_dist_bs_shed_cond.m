function pbs = clone_dist_bs_shed_cond(r, gamma, mu, ts, M0, N0)

% condition on b>0

p0 = generating_function2(r, gamma, ts, 0, 0);
pbs0 = clone_dist_bs_shed(r, gamma, mu, ts, M0, N0);

pbs = pbs0;
for i = 1:numel(ts)
    pbs(0+1,:, i) = 0;
    pbs(:,:, i) = pbs(:,:, i) / (1 - p0(i));
end

end