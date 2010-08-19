function condpbs = clone_size_distribution_cond(r, gamma, mu, ts, M0, N0)

pbs = clone_size_distribution(r, gamma, mu, ts, M0, N0);
p0 = prob_zero_basal(r, gamma, ts);

Z1 = 1 - p0' - squeeze(pbs(1+1,0+1,:)); % 1 - Z

condpbs = zeros(size(pbs));
for k = 1:numel(ts)
    condpbs(:,:,k) = pbs(:,:,k) ./ Z1(k);
end
condpbs(0+1,:,:)   = 0;
condpbs(1+1,0+1,:) = 0;

end
