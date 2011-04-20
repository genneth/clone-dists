function pb = clone_dist_dilution(r1, r2, gamma, ts, n, N0)

pb = inverse_z_transform(@(z) generating_function_dilution(r1, r2, gamma, ts, repmat(z, 1, n+1), z), N0, 1e-4, 1e-5);

for i = 1:numel(ts)
    pb(:,i) = pb(:,i) / (1 - pb(0+1,i));
    pb(0+1,i) = 0;
end

end
