function average = basal_average(r, gamma, lambda, ts)

theory = clone_dist_b(r, gamma, lambda*ts, 1);
raw_average = 1/gamma * exp(-gamma*lambda*ts) .* expm1(gamma*lambda*ts) + 1;
average = (raw_average - theory(1+1,:)) ./ (1 - theory(0+1,:) - theory(1+1,:));

end
