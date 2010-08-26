function result = monte_carlo_integrate3(integrand, samples)

% we assume that samples is a cell array, with structure:
% samples{i,:} = {prob, r, gamma, lambda, [log p2], [log p3]}
% integrand :: (r, gamma, lambda) -> Double

lps = arrayfun(@(p2s, p3s)(sum(p2s) + sum(p3s)), [samples{:,5}], [samples{:,6}]);
ps = exp(lps - max(lps)); % normalised to max(ps) == 1

result = 0;
for i = 1:numel(ps)
    result = result + ...
        integrand(samples{i,2}, samples{i,3}, samples{i,4}) * ...
        ps(i) / samples{i,1};
end

end
