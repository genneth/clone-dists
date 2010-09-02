function plot_posterior_dist(ah, samples)

coords = [[samples{:,2}]' [samples{:,3}]' [samples{:,4}]'];
coords(:,2) = coords(:,2) ./ (coords(:,2) + 1); % rho instead of gamma
lps = cellfun(@(p2s, p3s)(sum(p2s) + sum(p3s)), {samples{:,5}}, {samples{:,6}});
ps = exp(lps - max(lps)); % normalised to max(ps) == 1

contour3_tetra(ah, coords(:,1), coords(:,2), coords(:,3), ps, 0.5);

xlabel(ah, 'r');
ylabel(ah, '\rho');
zlabel(ah, '\lambda');

end
