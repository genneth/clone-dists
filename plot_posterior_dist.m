function plot_posterior_dist(ah, samples, pf)

coords = [[samples{:,1}]' [samples{:,2}]' [samples{:,3}]'];
coords(:,2) = coords(:,2) ./ (coords(:,2) + 1); % rho instead of gamma
lps = cellfun(pf, samples(:,4), samples(:,5)) - cellfun(@(l) log(l(1)), samples(:,3));
ps = exp(lps - max(lps)); % normalised to max(ps) == 1

contour3_tetra(ah, coords(:,1), coords(:,2), coords(:,3), ps, 0.3);

xlabel(ah, 'r');
ylabel(ah, '\rho');
zlabel(ah, '\lambda');

end
