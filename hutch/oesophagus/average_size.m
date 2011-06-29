function average_size()

oesophagus_data;

% condition the data and find the maximum and averages
M = 0;
for i=1:numel(ts2)
    basal{i}(1+1) = 0;
    M = max(M, numel(basal{i}));
    tot(i) = sum(basal{i});
    av(i) = dot(0:(numel(basal{i})-1), basal{i}/tot(i));
end
M = M-1;

r = 0.135; % ± 0.001
rho = 0.658; % ± 0.005
lambda = 1.32; % ± 0.00 / week
gamma = rho/(1-rho);

theory = clone_dist_b(r, gamma, lambda*ts2, 1);
raw_average = 1/gamma * exp(-gamma*lambda*ts2) .* expm1(gamma*lambda*ts2) + 1;
average = (raw_average - theory(1+1,:)) ./ (1 - theory(0+1,:) - theory(1+1,:));

fh = figure;
set(fh, 'PaperUnits', 'inches');
w = 7; h = 4.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

ah = newplot(fh);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 10);
plot(ah, ts2, av, 's');
errorbar(ah, ts2, average, average ./ sqrt(tot));
set(ah, 'XScale', 'log', 'YScale', 'log');
xlabel(ah, 'time (weeks)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
ylabel(ah, 'average basal clone size', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');

print(fh, '-depsc2', '-painters', 'basal-average-size-log');

set(ah, 'XScale', 'linear', 'YScale', 'linear');
set(ah, 'XLim', [0 60]);
print(fh, '-depsc2', '-painters', 'basal-average-size-linear');


close(fh);


end
