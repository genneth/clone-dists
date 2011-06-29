function average_size_exp()

% just plot the experiment, no fits

oesophagus_data;

% condition the data and find the maximum and averages
for i=1:numel(ts2)
    basal{i}(1+1) = 0;
    tot(i) = sum(basal{i});
    av(i) = dot(0:(numel(basal{i})-1), basal{i}/tot(i));
end

fh = figure;
set(fh, 'PaperUnits', 'inches');
w = 7; h = 4.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

ah = newplot(fh);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 10);
plot(ah, ts2, av, 's');
set(ah, 'XScale', 'log', 'YScale', 'log');
xlabel(ah, 'time (weeks)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
ylabel(ah, 'average basal clone size', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');

print(fh, '-depsc2', '-painters', 'basal-average-size-exp-log');

set(ah, 'XScale', 'linear', 'YScale', 'linear');
print(fh, '-depsc2', '-painters', 'basal-average-size-exp-linear');


close(fh);


end