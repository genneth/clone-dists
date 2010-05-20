function scaling_plot(gamma, q, ts)

fh = figure;
gh = newplot(fh);
set(gh, 'NextPlot', 'add');

for i = 1:numel(ts)
    av = average_size(gamma, q, ts(i));
    p = clone_dist(gamma, q, ts(i), ceil(av) * 5 * 4);
    semilogy(gh, (1:(ceil(av)*5)) / av, p(1:(ceil(av)*5)) * av);
    pause(0.1);
end

xlabel(gh, 'scaled size', 'FontName', 'Times', 'FontSize', 8);

set(gh, 'FontName', 'Times', 'FontSize', 8);

print(fh, '-dpng', '-painters', '-r300', 'scaling_plot');

end