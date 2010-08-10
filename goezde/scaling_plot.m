function scaling_plot(r1, r2, ts, n)

fh = figure;
gh = newplot(fh);
set(gh, 'NextPlot', 'add');
%set(gh, 'YScale', 'log');
set(gh, 'XLim', [0 5]);

ps = clone_dist(r1, r2, ts, n);

for i = 1:numel(ts)
    av = average_size(r1, r2, ts(i));
    p = ps(:,i);
    semilogy(gh, (0:(n-1)) ./ av, p .* av);
    drawnow;
end

xlabel(gh, 'scaled size', 'FontName', 'Times', 'FontSize', 8);

set(gh, 'FontName', 'Times', 'FontSize', 8);

print(fh, '-dpng', '-painters', '-r300', 'scaling_plot');

end
