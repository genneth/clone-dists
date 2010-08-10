function scaling_plot(r1, r2, ts, n)

fh = figure;
gh = newplot(fh);
set(gh, 'NextPlot', 'add');
%set(gh, 'YScale', 'log');
set(gh, 'XLim', [0 5]);
drawnow;

ps = clone_dist(r1, r2, ts, n);

for i = 1:numel(ts)
    av = average_size(r1, r2, ts(i));
    p = ps(:,i);
    c = 1 - cumsum(p);
    plot(gh, (1:n) ./ av, c);
    drawnow;
end

xlabel(gh, 'scaled size', 'FontName', 'Times', 'FontSize', 8);
ylabel(gh, 'cumulative distribution', 'FontName', 'Times', 'FontSize', 8);
set(gh, 'FontName', 'Times', 'FontSize', 8);

print(fh, '-dpng', 'scaling_plot');

end
