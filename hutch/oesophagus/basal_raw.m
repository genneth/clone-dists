function basal_raw

basal = {};
oesophagus_data;
colour_palette;

% pick out david's data and collate
david = ts2 ~= 1 & ts2 ~= 30/7;
ts = unique(ts2(david));
tots = cellfun(@(d) sum(d((2+1):end)), basal);
av = zeros(size(ts));
avdev = zeros(size(ts));
for i = 1:numel(ts)
    ind = find(ts2==ts(i));
    avs = zeros(size(ind));
    for j = 1:numel(ind)
        avs(j) = dot(2:(numel(basal{ind(j)})-1), basal{ind(j)}((2+1):end) / tots(ind(j)));
    end
    av(i) = mean(avs);
    avdev(i) = sqrt(var(avs) * numel(ind) / (numel(ind)-1)); % "sample standard deviation"
%     avdev(i) = std(avs); % "standard deviation of the sample"
end
data = arrayfun(@(t) collate(basal(ts2==t)), ts, 'UniformOutput', false);
m = max(cellfun(@max, data));

fh = figure;
set(fh, 'PaperUnits', 'inches');
w = 4; h = 4;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

gh = axes('OuterPosition', [0 0 1 0.7]);
set(gh, 'NextPlot', 'add');
set(gh, 'XLim', [0.5 numel(ts)+0.5]);
set(gh, 'YScale', 'log', 'YLim', [0.5 200]);
set(gh, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(gh, 'time post-induction', 'FontSize', 9, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(gh, 'YAxisLocation', 'right');
ylabel(gh, 'clone size', 'Rotation', -90, 'FontSize', 9, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(gh, 'Box', 'off');
set(gh, 'YGrid', 'on', 'YMinorGrid', 'off', 'GridLineStyle', '-');
set(gh, 'YTick', 5.^(0:6));
set(gh, 'XTick', 1:numel(ts));
set(gh, 'XTickLabel', {'3d', '10d', '3w', '6w', '3m', '6m', '1y'});

for i=1:numel(ts)
    for j=0:(numel(data{i})-1)
        if data{i}(j+1) ~= 0
            xs = linspace(-1/2,1/2,data{i}(j+1)) * data{i}(j+1)/m + i;
            plot(gh, xs, j*ones(size(xs)), '.', 'Color', colours(mod(i-1,7)+1,:), 'MarkerSize', 4.0);
        end
    end
end
% scale bar
plot(gh, linspace(-1/2,1/2,50)*50/m+4.5, 100*ones(1,50), '.', 'Color', [0 0 0], 'MarkerSize', 4.0);
text(4, 60, '50 clones', 'FontSize', 7, 'FontName', 'Helvetica', 'FontWeight', 'bold');
% remove padding on left
p = get(gh, 'Position');
p = [0.03 p(2) p(3)+p(1)-0.03 p(4)];
set(gh, 'Position', p);

% insets
ah = axes('OuterPosition', [0 0.7 0.45 0.3]);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 8);
% set(ah, 'XScale', 'log', 'YScale', 'log');
set(ah, 'XLim', [0 55]);
set(ah, 'XGrid', 'on', 'XMinorGrid', 'off');
set(ah, 'XTick', [4 12 26 52]);
set(ah, 'XTickLabel', {'1m', '3m', '6m', '1y'});
set(ah, 'YAxisLocation', 'right', 'YGrid', 'on', 'YMinorGrid', 'off');
ylabel(ah, 'average clone size', 'Rotation', -90, 'FontSize', 8, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(ah, 'YLim', [0 25], 'YTick', [0 10 20]);
plot(ah, ts, av, '+', 'MarkerSize', 2.0, 'Color', colours(1,:));
for i = 1:numel(ts)
    plot(ah, [ts(i) ts(i)], av(i) - avdev(i)*[-1 1], '-', 'Color', colours(1,:));
end
theory = [];
load basal-average
plot(ah, linspace(min(ts), max(ts), 100), theory, '-', 'Color', colours(2,:));
p = get(ah, 'Position');
p = [0.03 p(2) p(3)+p(1)-0.03 p(4)];
set(ah, 'Position', p);

ah = axes('OuterPosition', [0.55 0.7 0.45 0.3]);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 8);
set(ah, 'XLim', [0 55]);
set(ah, 'XGrid', 'on', 'XMinorGrid', 'off');
set(ah, 'XTick', [4 12 26 52]);
set(ah, 'XTickLabel', {'1m', '3m', '6m', '1y'});
set(ah, 'YGrid', 'on', 'YMinorGrid', 'off');
ylabel(ah, 'clone density', 'FontSize', 8, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(ah, 'YLim', [0 25], 'YTick', [0 10 20]);
p = get(ah, 'Position');
p = [0.58 p(2) p(3)+p(1)-0.58 p(4)];
set(ah, 'Position', p);

ah = axes('OuterPosition', [0.0 0.4 0.45 0.3]);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 8);
set(ah, 'XLim', [0 55]);
set(ah, 'XGrid', 'on', 'XMinorGrid', 'off');
set(ah, 'XTick', [4 12 26 52]);
set(ah, 'XTickLabel', {'1m', '3m', '6m', '1y'});
set(ah, 'YAxisLocation', 'right', 'YGrid', 'on', 'YMinorGrid', 'off');
ylabel(ah, 'labelling density', 'Rotation', -90, 'FontSize', 8, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(ah, 'YLim', [0 25], 'YTick', [0 10 20]);
p = get(ah, 'Position');
p = [0.03 p(2) p(3)+p(1)-0.03 p(4)];
set(ah, 'Position', p);

print(fh, '-dpdf', '-painters', 'basal-raw');

close(fh);


end