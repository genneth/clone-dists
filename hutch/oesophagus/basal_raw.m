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
% average size
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
ts = linspace(0, 60, 100);
% theory = basal_average(0.1, 0.65 / (1-0.65), 1.87, ts);
% save basal-average theory
theory = [];
load basal-average
plot(ah, ts, theory, '-', 'Color', colours(2,:));
p = get(ah, 'Position');
p = [0.03 p(2) p(3)+p(1)-0.03 p(4)];
set(ah, 'Position', p);

% clone density
times = [3 6 12 26 52];
densities = [
    0.458515284	0.611353712	0.742358079
    0.043668122	0.349344978	0.524017467
    0.262008734	0.305676856	0.131004367
    0.021834061	0.305676856	0.218340611
    0.139737991	0.043668122	0.021834061
    ];
mds = mean(densities, 2);
mdd = sqrt(var(densities, 1, 2) * 3 / 2);

ah = axes('OuterPosition', [0.55 0.7 0.45 0.3]);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 8);
set(ah, 'XLim', [0 55]);
set(ah, 'XGrid', 'on', 'XMinorGrid', 'off');
set(ah, 'XTick', [4 12 26 52]);
set(ah, 'XTickLabel', {'1m', '3m', '6m', '1y'});
set(ah, 'YGrid', 'on', 'YMinorGrid', 'off');
ylabel(ah, 'clone density', 'FontSize', 8, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(ah, 'YLim', [0 0.8], 'YTick', [0 0.3 0.6]);
plot(ah, times, mds, '+', 'MarkerSize', 2.0, 'Color', colours(1,:));
for i = 1:numel(times)
    plot(ah, [times(i) times(i)], mds(i) - mdd(i)*[-1 1], '-', 'Color', colours(1,:));
end
ts = linspace(0, 60, 100);
extinction = generating_function2(0.10, 0.65/(1-0.65), 1.87*ts, 0, 0);
plot(ah, ts, (1-extinction) * 0.8, '-', 'Color', colours(2,:));
p = get(ah, 'Position');
p = [0.58 p(2) p(3)+p(1)-0.58 p(4)];
set(ah, 'Position', p);

% labelled percentage
labelling = [
    1.4819869	1.671033479	2.176744875
    0.191467921	1.470160116	1.871490954
    1.505115864	1.344978166	0.741293002
    0.175722204	2.170614441	2.793475469
    3.28625207	0.894267398	0.432751092
    ];
ls = mean(labelling, 2);
ld = sqrt(var(labelling, 1, 2) / 2);
lm = mean(ls);
lmd = sqrt(var(reshape(labelling, numel(labelling), 1)) / (numel(labelling) - 1));
ah = axes('OuterPosition', [0.0 0.4 0.45 0.3]);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 8);
set(ah, 'XLim', [0 55]);
set(ah, 'XGrid', 'on', 'XMinorGrid', 'off');
set(ah, 'XTick', [4 12 26 52]);
set(ah, 'XTickLabel', {'1m', '3m', '6m', '1y'});
set(ah, 'YAxisLocation', 'right', 'YGrid', 'on', 'YMinorGrid', 'off');
ylabel(ah, '% labelled cells', 'Rotation', -90, 'FontSize', 8, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(ah, 'YLim', [0 2.5], 'YTick', [0 1 2]);
plot(ah, times, ls, '+', 'MarkerSize', 2.0, 'Color', colours(1,:));
for i = 1:numel(times)
    plot(ah, [times(i) times(i)], ls(i) - ld(i)*[-1 1], '-', 'Color', colours(1,:));
end
ts = linspace(0, 60, 100);
lm = repmat(lm, size(ts));
lmd = repmat(lmd, size(ts));
plot(ah, ts, lm, '-', 'Color', colours(2,:));
plot(ah, ts, lm+lmd, '-', 'Color', colours(2,:));
plot(ah, ts, lm-lmd, '-', 'Color', colours(2,:));
p = get(ah, 'Position');
p = [0.03 p(2) p(3)+p(1)-0.03 p(4)];
set(ah, 'Position', p);

print(fh, '-dpdf', '-painters', 'basal-raw');

close(fh);


end