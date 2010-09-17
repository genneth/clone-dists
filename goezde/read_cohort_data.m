function [ts ms binned] = read_cohort_data(file)

data = dlmread(file);

times = data(:,2);
sizes = data(:,3);

ts = unique(times);
ms = cell(numel(ts),1);
av = zeros(numel(ts),1);

for i=1:numel(ts)
    ms{i} = sizes(times==ts(i));
    av(i) = mean(ms{i});
end

xs = zeros(1,ceil(log2(max(cellfun(@max,ms))))+1);
exs = zeros(size(xs));
for i=1:numel(xs)
    xs(i) = 2^(i-1) * 1.5;
    exs(i) = 2^(i-1) * 0.5;
end

binned = zeros(numel(ts), ceil(log2(max(cellfun(@max,ms))))+1);
for i=1:numel(ts)
    for j=1:numel(ms{i})
        binned(i,floor(log2(ms{i}(j)))+1) = ...
            binned(i,floor(log2(ms{i}(j)))+1) + 1;
    end
end

colours = [
 0.996078, 0.360784, 0.027451;
% 0.996078, 0.988235, 0.0352941;
 0.541176, 0.713725, 0.027451;
 0.145098, 0.435294, 0.384314;
 0.00784314, 0.509804, 0.929412;
 0.152941, 0.113725, 0.490196;
 0.470588, 0.262745, 0.584314;
 0.890196, 0.0117647, 0.490196;
 0.905882, 0.027451, 0.129412
];

% non-scaled
fh = figure;
ah = newplot(fh);
set(ah, 'NextPlot', 'add');
for i=1:numel(ts)
    lh = herrorbar(1:numel(xs), ...
        binned(i,:) / sum(binned(i,:)), ...
        0.5*ones(1,numel(xs)));
    set(get(get(lh(1), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    for j=1:numel(lh)
        set(lh(j), 'Color', colours(i,:));
        set(lh(j), 'LineWidth', j);
    end
end

legend(ah, arrayfun(@(l)(sprintf('%.2f',l)), ts, 'UniformOutput', false));

set(ah, 'XTick', 1:numel(xs), ...
    'XTickLabel', arrayfun(@(i)(sprintf('%d-%d', 2^(i-1), 2^i-1)), 1:numel(xs), ...
    'UniformOutput', false));

set(ah, 'YLim', [0 1]);

title(ah, file, 'FontName', 'Palatino', 'FontWeight', 'bold', 'FontSize', 12);
xlabel(ah, 'clone size', ...
    'FontName', 'Palatino', 'FontWeight', 'bold', 'FontSize', 11);
set(ah, 'FontName', 'Palatino', 'FontSize', 10);

set(fh, 'PaperUnits', 'inches');
set(fh, 'PaperSize', [6 4]);
set(fh, 'PaperPosition', [0 0 6 4]);

print(fh, '-dpdf', '-painters', strcat(file, '.raw.pdf'));

% scaled
fh = figure;
ah = newplot(fh);
set(ah, 'NextPlot', 'add');
for i=1:numel(ts)
    lh = herrorbar(xs / av(i), ...
        (binned(i,:) / sum(binned(i,:))) ./ (2*exs) * av(i), ...
        exs / av(i));
    set(get(get(lh(1), 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    for j=1:numel(lh)
        set(lh(j), 'Color', colours(i,:));
        set(lh(j), 'LineWidth', j);
    end
end

legend(ah, arrayfun(@(l)(sprintf('%.2f',l)), ts, 'UniformOutput', false));

set(ah, 'XLim', [0 5], 'YLim', [10^-3 5]);
set(ah, 'YScale', 'log');

title(ah, file, 'FontName', 'Palatino', 'FontWeight', 'bold', 'FontSize', 12);
xlabel(ah, 'scaled clone size', ...
    'FontName', 'Palatino', 'FontWeight', 'bold', 'FontSize', 11);
set(ah, 'FontName', 'Palatino', 'FontSize', 10);

set(fh, 'PaperUnits', 'inches');
set(fh, 'PaperSize', [6 4]);
set(fh, 'PaperPosition', [0 0 6 4]);

print(fh, '-dpdf', '-painters', strcat(file, '.scaled.pdf'));

end
