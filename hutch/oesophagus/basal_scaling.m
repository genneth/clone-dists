function basal_scaling

% avoid warning
colours = [];

colour_palette;
oesophagus_data;

% remove the 1 week time point and collate the data together
ts = unique(ts2(ts2 ~= 1 & ts2 ~= 30/7));
data = arrayfun(@(t) collate(basal(ts2==t)), ts, 'UniformOutput', false);

% condition out singles
for i = 1:numel(ts)
%     data{i}(1+1) = 0;
end

kk = max(cellfun(@numel, data)) - 1;

tot = cellfun(@sum, data);
av = zeros(size(ts));
for i = 1:numel(ts)
    av(i) = dot(0:(numel(data{i})-1), data{i} / tot(i));
end

fh = figure;
gh = newplot(fh);
set(gh, 'NextPlot', 'add');
set(gh, 'XLim', [0 4], 'YLim', [0 1], 'YScale', 'linear');
set(gh, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(gh, 'scaled clone size', 'FontSize', 9, 'FontName', 'Helvetica', 'FontWeight', 'bold');
ylabel(gh, 'probability', 'FontSize', 9, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(gh, 'Box', 'off');
set(gh, 'YGrid', 'on', 'YMinorGrid', 'off', 'GridLineStyle', '--');
% set(gh, 'YTick', logspace(-8, 1, 10));
% set(gh, 'YTickLabel', '10 ');

% ps = clone_dist_b(r, gamma, ts*lambda, kk);
% Z = ps(0+1,:) + ps(1+1,:);
% ps = ps ./ repmat(1-Z, kk+1, 1);
% ps(0+1,:) = 0;
% ps(1+1,:) = 0;
% tav = basal_average(r, gamma, lambda, ts);
% cps = 1 - cumsum(ps, 1);
cps = [];
load basal-scaling

for i = 1:numel(ts)
%     plot(gh, (1:kk) / tav(i), cps((1+1):end,i), ...
%         '-', 'Color', colours(mod(i-1,7)+1,:), 'LineWidth', 1, 'HandleVisibility', 'off');
%     plot(gh, (1:kk) / tav(i), binoinv(0.975, tot(i), cps((1+1):end,i))/tot(i), ...
%         '--', 'Color', colours(mod(i-1,7)+1,:), 'LineWidth', 0.5, 'HandleVisibility', 'off');        
%     plot(gh, (1:kk) / tav(i), (1-binoinv(0.975, tot(i), 1-cps((1+1):end,i))/tot(i)), ...
%         '--', 'Color', colours(mod(i-1,7)+1,:), 'LineWidth', 0.5, 'HandleVisibility', 'off');        
    plot(gh, (1:(numel(data{i})-1)) / av(i), ...
        data{i}((1+1):end) * av(i) / tot(i), ...
        '^-', ...
        'MarkerEdgeColor', colours(mod(i-1,7)+1,:), ...
        'MarkerFaceColor', colours(mod(i-1,7)+1,:), ...
        'Color', colours(mod(i-1,7)+1,:), ...
        'MarkerSize', 2);

    drawnow;
end
plot(gh, linspace(0,5,50), exp(-linspace(0,5,50)), ...
    '-', 'Color', 'black', 'LineWidth', 1.0);

[~, ohs, ~, ~] = legend(gh, ...
    '3d', ...
    '10d', ...
    '3w', ...
    '6w', ...
    '3m', ...
    '6m', ...
    '12m', ...
    '{\infty} limit', ...
    'Location', 'EastOutside');

for i = 1:numel(ohs)
    if strcmp(get(ohs(i), 'Type'), 'text')
        set(ohs(i), 'Color', colours(mod(i-1,7)+1,:));
    end
end

set(fh, 'PaperUnits', 'inches');
w = 4.5; h = 2.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'basal-scaling.pdf');

close(fh);

end
