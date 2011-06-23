function basal_comparison

% avoid warning
colours = [];

colour_palette;
oesophagus_data;

% remove the 1 week time point and collate the data together
ts = unique(ts2(ts2 ~= 1 & ts2 ~= 30/7));
data = arrayfun(@(t) collate(basal(ts2==t)), ts, 'UniformOutput', false);

% condition out singles
for i = 1:numel(ts)
     data{i}(1+1) = 0;
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
set(gh, 'XLim', [0 4], 'YLim', [1e-8 10], 'YScale', 'log');
set(gh, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(gh, 'scaled clone size', 'FontSize', 9, 'FontName', 'Helvetica', 'FontWeight', 'bold');
ylabel(gh, 'scaled probability', 'FontSize', 9, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(gh, 'Box', 'off');
set(gh, 'YGrid', 'on', 'YMinorGrid', 'off');
set(gh, 'YAxisLocation', 'right');
% set(gh, 'YTick', logspace(-8, 1, 10));
% set(gh, 'YTickLabel', '10 ');

% ps = clone_dist_b(0.10, 0.65/(1-0.65), ts*1.87, kk);
% Z = ps(0+1,:) + ps(1+1,:);
% ps = ps ./ repmat(1-Z, kk+1, 1);
% ps(0+1,:) = 0;
% ps(1+1,:) = 0;
% tav = basal_average(0.10, 0.65/(1-0.65), 1.87, ts);
% save basal-scaling ps tav
ps = []; tav = [];
load basal-scaling

for i = 1:numel(ts)
    % bin data "intelligently". *gulp*
    lo = [2]; hi = []; theory = []; experiment = [];
    for j = 2:kk
        if binopdf(0, tot(i), sum(ps((lo(end):j)+1,i))) < 0.01
            hi(end+1) = j;
            theory(end+1) = sum(ps((lo(end):j)+1,i));
            if lo(end) < numel(data{i})-1
                experiment(end+1) = sum(data{i}((lo(end):min(j, numel(data{i})-1))+1));
            else
                experiment(end+1) = 0;
            end
            lo(end+1) = j+1;
        end
    end
    lo = lo(1:(end-1)); % trim the last element
%     disp(lo);
%     disp(hi);
%     disp(experiment ./ (hi-lo+1));
%     disp(theory ./ (hi-lo+1));
    
    plot(gh, (lo+hi)/2 / av(i), ...
        experiment ./ (hi-lo+1) * av(i) / tot(i) * 10^(i-7), ...
        '^', ...
        'MarkerEdgeColor', colours(mod(i-1,7)+1,:), ...
        'MarkerFaceColor', colours(mod(i-1,7)+1,:), ...
        'Color', colours(mod(i-1,7)+1,:), ...
        'MarkerSize', 2);
    
    plot(gh, (2:kk) / tav(i), ...
        ps((2+1):end,i) * tav(i) * 10^(i-7), ...
        '-', 'Color', colours(mod(i-1,7)+1,:), 'LineWidth', 0.4, ...
        'HandleVisibility', 'off');
    plot(gh, (lo+hi)/2 / tav(i), ...
        binoinv(0.842, tot(i), theory ./ (hi-lo+1))/tot(i) * tav(i) * 10^(i-7), ...
        '--', 'Color', colours(mod(i-1,7)+1,:), 'LineWidth', 0.4, ...
        'HandleVisibility', 'off');
    plot(gh, (lo+hi)/2 / tav(i), ...
        (1-binoinv(0.842, tot(i), 1 - theory ./ (hi-lo+1))/tot(i)) * tav(i) * 10^(i-7), ...
        '--', 'Color', colours(mod(i-1,7)+1,:), 'LineWidth', 0.4, ...
        'HandleVisibility', 'off');
    % the upper error line when it gets to 1
    plot(gh, (0:5), ...
        (1/tot(i) * tav(i) * 10^(i-7)) * ones(1,6), ...
        '--', 'Color', colours(mod(i-1,7)+1,:), 'LineWidth', 0.4, ...
        'HandleVisibility', 'off');

    drawnow;
end
plot(gh, linspace(0,5,50), exp(-linspace(0,5,50))*10, ...
    '-', 'Color', 'black', 'LineWidth', 1.0);

[~, ohs, ~, ~] = legend(gh, ...
    '3d', ...
    '10d', ...
    '3w', ...
    '6w', ...
    '3m', ...
    '6m', ...
    '12m', ...
    'limit', ...
    'Location', 'EastOutside');

for i = 1:numel(ohs)
    if strcmp(get(ohs(i), 'Type'), 'text')
        set(ohs(i), 'Color', colours(mod(i-1,7)+1,:));
    end
end

set(fh, 'PaperUnits', 'inches');
w = 4.5; h = 4;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'basal-comparison.pdf');

close(fh);

end
