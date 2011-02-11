function oesophagus_basal_scaling

ts = [3/7 10/7 3 6 12 26 52];

data{1} = [0 142 71 8 3];
data{2} = [0 111 109 35 15 4 2];
data{3} = [0 140 79 38 25 11 4 1 1 0 0 1];
data{4} = [0 73 60 44 26 22 10 3 4 0 2 4 0 2 1 0 0 0 0 0 1 0 0 0 0 0 0 0 1];
data{5} = [0 87 87 50 52 26 24 14 8 9 8 6 7 4 3 3 2 2 0 0 2 0 1 0 0 0 1 0 0 1];
data{6} = [0 40 59 43 18 27 23 16 12 14 10 11 5 8 7 6 3 4 3 3 4 1 2 1 3 ...
    4 0 1 0 0 1 0 0 1 1 5 1 2 0 1 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 1 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1];
data{7} = [0 12 18 15 14 6 11 8 6 9 9 3 4 4 2 4 4 2 4 2 2 3 4 1 5 0 2 4 ...
    3 1 3 3 3 2 4 2 1 1 1 0 1 1 1 1 1 1 0 2 1 1 0 0 0 0 1 3 0 0 1 0 2 0 ...
    1 1 0 0 0 0 0 0 2 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ...
    1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
    0 0 0 0 0 0 0 0 0 0 0 1];

kk = max(cellfun(@numel, data)) - 1;

tot = cellfun(@sum, data);
av = zeros(size(ts));
for i = 1:numel(ts)
    av(i) = dot(0:(numel(data{i})-1), data{i} / tot(i));
end

lambda = 0.7676;
r = 0.1926;
rho = 0.5164;
gamma = rho / (1-rho);
tau = rho / (r * lambda);

fh = figure;
gh = newplot(fh);
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
set(gh, 'NextPlot', 'add');
set(gh, 'XLim', [0 5], 'YLim', [1e-8 10], 'YScale', 'log');
set(gh, 'YAxisLocation', 'right');
set(gh, 'Box', 'on');
set(gh, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(gh, 'scaled clone size', 'FontSize', 9, 'FontName', 'Helvetica', 'FontWeight', 'bold');
set(gh, 'YTick', logspace(-8, 1, 10));
% set(gh, 'YTickLabel', '10 ');

ps = clone_dist_b(r, gamma, ts*lambda, kk);
Z = ps(0+1,:);
ps(0+1,:) = 0;
for i = 1:numel(ts);
    ps(:,i) = ps(:,i) ./ (1-Z(i));
end
tav = (0:kk) * ps;
% plot(newplot(figure), ts, tav, ts, av);

% bin data according to theoretical probs
lo = {}; hi = {}; pb = {}; db = {};
for i = 1:numel(ts)
    lo{i} = {0 1}; hi{i} = {0}; pb{i} = {0}; db{i} = {0};
    acc = 0;
    for j = 1:numel(ps(:,i))
        acc = acc + ps(j,i);
        if (1-acc) ^ tot(i) < 0.01
            hi{i}{end+1} = j-1;
            pb{i}{end+1} = sum(ps((lo{i}{end}+1):(hi{i}{end}+1),i));
            if lo{i}{end} <= numel(data{i}) && hi{i}{end} <= numel(data{i})
                db{i}{end+1} = sum(data{i}((lo{i}{end}+1):(hi{i}{end}+1)));
            elseif lo{i}{end} <= numel(data{i}) && hi{i}{end} > numel(data{i})
                db{i}{end+1} = sum(data{i}((lo{i}{end}+1):end));
            else
                db{i}{end+1} = 0;
            end
            lo{i}{end+1} = j;
            acc = 0;
        end
    end

    lo{i} = cell2mat({lo{i}{1:(end-1)}});
    hi{i} = cell2mat(hi{i});
    pb{i} = cell2mat(pb{i});
    db{i} = cell2mat(db{i});
    
end

for i = 1:numel(ts)
    plot(gh, (1:kk) / tav(i), ps(2:end,i) * tav(i) * 10^(i-7), ...
        '-', 'Color', colours(i,:), 'LineWidth', 1, 'HandleVisibility', 'off');
%     plot(gh, (lo{i}+hi{i}) / 2 / tav(i), ...
%         pb{i} ./ (hi{i}-lo{i}+1) * tav(i) * 10^(i-8), ...
%         '-', 'Color', colours(i,:), 'LineWidth', 1, 'HandleVisibility', 'off');
    
%     error = 0.317310507863; % 2-tailed
%     lower = arrayfun(@(p)(max(0,binoinv(error/2, tot(i), p)-1)), pb{i});
%     upper = arrayfun(@(p)(min(tot(i),binoinv(1-error/2, tot(i), p))), pb{i});
%     plot(gh, (lo{i}+hi{i}) / 2 / tav(i), ...
%                   lower ./ (hi{i}-lo{i}+1) / tot(i) * tav(i) * 10^(i-8), '--', ...
%                   (lo{i}+hi{i}) / 2 / tav(i), ...
%                   upper ./ (hi{i}-lo{i}+1) / tot(i) * tav(i) * 10^(i-8), '--', ...
%             'Color', colours(i,:), 'LineWidth', 1, 'HandleVisibility', 'off');

%     plot(gh, (lo{i}+hi{i}) / 2 / av(i), ...
%         db{i} ./ (hi{i}-lo{i}+1) / tot(i) * av(i) * 10^(i-8), '^', ...
%         'MarkerEdgeColor', colours(i,:), ...
%         'MarkerFaceColor', colours(i,:), ...
%         'MarkerSize', 2);
    errorbar(gh, (lo{i}+hi{i}) / 2 / av(i), ...
        db{i} ./ (hi{i}-lo{i}+1) / tot(i) * av(i) * 10^(i-7), ...
        sqrt(db{i}) ./ (hi{i}-lo{i}+1) / tot(i) * av(i) * 10^(i-7), ...
        '^', ...
        'MarkerEdgeColor', colours(i,:), ...
        'MarkerFaceColor', colours(i,:), ...
        'Color', colours(i,:), ...
        'MarkerSize', 2);

    drawnow;
end
plot(gh, linspace(0,5,50), exp(-linspace(0,5,50)) * 10, ...
    '-', 'Color', 'black', 'LineWidth', 1.5);

[~, ohs, ~, ~] = legend(gh, ...
    '3d', ...
    '10d', ...
    '3w', ...
    '6w', ...
    '3m', ...
    '6m', ...
    '12m', ...
    '{\infty} limit', ...
    'Location', 'NorthOutside');

for i = 1:numel(ohs)
    if strcmp(get(ohs(i), 'Type'), 'text')
        set(ohs(i), 'Color', colours(i,:));
    end
end

set(fh, 'PaperUnits', 'inches');
w = 4; h = 6;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'oesophagus-basal-scaling.pdf');

close(fh);

end
