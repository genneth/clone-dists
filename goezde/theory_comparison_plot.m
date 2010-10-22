function theory_comparison_plot(ts, raft, maxm, file)

[m,~] = size(raft);
raft(1,:) = 0;
tot = sum(raft);
av = ((1:m) * raft) ./ tot;

r1 = 0.2; r2 = 0.1;
times = arrayfun(@(s)(fzero(@(t)(average_size(r1,r2,t)-s),1)), av);

ps = clone_dist(r1,r2,times,m);
ps = ps ./ (1 - repmat(ps(1+1,:), m+1, 1));
ps(1+1,:) = 0;

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

fh = figure;
ah = newplot(fh);
set(ah, 'NextPlot', 'add');
for i=1:numel(ts)
    lo = {2};
    hi = {};
    counts = {0};
    theory = {0};
    acc = 0;
    for j = 2:m
        acc = acc + ps(j+1,i);
        counts{end} = counts{end} + raft(j,i);
        theory{end} = theory{end} + ps(j+1,i);
        if binopdf(0,tot(i),acc) < 0.01
            hi{end+1} = j;
            lo{end+1} = j+1;
            counts{end+1} = 0;
            theory{end+1} = 0;
            acc = 0;
        end
    end
    lo = [lo{1:(end-1)}];
    hi = cell2mat(hi);
    counts = [counts{1:(end-1)}];
    theory = [theory{1:(end-1)}];
%     disp([lo;hi;counts;theory]);
    
    plot(ah, (hi+lo)/2, (counts / tot(i)) ./ (hi-lo+1), 's', ...
        'Color', colours(i,:), 'HandleVisibility','off');
    plot(2:m, ps((2:m)+1,i), '-', ...
        'LineWidth', 1, 'Color', colours(i,:));
    errorbar((hi+lo)/2, theory ./ (hi-lo+1), ...
        (theory - (max(0, binoinv(0.025, tot(i), theory) - 1) / tot(i))) ./ (hi-lo+1), ...
        (binoinv(1-0.025, tot(i), theory) / tot(i) - theory) ./ (hi-lo+1), ...
        'Marker', 'none', 'LineStyle', 'none', 'Color', colours(i,:), 'HandleVisibility', 'off');
end

set(ah, 'FontName', 'Palatino', 'FontSize', 10);
xlabel(ah, 'Clone size', 'FontName', 'Palatino', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(ah, 'Proportion', 'FontName', 'Palatino', 'FontSize', 11, 'FontWeight', 'bold');
set(ah, 'XLim', [0 maxm]);
% set(ah, 'YScale', 'log', 'YLim', [10^-3 1]);

legend(ah, arrayfun(@(t)(sprintf('%.1f days',t)), ts, 'UniformOutput', false));

set(fh, 'PaperUnits', 'inches');
w = 6; h = 4;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', file);

close(fh);

end
