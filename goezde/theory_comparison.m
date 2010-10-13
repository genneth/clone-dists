function theory_comparison

ts1 = [1 2 4 8 16];
raft1 = [
    31	21	15	12	10;
    33	30	11	9	7;
    14	10	11	8	10;
    10	11	9	7	5;
    4	0	2	9	6;
    2	4	5	5	12;
    0	1	4	4	3;
    0	0	1	5	4;
    0	1	2	3	1;
    0	0	1	1	3;
    0	0	1	1	3;
    0	1	1	4	5;
    0	0	1	4	0;
    0	0	1	1	1;
    0	0	0	1	0;
    0	0	0	0	1;
    0	0	0	1	2;
    0	0	1	1	1;
    0	0	0	0	1;
    0	0	0	1	0;
    0	0	0	0	1;
    0	0	0	0	0;
    0	0	0	0	0;
    0	0	0	0	0;
    0	0	0	0	0;
    0	0	0	1	1;
    0	0	0	0	2;
    0	0	0	0	0;
    0	0	0	0	1];
[m1,~] = size(raft1);
raft1(1,:) = 0;
tot1 = sum(raft1);
av1 = ((1:m1) * raft1) ./ tot1;

r1 = 0.2; r2 = 0.1;
ts = arrayfun(@(s)(fzero(@(t)(average_size(r1,r2,t)-s),1)), av1);

theory1 = clone_dist(r1,r2,ts,m1);
theory1 = theory1 ./ (1 - repmat(theory1(1+1,:), m1+1, 1));
theory1(1+1,:) = 0;
% tav1 = ((0:m1) * theory1) ./ sum(theory1)

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
for i=1:numel(ts1)
    
    lo = {2};
    hi = {};
    counts = {0};
    theory = {0};
    acc = 0;
    for j = 2:m1
        acc = acc + theory1(j+1,i);
        counts{end} = counts{end} + raft1(j,i);
        theory{end} = theory{end} + theory1(j+1,i);
        if binopdf(0,tot1(i),acc) < 0.01
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
    
    plot(ah, (hi+lo)/2, (counts / tot1(i)) ./ (hi-lo+1), 's', ...
        'Color', colours(i,:), 'HandleVisibility','off');
    plot(2:m1, theory1((2:m1)+1,i), '-', ...
        'LineWidth', 1, 'Color', colours(i,:));
    errorbar((hi+lo)/2, theory ./ (hi-lo+1), ...
        (theory - (max(0, binoinv(0.025, tot1(i), theory) - 1) / tot1(i))) ./ (hi-lo+1), ...
        (binoinv(1-0.025, tot1(i), theory) / tot1(i) - theory) ./ (hi-lo+1), ...
        'Marker', 'none', 'LineStyle', 'none', 'Color', colours(i,:), 'HandleVisibility', 'off');
end

set(ah, 'FontName', 'Palatino', 'FontSize', 10);
xlabel(ah, 'Clone size', 'FontName', 'Palatino', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(ah, 'Proportion', 'FontName', 'Palatino', 'FontSize', 11, 'FontWeight', 'bold');
% set(ah, 'YScale', 'log', 'YLim', [10^-3 1]);

legend(ah, arrayfun(@(t)(sprintf('%d',t)), ts1, 'UniformOutput', false));

set(fh, 'PaperUnits', 'inches');
w = 6; h = 4;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'clone-dist');

close(fh);

end
