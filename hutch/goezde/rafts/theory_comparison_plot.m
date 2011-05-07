function theory_comparison_plot(ts, raft, maxm, file)

[m,~] = size(raft);
raft(1,:) = 0; % condition out single cell clones
tot = sum(raft);
av = ((1:m) * raft) ./ tot;

r1 = 0.2; r2 = 0.1;
times = arrayfun(@(s)(fzero(@(t)(average_size(r1,r2,t)-s),1)), av);

ps = inverse_z_transform(@(z) generating_function_unbalanced(r1,r2,times,z,z), m, 1e-6, 1e-6);
ps = ps((1+1):end,:);
ps = ps ./ (1 - repmat(ps(1,:), m, 1));
ps(1,:) = 0;

fh = figure;
set(fh, 'PaperUnits', 'inches');
w = 7; h = 4.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);
margin_left = 0.1; margin_bottom = 0.1;
hpadding = 0.03/w*h; vpadding = 0.03;

for i = 1:numel(times)
    i_ = mod(i-1,3); j_ = floor((i-1)/3);
    ah = subplot('Position', [
        margin_left+i_/3*(1-margin_left-2*hpadding)+i_*hpadding
        margin_bottom+(1-j_)/2*(1-margin_bottom-2*vpadding)+(1-j_)*vpadding
        (1-margin_left-2*hpadding)/3 
        (1-margin_bottom-2*vpadding)/2]);
    set(ah, 'NextPlot', 'add');
    set(ah, 'FontName', 'Helvetica', 'FontSize', 10);
    set(ah, 'Box', 'on');
    
    lo = {2};
    hi = {};
    counts = {0};
%     theory = {0};
    acc = 0;
    for j = 2:m
        acc = acc + ps(j,i);
        counts{end} = counts{end} + raft(j,i);
%         theory{end} = theory{end} + ps(j,i);
        if binopdf(0,tot(i),acc) < 0.01
            hi{end+1} = j;
            lo{end+1} = j+1;
            counts{end+1} = 0;
%             theory{end+1} = 0;
            acc = 0;
        end
    end
    lo = [lo{1:(end-1)}];
    hi = cell2mat(hi);
    counts = [counts{1:(end-1)}];
%     theory = [theory{1:(end-1)}];
    
    plot(ah, (hi+lo)/2, (counts / tot(i)) ./ (hi-lo+1), 's');
    plot(ah, 2:m, ps(2:m,i), '-');
    plot(ah, 2:m, 1 - binoinv(1-0.158655254, tot(i), 1-ps(2:m,i)) / tot(i), '--');
    plot(ah, 2:m, binoinv(1-0.158655254, tot(i), ps(2:m,i)) / tot(i), '--');
    
    set(ah, 'YLim', [0 0.6]);
    if j_ == 0 % top row
        set(ah, 'XTickLabel', {});
    else
    end
    if i_ == 0 % left column
        ylabel(ah, 'proportion', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    else
        set(ah, 'YTickLabel', {});
    end
    if i == 5
        xlabel(ah, 'clone size', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    end

    set(ah, 'XLim', [0 maxm]);
    
    th = title(ah, sprintf('%.1f days\n%.1f divisions',ts(i), times(i)), ...
        'FontName', 'Helvetica', 'FontSize', 10, 'Units', 'normalized');
    ex = get(th, 'Extent');
    set(th, 'HorizontalAlignment', 'right');
    set(th, 'VerticalAlignment', 'top');
    set(th, 'Position', [1-ex(4)/4 1-ex(4)/4]);
    
    % set(ah, 'YScale', 'log', 'YLim', [10^-3 1]);
end


% legend(ah, arrayfun(@(t)(sprintf('%.1f days',t)), ts, 'UniformOutput', false));

print(fh, '-dpdf', '-painters', file);

close(fh);

end
