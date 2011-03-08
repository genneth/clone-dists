function k14

times = [3.5/7 1 2 4 8]; % weeks
data = {
    sparse([
0	0	0	0
89	15	5	1
67	6	4	0
19	2	1	0
6	1	0	0
1	0	0	0
    ]), ...
    sparse([
0	0	0
15	6	0
13	3	2
7	1	1
1	2	0
1	1	0
    ]), ...
    sparse([
0	0	0	0	0	0	0	0	0	0	0
2	6	1	1	0	0	0	0	0	0	0
7	7	0	0	0	0	0	0	0	0	0
1	0	1	0	0	0	0	0	0	0	0
1	2	1	1	0	0	0	0	0	0	0
0	1	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	1
0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	1
    ]), ...
    sparse([
0	0	0	0	0	0	0	0	0	0	0	0	0
11	11	6	2	1	1	0	1	0	0	0	0	0
8	7	7	2	2	1	1	0	0	0	0	0	0
3	3	3	2	2	0	0	0	0	0	0	0	0
1	2	4	0	0	0	0	0	0	0	0	0	0
1	0	1	0	0	0	2	0	1	0	0	0	0
0	0	1	0	0	1	0	0	0	0	0	0	1
0	0	0	0	1	1	0	1	0	0	0	0	0
0	0	0	0	1	0	0	0	0	0	0	0	0
0	0	0	1	0	0	0	0	0	0	0	0	0
    ]), ...
    sparse([
0	0	0	0	0	0	0	0	0
3	0	3	4	1	0	0	0	0
10	7	7	3	3	0	0	0	0
2	5	1	2	2	1	0	0	0
0	0	2	5	1	1	0	0	0
1	1	3	1	0	0	0	0	0
0	0	1	2	0	0	0	0	0
0	0	0	0	0	0	0	0	0
0	0	0	1	1	0	0	0	0
0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	1
    ])};

% parameters
r = 0.1; gamma = 2; mu = 1/3; lambda = 1; theta = 0.5;

% prepare the theoretical predictions
max_size = max(cell2mat(...
   cellfun(@(d)(size(d)'), ...
       data, 'UniformOutput', false))');
theory_raw = clone_dist_bs_shed(r, gamma, mu, lambda*times, max_size(1), max_size(2));
basal_raw = clone_dist_b(r, gamma, lambda*times, max_size(1));
extinction_prob = basal_raw(0+1, :);

% we need to account for theta dependence
theory = theory_raw; theory(0+1,:,:) = 0;
basal = basal_raw; basal(0+1,:) = 0;
for i = 1:numel(times)
    theory(:,:,i) = theory(:,:,i) / (1 - extinction_prob(i));
    basal(:,i) = basal(:,i) / (1 - extinction_prob(i));
end
theory = theory * (1-theta);
theory(2:end,:,:) = theory(2:end,:,:) + theta * theory_raw(1:(end-1),:,:);
basal = basal * (1-theta);
basal(2:end,:) = basal(2:end,:) + theta * basal_raw(1:(end-1),:);

% first up, our home-grown "how many experiments" test
f = fopen('k14-theory-comparison.txt', 'w');
fclose(f);
for i = 1:numel(times)
    [m,n] = size(data{i});
    chances = binocompare(data{i}, theory(1:m,1:n,i));
    dlmwrite('k14-theory-comparison.txt', 1 ./ chances, '-append', ...
        'delimiter', '\t', 'roffset', 2);
end

% Ben's battery: clone size distributions for basal, suprabasal and
% total; in both linear and log scale, and both pdf and cdf. Although we
% don't have 6 graphs, 5 will do :-)
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
        margin_bottom+(1-j_)/2*(1-margin_bottom-vpadding)+(1-j_)*vpadding
        (1-margin_left-2*hpadding)/3 
        (1-margin_bottom-vpadding)/2]);
    set(ah, 'NextPlot', 'add');
    set(ah, 'FontName', 'Helvetica', 'FontSize', 10);
    set(ah, 'Box', 'on');
    total = sum(sum(data{i}));
    exp_wz = sum(data{i},2); exp_wz(max_size(1)+1) = 0;
    plot(ah, 1:(max_size(1)-1), exp_wz(2:max_size(1)) / total, 's');
    plot(ah, 1:(max_size(1)-1), basal(2:max_size(1),i), '-');
    plot(ah, 1:(max_size(1)-1), binoinv(0.158655254, total, basal(2:max_size(1),i)) / total, '--');
    plot(ah, 1:(max_size(1)-1), binoinv(1-0.158655254, total, basal(2:max_size(1),i)) / total, '--');
    set(ah, 'XLim', [0 max_size(1)]);
    if i <= 3
        set(ah, 'YLim', [0 1]);
        set(ah, 'XTickLabel', {});
    else
        set(ah, 'YLim', [0 0.5]);
    end
    if mod(i,3)==1
        ylabel(ah, 'probability', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    else
        set(ah, 'YTickLabel', {});
    end
    if i == 5
        xlabel(ah, 'clone size', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    end
end

print(fh, '-depsc2', '-painters', 'k14-theory-comparison-basal-pdf-linear');
close(fh);

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
        margin_bottom+(1-j_)/2*(1-margin_bottom-vpadding)+(1-j_)*vpadding
        (1-margin_left-2*hpadding)/3 
        (1-margin_bottom-vpadding)/2]);
    set(ah, 'NextPlot', 'add');
    set(ah, 'FontName', 'Helvetica', 'FontSize', 10);
    set(ah, 'Box', 'on');
    [m,n] = size(data{i});
    total = sum(sum(data{i}));
    exp_wz = arrayfun(@(d) sum(diag(flipud(full(data{i})), d)), -(m-1):(n-1));
    exp_wz(sum(max_size)) = 0;
    plot(ah, 1:(sum(max_size)-2), exp_wz(2:(sum(max_size)-1)) / total, 's');
    tot_wz = arrayfun(@(d) sum(diag(flipud(theory(1:m,1:n,i)), d)), -(m-1):(n-1));
    tot_wz(sum(max_size)) = 0;
    plot(ah, 1:(sum(max_size)-2), tot_wz(2:(sum(max_size)-1)), '-');
    plot(ah, 1:(sum(max_size)-2), binoinv(0.158655254, total, tot_wz(2:(sum(max_size)-1))) / total, '--');
    plot(ah, 1:(sum(max_size)-2), binoinv(1-0.158655254, total, tot_wz(2:(sum(max_size)-1))) / total, '--');
    set(ah, 'XLim', [0 20]);
    if i <= 3
        set(ah, 'YLim', [0 1]);
        set(ah, 'XTickLabel', {});
    else
        set(ah, 'YLim', [0 0.5]);
    end
    if mod(i,3)==1
        ylabel(ah, 'probability', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    else
        set(ah, 'YTickLabel', {});
    end
    if i == 5
        xlabel(ah, 'clone size', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    end
end

print(fh, '-depsc2', '-painters', 'k14-theory-comparison-total-pdf-linear');
close(fh);

end
