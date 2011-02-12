function invol

times = [3.5/7 1 2 4 8 12]; % weeks
data = {
    sparse([
        0	0;
        32	14;
        5	2
    ]), ...
    sparse([
        0	0	0;
        9	4	2;
        2	1	1;
        0	2	0;
        1	0	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0;
        9	7	8	2	0	4	2;
        10	4	4	0	0	0	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0;
        4	11	1	4	3	5	0	0	0	1;
        10	1	1	1	2	1	0	0	0	0;
        0	1	2	1	1	1	0	1	0	0;
        1	0	1	0	0	2	0	0	0	0;
        0	1	0	0	0	0	0	0	0	0;
        0	1	1	0	0	0	0	0	0	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0	0	0;
        1	3	2	1	3	1	0	0	0	1	0	0;
        0	1	0	2	1	0	0	1	0	2	1	0;
        0	0	1	1	1	1	1	0	1	0	0	0;
        1	0	1	0	0	0	0	0	1	0	0	0;
        1	0	0	0	0	1	0	0	0	0	0	0;
        0	1	1	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	1;
        0	0	0	0	0	0	0	1	0	0	0	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        1	3	2	1	3	1	0	0	0	1	0	0	0	0	0	0	0	0;
        0	1	0	2	1	0	0	1	0	2	1	0	0	0	0	0	0	0;
        0	0	1	1	1	2	1	0	1	0	0	0	0	0	0	0	0	0;
        1	0	1	3	1	1	0	0	1	0	0	0	0	0	0	0	0	0;
        1	0	3	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0;
        0	1	2	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0;
        0	0	0	1	1	0	0	1	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ])};

% if we're doing inference...
if 0
    % remove single cells from consideration --- don't know if we're inducing
    % post-mitotic cells
    for i = 1:numel(times)
        data{i}(1+1,0+1) = 0;
    end

    samples = sample3_no_shed(...
        @()(random('beta', 1, 3) / 2), ...
        @()(random('logn', log(2), log(1.8))), ...
        @()(random('logn', 0, log(2))), ...
        [], {}, times, data, 1000, 30);

    save invol_samples.mat samples
else
    % if we're doing fitting tests...

    % parameters
    r = 0.1; gamma = 2; mu = 1/3; lambda = 1;

    % prepare the theoretical predictions
    max_size = max(cell2mat(...
       cellfun(@(d)(size(d)'), ...
           data, 'UniformOutput', false))');
    theory_raw = clone_dist_bs_shed(r, gamma, mu, lambda*times, max_size(1), max_size(2));
    basal_raw = clone_dist_b(r, gamma, lambda*times, max_size(1));
    extinction_prob = basal_raw(0+1, :);

    theory = theory_raw; theory(0+1,:,:) = 0;
    basal = basal_raw; basal(0+1,:) = 0;
    for i = 1:numel(times)
        theory(:,:,i) = theory(:,:,i) / (1 - extinction_prob(i));
        basal(:,i) = basal(:,i) / (1 - extinction_prob(i));
    end
                
    % first up, our home-grown "how many experiments" test
    f = fopen('invol-theory-comparison.txt', 'w');
    fclose(f);
    for i = 1:numel(times)
        [m,n] = size(data{i});
        chances = binocompare(data{i}, theory(1:m,1:n,i));
        dlmwrite('invol-theory-comparison.txt', 1 ./ chances, '-append', ...
            'delimiter', '\t', 'roffset', 2);
    end
    
    % Ben's battery: clone size distributions for basal, suprabasal and
    % total; in both linear and log scale, and both pdf and cdf. Here we
    % want to exploit the fact that we have 6 time points, and that 6 =
    % 2x3, to layout the graphs
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
    
    print(fh, '-depsc2', '-painters', 'invol-theory-comparison-basal-pdf-linear');
    close(fh);

end

end