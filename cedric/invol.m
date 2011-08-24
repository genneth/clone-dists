function invol

times = [3.5/7 1 2 4 8 12 24]; % weeks
data = {
    sparse([
        0	0
        52	27
        7	2
    ]), ...
    sparse([
        0	0	0	0
        20	11	6	2
        4	3	1	1
        0	2	0	0
        1	0	0	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0
        10	7	9	3	2	4	2
        11	6	4	0	0	0	0
        0	1	0	0	0	0	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0
        5	16	3	12	5	6	0	0	0	1
        14	2	7	2	3	1	0	2	0	0
        3	3	4	2	1	1	0	1	1	0
        1	1	1	0	0	2	0	0	0	0
        0	1	0	0	1	0	0	1	0	0
        0	1	2	0	0	1	0	0	0	0
        0	0	1	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	1	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0	0	0
        1	3	2	1	3	1	0	0	0	1	0	0
        0	1	0	2	1	0	0	1	0	2	1	0
        0	0	1	1	1	1	1	0	1	0	0	0
        1	0	1	0	0	0	0	0	1	0	0	0
        1	0	0	0	0	1	0	0	0	0	0	0
        0	1	1	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	1
        0	0	0	0	0	0	0	1	0	0	0	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	8	2	6	4	0	0	0	2	0	0	0	0	0	0	0	0	0
        0	0	4	5	2	0	1	2	0	0	0	0	1	1	0	0	0	2
        0	0	2	4	2	0	2	4	0	0	0	0	0	0	0	0	0	0
        0	2	2	4	0	0	2	0	0	2	0	0	0	0	0	0	0	0
        0	0	0	4	0	0	0	0	0	2	0	0	0	0	0	0	0	0
        0	2	0	0	1	0	2	0	4	0	1	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0
        0	0	0	0	0	0	0	0	0	2	1	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
        0	0	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1
        0	0	1	0	1	1	1	0	0	0	0	0	0	0	0	0	0
        0	0	1	0	0	0	1	0	0	0	0	0	0	0	0	0	0
        0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
    ])};

if 1
    % if we're doing inference...
    % remove single cells from consideration --- don't know if we're inducing
    % post-mitotic cells
    for i = 1:numel(times)
        data{i}(1+1,0+1) = 0;
    end

    stream = RandStream('mt19937ar','seed',sum(100*clock));
    RandStream.setDefaultStream(stream);
    sample3_shed(...
        @()(random('beta', 1.1, 8.9) / 2), ... % r
        @()(random('logn', log(3), log(1.8))), ... % gamma
        @()(random('logn', 0, log(2))), ... % lambda
        5/4, ... % suprabasal:basal ratio (m)
        times, data, 1000, ...
        'invol_samples.mat');

else
    % if we're doing fitting tests...

    % parameters
    r = 0.4; gamma = 3; mu = (3/4)/(5/4); lambda = 1; delta = -1/13;

    % prepare the theoretical predictions
    max_size = max(cell2mat(...
       cellfun(@(d)(size(d)'), ...
           data, 'UniformOutput', false))');
    theory_raw = inverse_z_transform2(@(x,y) generating_function3_shed_unbalance(r, delta, gamma, mu, lambda*times, x, x, y), max_size(1), max_size(2), 1e-4, 1e-4);
    basal_raw = inverse_z_transform(@(x) generating_function2_unbalance(r, delta, gamma, lambda*times, x, x), max_size(1), 1e-7, 1e-7);
    extinction_prob = basal_raw(0+1, :);

    theory = theory_raw; theory(0+1,:,:) = 0;
    basal = basal_raw; basal(0+1,:) = 0;
    for i = 1:numel(times)
        theory(:,:,i) = theory(:,:,i) / (1 - extinction_prob(i));
        basal(:,i) = basal(:,i) / (1 - extinction_prob(i));
    end
    
    save('invol', 'theory', 'basal');
%     theory = []; basal = [];
%     load('invol');
                
    % first up, our home-grown "how many experiments" test
%     f = fopen('invol-theory-comparison.txt', 'w');
%     fclose(f);
%     for i = 1:numel(times)
%         [m,n] = size(data{i});
%         chances = binocompare(data{i}, theory(1:m,1:n,i));
%         dlmwrite('invol-theory-comparison.txt', 1 ./ chances, '-append', ...
%             'delimiter', '\t', 'roffset', 2);
%     end
    
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
    
    print(fh, '-depsc2', '-painters', 'invol-theory-comparison-total-pdf-linear');
    close(fh);
end

end
