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

    max_size = max(cell2mat(...
       cellfun(@(d)(size(d)'), ...
           data, 'UniformOutput', false))');
    theory_raw = clone_dist_bs_shed(r, gamma, mu, lambda*times, max_size(1), max_size(2));
    extinction_prob = generating_function2(r, gamma, lambda*times, 0, 0);

    theory = theory_raw;
    theory(0+1,:,:) = 0;
    f = fopen('invol-theory-comparison.txt', 'w');
    fclose(f);
    for i = 1:numel(times)
        theory(:,:,i) = theory(:,:,i) / (1 - extinction_prob(i));
        
        [m,n] = size(data{i});
        chances = binocompare(data{i}, theory(1:m,1:n,i));
        dlmwrite('invol-theory-comparison.txt', 1 ./ chances, '-append', ...
            'delimiter', '\t', 'roffset', 2);
    end
end

end