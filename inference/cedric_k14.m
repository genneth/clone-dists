function cedric_k14

times = [3.5/7 1 2 4 8]; % weeks
data = {
    sparse([
        0	0	0	0;
        89	15	5	1;
        67	6	4	0;
        19	2	1	0;
        6	1	0	0;
        1	0	0	0
    ]), ...
    sparse([
        0	0	0;
        15	6	0;
        13	3	2;
        7	1	1;
        1	2	0;
        1	1	0
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0	0;
        2	6	1	1	0	0	0	0	0	0	0;
        7	7	0	0	0	0	0	0	0	0	0;
        1	0	1	0	0	0	0	0	0	0	0;
        1	2	1	1	0	0	0	0	0	0	0;
        0	1	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	1;
        0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	1
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        11	11	6	2	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        8	7	7	2	2	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        3	3	3	2	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        1	2	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        1	0	1	0	0	0	2	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	1	0	0	1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
    ]), ...
    sparse([
        0	0	0	0	0	0	0	0	0;
        3	0	3	4	1	0	0	0	0;
        10	7	7	3	3	0	0	0	0;
        2	5	1	2	2	1	0	0	0;
        0	0	2	5	1	1	0	0	0;
        1	1	3	1	0	0	0	0	0;
        0	0	1	2	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0;
        0	0	0	1	1	0	0	0	0;
        0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	1;
    ])};

% remove single cells from consideration --- don't know if we're inducing
% post-mitotic cells
for i = 1:numel(times)
    data{i}(1+1,0+1) = 0;
end

    function [rp, r] = rfun
        r = random('beta', 2, 2);
        rp = betapdf(r, 2, 2);
        r = r / 2; % actually between 0 and 1/2
    end

    function [lp, l] = lambdafun
        l = random('logn', log(1/4), log(2));
        lp = lognpdf(l, log(1/4), log(2));
    end

    function [gp, g] = gammafun
        g = random('logn', 0, log(2));
        gp = lognpdf(g, 0, log(2));
    end

samples = sample3_no_shedding(@rfun, @gammafun, @lambdafun, [], {}, times, data, 1000, 30);

save cedric_k14_samples.mat samples

end