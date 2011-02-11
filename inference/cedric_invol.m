function cedric_invol

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

% remove single cells from consideration --- don't know if we're inducing
% post-mitotic cells
for i = 1:numel(times)
    data{i}(1+1,0+1) = 0;
end

    function [rp, r] = rfun
        r = random('beta', 1, 3);
        rp = betapdf(r, 1, 3);
        r = r / 2; % actually between 0 and 1/2
    end

    function [lp, l] = lambdafun
        l = random('logn', 0, log(2));
        lp = lognpdf(l, 0, log(2));
    end

    function [gp, g] = gammafun
        g = random('logn', log(2), log(1.8));
        gp = lognpdf(g, log(2), log(1.8));
    end

samples = sample3_no_shed(@rfun, @gammafun, @lambdafun, [], {}, times, data, 1000, 30);

save cedric_invol_samples.mat samples

end