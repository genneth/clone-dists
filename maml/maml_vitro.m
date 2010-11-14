function maml_vitro

ts2 = [3 7 10];
data2{1} = sparse([0 25 111 62 82 6 11 2 1]);
data2{2} = sparse([0 1 16 33 38 39 30 30 23 23 23 14 12 6 3 5 1 1 1 0 ...
    0 0 1]);
data2{3} = sparse([0 1 6 10 13 18 23 21 23 20 15 15 15 12 9 11 4 4 6 4 ...
    8 1 0 1 6 1 1 0 0 0 0 1 0 0 0 0 0 0 0 0 1]);

    function [rp, r] = rfun
        r = random('beta', 1, 5);
        rp = betapdf(r, 1, 5);
        r = r / 2; % actually between 0 and 1/2
    end

    function [lp, l] = lambdafun
        l = random('logn', log(1/3), log(2));
        lp = lognpdf(l, log(1/3), log(2));
    end

    function [gp, g] = gammafun
        g = random('logn', log(0.3), log(2));
        gp = lognpdf(g, log(0.3), log(2));
    end

rfunh = @rfun;
gfunh = @gammafun;
lfunh = @lambdafun;

samples = sample2(rfunh, gfunh, lfunh, ts2, data2, 10, 3);

save maml_vitro_samples.mat samples;

end
