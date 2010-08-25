function ear

ts3(1) = 6; % 6 weeks
data3{1} = sparse([
    0	133	33	5	3;
    40	46	7	2	0;
    42	38	12	0	0;
    12	20	3	1	1;
    7	7	4	2	0;
    0	2	2	0	0;
    1	0	0	0	0
	]);

% take out unreliable entries
for i=1:numel(ts3)
    data3{i}(0+1,:)   = 0;
    data3{i}(1+1,0+1) = 0;
end

    function [rp, r] = rfun
        r = random('beta', 4, 4);
        rp = betapdf(r, 4, 4);
    end

    function [lp, l] = lambdafun
        l = random('logn', log(1/4), log(2));
        lp = lognpdf(l, log(1/4), log(2));
    end

    function [gp, g] = gammafun
        g = random('logn', 0, log(2));
        gp = lognpdf(g, 0, log(2));
    end

samples = sample3(@rfun, @gammafun, 0.0, @lambdafun, [], {}, ts3, data3, 900);

save('ear_samples.mat', 'samples');

end
