function oesophagus

ts3(1) = 3.0+1/7; % 22 days
data3{1} = sparse([
	 0  0 39 14  7  4  1  0;
	 0 59 24  6  3  0  0  0;
	32 33 17 15  1  1  0  0;
	 5 11 10  2  3  1  0  0;
	 4  6  3  1  1  0  1  0;
	 5  1  1  2  1  0  1  0;
	 1  4  3  2  0  1  0  1;
	 0  1  0  0  0  0  1  0;
	 0  0  0  0  1  0  0  0;
	 0  0  1  0  1  0  0  0;
	 0  0  0  1  0  0  0  0;
	 0  0  0  0  1  0  0  0
	]);

ts3(2) = 3.0+9/7; % 30 days
data3{2} = sparse([
	 0  0 26  9  3  2  1  0  0  0;
	 0 35 32  8  9  1  1  0  0  0;
	12 26 23 12  8  2  2  1  0  0;
	 3  8  6 10  7  2  0  0  1  0;
	 3  5  9  6  7  3  0  0  0  1;
	 1  0  2  9  6  2  1  0  0  0;
	 0  2  1  1  0  2  1  0  0  0;
	 0  1  0  1  1  2  1  3  1  0;
	 0  0  1  1  1  2  1  1  0  1;
	 0  0  0  0  2  1  0  1  0  0
	]);
data3{2}(16, 7) = 1;
data3{2}(17,10) = 1;

% take out unreliable entries
for i=1:numel(ts3)
    data3{i}(0+1,:)   = 0;
    data3{i}(1+1,0+1) = 0;
end

    function [rp, r] = rfun
        rp = 1;
        r = random('unif', 0, 1);
    end

    function [lp, l] = lambdafun
        l = random('logn', log(0.77), 0.5);
        lp = lognpdf(l, log(0.77), 0.5);
    end

    function [gp, g] = gammafun
        g = random('logn', 0, log(2));
        gp = lognpdf(g, 0, log(2));
    end

samples = sample3(@rfun, @gammafun, 0.0, @lambdafun, [], {}, ts3, data3, 1);

save('oesophagus_samples.mat', 'samples');

end
