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

ts2 = [3/7 3 6 13 26 52];
data2{1} = sparse([0 74 5]);
data2{2} = sparse([0 65 65 19 8 1 1]);
data2{3} = sparse([0 118 112 53 27 17 6 1 3]);
data2{4} = sparse([0 81 56 35 15 14 4 7 5 2 2 1 2 0 0 1]);
data2{5} = sparse([0 50 57 38 20 24 15 8 2 8 6 2 2 4 0 0 1 0 0 0 1 0 1 0]);
data2{6} = sparse([0 25 33 26 21 15 10 9 6 8 6 1 5 6 2 3 7 2 0 1 1 1 0 ...
    0 0 2 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1]);

% take out unreliable entries
for i = 1:numel(ts2)
    data2{i}(0+1) = 0;
    data2{i}(1+1) = 0;
end

samples = sample3_no_shed(...
    @()(random('beta', 2, 2) / 2), ...
    @()(random('logn', log(1/4), log(2))), ...
    @()(random('logn', 0, log(2))), ...
    ts2, data2, ts3, data3, 1000, 30);

save ear_samples.mat samples;

end
