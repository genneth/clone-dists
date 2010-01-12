function infer_oes_atra

% atar "pretreat", i.e. atra homeostatis

t1 = 4/7;
data1bs = [
     0  0 57 13  0  0;
     0 73 30  8  0  1;
    26 54 14  4  0  0;
     8  6  3  1  0  0;
     4  3  1  0  0  0;
     0  1  1  0  0  0;
     0  0  0  1  0  1;
     0  0  0  0  0  0;
     0  0  0  0  0  0;
     0  0  1  0  0  0
	];
data1bs(0+1,:) = 0; % condition out suprabasal only
data1 = Pbs2lst(data1bs);

t2 = 11/7;
data2bs = [
     0  0 14  6  5  1  2  1;
     0 36 33 17  7  2  0  0;
    14 38 33 18  6  2  3  1;
     6 18 11 22  8  3  1  2;
     2 10 14 12  6  2  0  0;
     3  4  7  5  2  1  2  0;
     0  2  4  5  4  0  1  1;
     0  0  1  0  4  0  1  2;
     0  0  1  1  2  0  0  0;
     0  0  0  0  1  0  0  2;
     0  0  0  1  0  1  0  0;
     0  0  0  0  0  1  0  0;
     0  0  0  0  1  0  0  0;
     0  0  0  0  1  0  0  0
    ];
data2bs(0+1,:) = 0; % condition out suprabasal only
data2 = Pbs2lst(data2bs);
data2(end+1,:) = [ 7  9 1];
data2(end+1,:) = [10  9 1];
data2(end+1,:) = [17  7 1];

t3 = 3;
data3bs = [
     0  0 23 21 10  7  1  0  0  1  0  0;
     0 21 24 10 10  5  3  0  1  0  0  0;
     8 25 16 14  4  1  2  1  1  0  0  0;
     4  2 14  5  5  6  0  2  1  0  0  1;
     3  2  7  7  4  2  3  3  1  0  0  0;
     0  2  2  1  1  3  0  0  1  0  0  0;
     0  1  1  2  1  2  2  1  0  1  0  0;
     1  0  0  1  0  5  0  0  0  1  0  0;
     0  0  0  0  0  1  3  0  0  0  0  1;
     0  0  0  0  1  0  2  0  0  0  0  0;
     0  0  0  0  0  1  1  0  0  0  0  1;
     0  0  0  0  0  0  0  0  0  0  0  1;
     0  0  0  0  0  0  1  0  0  0  1  0
    ];
data3bs(0+1,:) = 0; % condition out suprabasal only
data3 = Pbs2lst(data3bs);
data3(end+1,:) = [15 12 1];
data3(end+1,:) = [17  7 1];
data3(end+1,:) = [19 10 1];
data3(end+1,:) = [19 20 1];

rhos = linspace(0.3, 0.6, 10);
rs = linspace(0.0, 0.3, 11);
lambdas = linspace(1, 1.6, 12);

infer_eyp3(rhos, rs, lambdas, [t1], {data1}, 'oes_eyp_bs_day_4_atra');
infer_eyp3(rhos, rs, lambdas, [t2], {data2}, 'oes_eyp_bs_day_11_atra');
infer_eyp3(rhos, rs, lambdas, [t3], {data3}, 'oes_eyp_bs_week_3_atra');

end
