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

t4 = 6;
data4bs = [
     0  0 12  2  4  2  0  1  0  0  0  0  0  0  0;
     0 15 10  7  9  1  0  0  0  0  0  0  0  0  0;
     4 26 12  9  4  1  1  0  1  0  0  0  0  0  0;
     4  4 17 10  7  2  2  1  0  0  0  0  1  0  0;
     1  5  5  8  6  5  3  0  1  0  0  0  1  0  0;
     0  1  3  4  2  6  5  1  1  0  1  0  0  0  0;
     0  0  2  1  2  3  2  3  0  1  0  1  0  0  0;
     0  0  1  3  3  1  0  0  2  0  1  0  0  0  0;
     0  0  0  1  0  0  2  0  1  0  0  1  0  0  0;
     0  0  0  0  0  0  5  1  0  1  0  0  1  0  0;
     0  0  0  0  2  0  0  1  0  3  0  0  1  1  0;
     0  0  0  0  0  0  3  1  1  0  0  0  0  0  0;
     0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
     0  0  0  0  1  0  0  0  0  0  0  0  0  0  0;
     0  0  0  0  0  0  0  0  0  0  0  0  1  0  1;
     0  0  0  0  0  0  1  0  1  0  0  0  0  0  0
    ];
data4bs(0+1,:) = 0;
data4 = Pbs2lst(data4bs);
data4(end+1,:) = [20 20 1];
data4(end+1,:) = [21 13 1];
data4(end+1,:) = [25 21 1];
data4(end+1,:) = [46 40 1];
data4(end+1,:) = [51 26 1];

infer_eyp3(linspace(0.6,0.8,30),   linspace(0,0.25,31),   linspace(2,3,32),      [t1], {data1}, 'oes_eyp_bs_day_4_atra');
infer_eyp3(linspace(0.3,0.5,30),   linspace(0.05,0.2,31), linspace(2.2,2.7,32),  [t2], {data2}, 'oes_eyp_bs_day_11_atra');
infer_eyp3(linspace(0.35,0.6,30),  linspace(0.1,0.3,31),  linspace(1.2,1.6,32),  [t3], {data3}, 'oes_eyp_bs_week_3_atra');
infer_eyp3(linspace(0.2,0.4,30),   linspace(0.1,0.25,31), linspace(0.8,1.1,32),  [t4], {data4}, 'oes_eyp_bs_week_6_atra');

end
