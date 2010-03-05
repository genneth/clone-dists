function infer_oes_atra

% atar "pretreat", i.e. atra homeostatis

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
data3b = sum(data3bs,2);
data3 = [(0:numel(data3b)-1)' data3b];
data3(end+1,:) = [15 1];
data3(end+1,:) = [17 1];
data3(end+1,:) = [19 2];
data3 = data3(3:end,:); % remove 0 and 1

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
data4b = sum(data4bs,2);
data4 = [(0:numel(data4b)-1)' data4b];
data4(end+1,:) = [20 1];
data4(end+1,:) = [21 1];
data4(end+1,:) = [25 1];
data4(end+1,:) = [46 1];
data4(end+1,:) = [51 1];
data4 = data4(3:end,:); % remove the 0 and 1

infer_eyp(linspace(0.1,0.9,30), linspace(0.1,0.4,31), linspace(0.5,1.5,32),  [t3], {data3}, 'oes_eyp_b_week_3_atra');
infer_eyp(linspace(0.1,0.9,30), linspace(0.1,0.4,31), linspace(0.5,1.5,32),  [t4], {data4}, 'oes_eyp_b_week_6_atra');

end
