function infer_oes_eyp

%t1 = 3.0+1/7; % 22 days
%data1 = [2 0 16; 1 1 18; 0 2 20; 3 0 4; 2 1 19; 1 2 10; 0 3 7; 4 0 3; 3 1 8; 2 2 5; 1 3 1; 0 4 4; 4 1 5; 3 2 10; 2 3 3; 1 4 1; 0 5 1; 6 0 1; 5 1 2; 4 2 4; 3 3 3; 2 4 1; 1 5 2; 6 1 1; 5 2 2; 3 4 1; 6 2 3; 5 3 1; 4 4 4; 6 3 1; 5 4 1; 4 5 2; 6 4 1; 9 2 1; 8 3 2; 10 3 1; 9 4 1; 8 5 1; 10 4 1; 11 4 1; 12 4 1; 16 6 1; 18 5 1];

%t2 = 10/7; % 10 days
%data2 = [2 0 6; 1 1 26; 3 0 1; 2 1 15; 1 2 12; 4 0 1; 3 1 3; 2 2 9; 1 3 2; 4 1 1; 3 2 1; 2 3 3; 1 4 1; 6 0 1; 5 3 1; 6 8 1];

% data from ATRA control
t1 = 3.0+1/7; % 22 days
data1bs = [
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
	];
data1 = Pbs2lst(data1bs);

t2 = 3.0+9/7; % 30 days
data2bs = [
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
	];
data2 = Pbs2lst(data2bs);
data2(end+1,:) = [16  7 1];
data2(end+1,:) = [17 10 1];

rhos = linspace(0.3, 0.65, 30);
rs = linspace(0.1, 0.3, 31);
lambdas = linspace(0.6, 1.1, 32); % around 0.77/week

infer_eyp(rhos, rs, lambdas, [t1], {data1}, 'oes_eyp_bs_day_22_atra_control');
infer_eyp(rhos, rs, lambdas, [t2], {data2}, 'oes_eyp_bs_day_30_atra_control');

end
