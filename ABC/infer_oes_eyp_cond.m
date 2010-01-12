function infer_oes_eyp_cond

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
data1bs(0+1,:) = 0; % condition out suprabasal only
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
data2bs(0+1,:) = 0; % condition out suprabasal only clones
data2 = Pbs2lst(data2bs);
data2(end+1,:) = [16  7 1];
data2(end+1,:) = [17 10 1];

rhos = linspace(0.3, 0.65, 30);
rs = linspace(0.1, 0.35, 31);
lambdas = linspace(0.6, 1.0, 32); % around 0.77/week

infer_eyp3(rhos, rs, lambdas, [t1], {data1}, 'oes_eyp_bs_day_22_atra_control_cond');
infer_eyp3(rhos, rs, lambdas, [t2], {data2}, 'oes_eyp_bs_day_30_atra_control_cond');

end
