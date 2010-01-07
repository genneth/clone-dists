function infer_oes_eyp

t1 = 3.0+1/7; % 22 days
data1 = [2 0 16; 1 1 18; 0 2 20; 3 0 4; 2 1 19; 1 2 10; 0 3 7; 4 0 3; 3 1 8; 2 2 5; 1 3 1; 0 4 4; 4 1 5; 3 2 10; 2 3 3; 1 4 1; 0 5 1; 6 0 1; 5 1 2; 4 2 4; 3 3 3; 2 4 1; 1 5 2; 6 1 1; 5 2 2; 3 4 1; 6 2 3; 5 3 1; 4 4 4; 6 3 1; 5 4 1; 4 5 2; 6 4 1; 9 2 1; 8 3 2; 10 3 1; 9 4 1; 8 5 1; 10 4 1; 11 4 1; 12 4 1; 16 6 1; 18 5 1];

t2 = 10/7; % 10 days
data2 = [2 0 6; 1 1 26; 3 0 1; 2 1 15; 1 2 12; 4 0 1; 3 1 3; 2 2 9; 1 3 2; 4 1 1; 3 2 1; 2 3 3; 1 4 1; 6 0 1; 5 3 1; 6 8 1];

rhos = linspace(0.2, 0.6, 30);
rs = linspace(0.1, 0.3, 30);
lambdas = linspace(0.5, 1.0, 30); % around 1/week

infer_eyp(rhos, rs, lambdas, [t1 t2], {data1 data2}, 'infer_oes_eyp_bs');

end
