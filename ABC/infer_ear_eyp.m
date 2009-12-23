function infer_ear_eyp

t1 = 3.0;
data1 = [2 0 11; 1 1 4; 0 2 6; 3 0 1; 2 1 3; 0 3 1; 4 0 1; 4 1 1; 6 1 1];

t2 = 6.0;
data2 = [2 0 42; 3 0 12; 4 0 7; 6 0 1; 1 1 46; 2 1 38; 3 1 20; 4 1 7; 5 1 2; 0 2 33; 1 2 7; 2 2 12; 3 2 3; 4 2 4; 5 2 2; 0 3 5; 1 3 2; 3 3 1; 4 3 2; 0 4 3; 3 4 1];

rhos = linspace(0.5, 0.7, 30);
rs = linspace(0.15, 0.35, 31);
lambdas = linspace(0.2, 0.3, 29);

infer_eyp(rhos, rs, lambdas, [t1 t2], {data1 data2}, 'ear_eyp_bs');

end