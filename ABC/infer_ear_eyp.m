function infer_ear_eyp

t = 6.0;
databs = [
    0	133	33	5	3;
    40	46	7	2	0;
    42	38	12	0	0;
    12	20	3	1	1;
    7	7	4	2	0;
    0	2	2	0	0;
    1	0	0	0	0
    ];
databs(0+1, :)   = 0;
databs(1+1, 0+1) = 0;
data = Pbs2lst(databs);

rhos = linspace(0.3, 0.7, 30);
rs = linspace(0.01, 0.35, 31);
lambdas = linspace(0.2, 0.4, 29);

infer_eyp3(rhos, rs, lambdas, [t], {data}, 'ear_eyp_bs_cond');

end
