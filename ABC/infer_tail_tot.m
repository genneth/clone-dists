function infer_tail_tot

% todo: data

rhos = linspace(0.3, 0.65, 30);
rs = linspace(0.1, 0.3, 31);
lambdas = linspace(0.6, 1.1, 32); % around 0.77/week

infer_eyp_tot(rhos, rs, lambdas, [t1], {data1}, 'tail_tot');

end
