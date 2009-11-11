function infer_fake_eyp

% generated ear data

% 50 samples (including ones to be thrown away)
t1 = 6.0;
data1 = [0 2 5; 0 3 3; 1 1 8; 1 2 2; 1 3 1; 2 0 5; 2 1 2; 2 2 2; 3 0 3; 3 1 1; 3 2 1; 4 0 1; 4 1 2; 4 2 1];

% several 400 sample runs
t2 = 6.0;
data2a = [0 2 33; 0 3 6; 1 1 65; 1 2 22; 1 3 5; 1 5 1; 2 0 59; 2 1 27; 2 2 11; 2 3 2; 3 0 17; 3 1 17; 3 2 7; 3 3 3; 3 4 2; 3 5 1; 4 0 10; 4 1 6; 4 2 3; 4 3 1; 4 4 2; 5 0 5; 5 1 1; 6 0 1; 7 0 1; 7 2 1; 8 2 1; 8 4 1];
data2b = [0 2 41; 0 3 1; 1 1 76; 1 2 18; 1 3 4; 1 5 1; 1 6 1; 2 0 56; 2 1 35; 2 2 10; 2 3 5; 3 0 25; 3 1 11; 3 3 4; 3 4 2; 4 0 4; 4 1 3; 4 2 3; 4 3 2; 5 0 4; 5 1 3; 6 0 2; 6 1 2; 7 2 1; 8 0 1; 9 0 1; 10 2 1];
data2c = [0 2 36; 0 3 7; 0 4 1; 0 5 1; 1 1 71; 1 2 19; 1 3 6; 2 0 58; 2 1 24; 2 2 15; 2 3 4; 2 4 1; 3 0 25; 3 1 16; 3 2 5; 3 3 4; 4 0 8; 4 1 8; 4 3 1; 5 0 4; 5 1 3; 5 2 1; 6 0 1; 6 2 1];
data2d = [0 2 39; 0 3 1; 0 5 1; 1 1 56; 1 2 24; 1 3 5; 1 5 1; 2 0 44; 2 1 37; 2 2 9; 2 3 3; 3 0 15; 3 1 19; 3 2 7; 3 3 2; 3 4 1; 4 0 10; 4 1 6; 4 2 3; 5 0 1; 5 2 1; 6 0 2; 6 5 1; 7 0 1; 7 1 2];
data2e = [0 2 39; 0 3 7; 0 4 1; 1 1 68; 1 2 18; 1 3 6; 1 4 1; 2 0 54; 2 1 26; 2 2 15; 2 4 2; 3 0 21; 3 1 11; 3 2 1; 3 3 2; 4 0 9; 4 1 4; 4 2 2; 5 0 3; 5 1 3; 5 2 1; 5 3 1; 6 0 1; 6 1 2];
data2f = [0 2 41; 0 3 12; 1 1 58; 1 2 22; 1 3 4; 2 0 59; 2 1 29; 2 2 9; 2 3 6; 2 4 1; 2 5 1; 3 0 17; 3 1 6; 3 2 3; 3 3 2; 4 0 6; 4 1 8; 4 2 4; 4 3 1; 4 4 1; 5 0 3; 5 1 4; 5 2 1; 5 3 1; 6 0 2; 6 2 1; 6 4 1; 7 0 1; 9 3 1];
data2g = [0 2 35; 0 3 5; 1 1 65; 1 2 12; 1 3 4; 1 4 1; 2 0 63; 2 1 39; 2 2 14; 2 3 2; 3 0 17; 3 1 14; 3 2 8; 3 3 4; 3 5 1; 4 0 6; 4 1 6; 4 2 2; 5 0 3; 5 2 2; 5 4 1; 6 0 1; 6 1 1; 7 0 1; 7 3 1];
data2h = [0 2 51; 0 3 4; 0 5 1; 1 1 71; 1 2 17; 1 4 2; 1 5 1; 2 0 65; 2 1 28; 2 2 12; 2 3 3; 3 0 20; 3 1 8; 3 2 7; 3 3 2; 3 5 1; 4 0 7; 4 1 2; 4 2 5; 4 3 2; 5 0 1; 5 1 3; 5 2 1; 5 3 1; 6 0 2; 6 1 1; 6 2 1; 6 4 1];

rhos = linspace(0.5, 0.7, 30);
rs = linspace(0.15, 0.35, 30);
lambda = 0.25;
t = t2;
datas = {data2a data2b data2c data2d data2e data2f data2g data2h};

wh = waitbar(0, 'setting up...');

n = numel(datas);
for i = [1:n]
    ms(i) = max(datas{i}(:,1) + datas{i}(:,2));
end

k = max(ms);
fprintf(1, 'tracking a maximal clone size of %d\n', k);
[Tl Tr Tg] = generate_transition_matrix(k);
P0 = initial_eyp(k);

    % integrate p as a distribution over rho space
    function y = integrate(p)
        y = trapz(rhos, trapz(rs, p));
    end

prior = ones(numel(rs), numel(rhos), numel(datas)); % notice that we swap the first two indices

waitbar(0, wh, 'inferring...');
total_evals = numel(rhos) * numel(rs);
curr_evals = 0;
timer_start = 0;

    % need to control for underflow; before each round, scale things by a
    % suitable factor (staying within overflow) and renormalise afterwards
    function pxd = PXD(px, t, data)
        pxd = zeros(size(px));
        scale = 2^1022 / max(max(max(px))); % 1022 instead of 1023 because we're not greedy
        for i = 1:numel(rhos)
            for j = 1:numel(rs)
                rho = rhos(i);
                r = rs(j);
                pbs = condPbs(Pbs(...
                        population(Tl, Tr, Tg, ...
                            lambda, r, lambda * rho / (1-rho), t, P0)...
                    ));
                for h = 1:numel(datas)
                    pxd(j,i,h) =   PDX(px(j,i,h) * scale, pbs, datas{h}');
                end
                % output some progress
                curr_evals = curr_evals+1;
                timer_curr = toc(timer_start)/60;
                waitbar(curr_evals/total_evals, wh, ...
                    sprintf('inferring... %d/%d, %0.2fmin elapsed, %0.2fmin to go', ...
                    curr_evals, total_evals, timer_curr, ...
                    timer_curr/(curr_evals/total_evals) - timer_curr));
            end
        end
    end

pxd = prior;
timer_start = tic;
pxd = PXD(pxd, t, datas{i});
waitbar(1.0, wh, 'plotting and saving...');

% plot the resulting distribution
gh = newplot(figure);
hold(gh, 'on');
for i = [1:numel(datas)]
    contour(gh, rhos, rs, pxd(:,:,i), 1);
end
hold(gh, 'off');
xlabel(gh, '$\rho$', 'Interpreter', 'latex');
ylabel(gh, '$r$', 'Interpreter', 'latex');

saveas(gh, 'fake_eyp_bs_400', 'pdf');

close(wh);

end
