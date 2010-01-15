function infer_eyp(rhos, rs, lambdas, ts, datas, filename)

if numel(ts) ~= numel(datas)
    error 'need the same number of elements in ts and datas'
end

sz = size(datas{1});
if sz(2) ~= 2
    error 'datas does not have the right dimensions'
end

wh = waitbar(0, 'setting up...');

n = numel(ts);
for i = [1:n]
    ms(i) = max(datas{i}(:,1));
end

% make an estimate of the maximum number of cells we need to track. Clearly
% we need at least as many as the data says. But we also need m_max >
% lambda t. 
m_max = max(ms);

fprintf(1, 'tracking a maximal clone size of %d\n', m_max);

    % integrate p as a distribution over rho-r-lambda space
    function y = integrate(p)
        if numel(lambdas) == 1
            y = trapz(rhos, trapz(rs, p));
        else
            y = trapz(lambdas, trapz(rhos, trapz(rs, p)));
        end
    end

prior = ones(numel(rs), numel(rhos), numel(lambdas)); % notice that we swap the first two indices
prior = prior ./ integrate(prior); % normalise!

waitbar(0, wh, 'inferring...');
total_evals = n * numel(rhos) * numel(rs) * numel(lambdas);
curr_evals = 0;
timer_start = 0;

    % need to control for underflow; before each round, scale things by a
    % suitable factor (staying within overflow) and renormalise afterwards
    function pxd = PXD(px, t, data)
        pxd = zeros(size(px));
        scale = 2^1022 / max(max(max(px))); % 1022 instead of 1023 because we're not greedy
        for h = 1:numel(lambdas)
            for j = 1:numel(rs)
                parfor i = 1:numel(rhos)
                    rho = rhos(i);
                    r = rs(j);
                    lambda = lambdas(h);
                    pxd(j,i,h) = ...
                        PDX(px(j,i,h) * scale, condPb(...
                            exact_pops(lambda, r, lambda * rho / (1-rho), t, m_max+1)...
                        ), data');
                    % output some progress
                    curr_evals = curr_evals+1;
                end
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
for i = [1:n]
    pxd = PXD(pxd, ts(i), datas{i});
end
pxd = pxd ./ integrate(pxd);
waitbar(1.0, wh, 'plotting and saving...');

% plot the resulting distribution
gf = newplot;
if numel(lambdas) == 1
    save(filename, 'rhos', 'rs', 'pxd');
    contour3(gf, rhos, rs, pxd, 20);
    surface(rhos, rs, pxd, 'EdgeColor', [.8 .8 .8], 'FaceColor', 'none');
else
    save(filename, 'rhos', 'rs', 'lambdas', 'pxd');
    contourslice(gf, rhos, rs, lambdas, pxd, [],[],lambdas);
    zlabel(gf, '$\lambda$', 'Interpreter', 'latex');
end
xlabel(gf, '$\rho$', 'Interpreter', 'latex');
ylabel(gf, '$r$', 'Interpreter', 'latex');
saveas(gf, filename, 'fig');

close(wh);

end
