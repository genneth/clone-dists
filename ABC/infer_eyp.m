function infer_eyp(rhos, rs, lambdas, ts, datas, filename)

if numel(ts) ~= numel(datas)
    error 'need the same number of elements in ts and datas'
end

sz = size(datas{1});
if sz(2) ~= 3
    error 'datas does not have the right dimensions'
end

wh = waitbar(0, 'setting up...');

n = numel(ts);
for i = [1:n]
    ms(i) = max(datas{i}(:,1) + datas{i}(:,2));
end

k = max(ms);
fprintf(1, 'tracking a maximal clone size of %d\n', k);
[Tl Tr Tg] = generate_transition_matrix(k);
P0 = initial_eyp(k);

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
        for i = 1:numel(rhos)
            for j = 1:numel(rs)
                for h = 1:numel(lambdas)
                    rho = rhos(i);
                    r = rs(j);
                    lambda = lambdas(h);
                    pxd(j,i,h) =  ...
                        PDX(px(j,i,h) * scale, condPbs(Pbs(...
                            population(Tl, Tr, Tg, ...
                                lambda, r, lambda * rho / (1-rho), t, P0)...
                        )), data');
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
    end

pxd = prior;
timer_start = tic;
for i = [1:n]
    pxd = PXD(pxd, ts(i), datas{i});
end
pxd = pxd ./ integrate(pxd);
waitbar(1.0, wh, 'plotting and saving...');

% extract means and variances
fprintf(1, 'estimates:\n');

rho1 = integrate(repmat(rhos, [numel(rs), 1, numel(lambdas)]) .* pxd);
rho2 = integrate(repmat(rhos, [numel(rs), 1, numel(lambdas)]).^2 .* pxd);
fprintf(1, 'rho mean: %g; std: %g\n', rho1, sqrt(rho2 - rho1^2));

r1 = integrate(repmat(rs', [1, numel(rhos), numel(lambdas)]) .* pxd);
r2 = integrate(repmat(rs', [1, numel(rhos), numel(lambdas)]).^2 .* pxd);
fprintf(1, 'r mean: %g; std: %g\n', r1, sqrt(r2 - r1^2));

if numel(lambdas) > 1
    ls = reshape(lambdas, [1 1 numel(lambdas)]);
    l1 = integrate(repmat(ls, numel(rs), numel(rhos)) .* pxd);
    l2 = integrate(repmat(ls, numel(rs), numel(rhos)).^2 .* pxd);
    fprintf(1, 'lambda mean: %g; std: %g\n', l1, sqrt(l2 - l1^2));
end

% plot the resulting distribution
gf = newplot(figure);
if numel(lambdas) == 1
    save(filename, 'rhos', 'rs', 'pxd');
    contour3(gf, rhos, rs, pxd, 20);
    surface(rhos, rs, pxd, 'EdgeColor', [.8 .8 .8], 'FaceColor', 'none');
else
    save(filename, 'rhos', 'rs', 'lambdas', 'pxd');
    contourslice(gf, rhos, rs, lambdas, pxd, rho1,r1,l1);
    zlabel(gf, '$\lambda$', 'Interpreter', 'latex');
end
xlabel(gf, '$\rho$', 'Interpreter', 'latex');
ylabel(gf, '$r$', 'Interpreter', 'latex');
saveas(gf, filename, 'fig');

close(wh);

end
