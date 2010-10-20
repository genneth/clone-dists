function clayton_plot

% The goal is to reproduce the theoretical curves for clone size
% distributions as found in Clayton et al. Nature 2007, Figure 2d. Only
% *survining* clones with more than 1 cell are counted.

nn = 5;
N = 2^nn;
lambda = 1.1;
r = 0.08;
rho = 0.22;
gamma = rho / (1-rho);

ts = logspace(0,2,40);

% first up, using the generating function solution for basal counts
    function bps = bin2(ps)
        bps = zeros(nn,numel(ts));
        ps = ps ./ repmat(1 - ps(0+1,:) - ps(1+1,:), N+1, 1);
        ps(0+1,:) = 0;
        ps(1+1,:) = 0;
        for i = 1:nn
            bps(i,:) = sum(ps((2^(i-1)+2):(2^i+1),:),1);
        end
    end

pops = bin2(clone_dist_b(r, gamma, lambda*ts, N));
fh = figure;
ah = newplot(fh);
semilogx(ah, ts, pops);
set(ah, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(ah, 'time (weeks)', 'FontName', 'Helvetica', 'FontSize', 9, 'FontWeight', 'bold');
ylabel(ah, 'proportion', 'FontName', 'Helvetica', 'FontSize', 9, 'FontWeight', 'bold');

legend(ah, '2', '3--4', '5--8', '9--16', '17--32');

set(fh, 'PaperUnits', 'inches');
w = 4; h = 3;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'clayton-plot-2');
close(fh);

% next up, using the basal+suprabasal exponentiation routine
    function bps = bin3(ps)
        ps = squeeze(sum(ps, 2));
        bps = bin2(ps(1:(N+1),:));
    end

[Tl Tr Tg P0] = clone_dist_bs_expv_setup(4*N);
pops = bin3(clone_dist_bs_expv(Tl,Tr,Tg,P0, r,gamma,lambda*ts));
fh = figure;
ah = newplot(fh);
semilogx(ah, ts, pops);
set(ah, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(ah, 'time (weeks)', 'FontName', 'Helvetica', 'FontSize', 9, 'FontWeight', 'bold');
ylabel(ah, 'proportion', 'FontName', 'Helvetica', 'FontSize', 9, 'FontWeight', 'bold');

legend(ah, '2', '3--4', '5--8', '9--16', '17--32');

set(fh, 'PaperUnits', 'inches');
w = 4; h = 3;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'clayton-plot-3');
close(fh);

end
