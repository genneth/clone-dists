function bubble_compare(t, databs, r, gamma, file)

[rows,cols] = find(databs);
m = max(rows); n = max(cols);
k = m+n - 2;
kk = 6;
[Tl Tr Tg P0] = clone_dist_bs_expv_setup(k);
theory = clone_dist_bs_expv_cond(Tl,Tr,Tg,P0, r,gamma, t);

experiment = databs;
experiment(0+1,:) = 0;
experiment(1+1,0+1) = 0;
experiment = full(experiment);
count = sum(sum(experiment));

chances = binocompare(experiment(2:m,1:n), theory(2:m,1:n));
fh = fopen(strcat(file, '.tex'), 'w');
fprintf(fh, '~');
for i = 1:kk
    fprintf(fh, ' & %d', i);
end
fprintf(fh, ' \\\\\n');
fprintf(fh, '\\hline\n');
for j = 1:kk+1
    fprintf(fh, '%d', j-1);
    for i = 2:kk+1
        if chances(i-1,j) < 1/15787.192684
            fprintf(fh, ' & \\multicolumn{1}{>{\\color{white}\\columncolor{foursd}}c}{${%d\\,}^{%.1f}_{%.0f}$}', ...
                experiment(i,j), theory(i,j)*count, 1/chances(i-1,j));
        elseif chances(i-1,j) < 1/370.398347380
            fprintf(fh, ' & \\multicolumn{1}{>{\\color{white}\\columncolor{threesd}}c}{${%d\\,}^{%.1f}_{%.0f}$}', ...
                experiment(i,j), theory(i,j)*count, 1/chances(i-1,j));
        else
            fprintf(fh, ' & {${%d\\,}^{%.1f}_{%.0f}$}', ...
                experiment(i,j), theory(i,j)*count, 1/chances(i-1,j));
        end
    end
    fprintf(fh, ' \\\\\n');
end
fclose(fh);
fprintf(1, 'critical repetitions:\n');
disp(1 ./ chances);

fh = figure;
ah = newplot(fh);

colours = [
 0.996078, 0.360784, 0.027451;
% 0.996078, 0.988235, 0.0352941;
 0.541176, 0.713725, 0.027451;
 0.145098, 0.435294, 0.384314;
 0.00784314, 0.509804, 0.929412;
 0.152941, 0.113725, 0.490196;
 0.470588, 0.262745, 0.584314;
 0.890196, 0.0117647, 0.490196;
 0.905882, 0.027451, 0.129412
];

set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(ah, 'Basal', 'FontName', 'Helvetica', 'FontSize', 9);
ylabel(ah, 'Suprabasal', 'FontName', 'Helvetica', 'FontSize', 9);
set(ah, 'DataAspectRatio', [1 1 1]);
set(ah, 'XLim', [0.5 kk+0.5], 'YLim', [-0.5 kk+0.5]);
set(ah, 'XAxisLocation', 'top', 'YDir', 'reverse', 'Box', 'on');

for i=2:kk+1
    for j=1:kk+1
        % experiment
        if experiment(i,j) > 0
            r = sqrt(experiment(i,j)/count);
            rectangle('Position', [i-1-r j-1-r 2*r 2*r], 'Curvature', [1 1], ...
                'FaceColor', colours(1,:), 'EdgeColor', 'none');
        end
        
        % theory
        if theory(i,j) > 0
            r = sqrt(theory(i,j));
            rectangle('Position', [i-1-r j-1-r 2*r 2*r], 'Curvature', [1 1], ...
                'FaceColor', 'none', 'EdgeColor', colours(4,:), ...
                'LineWidth', 1.2);
        end
    end
end

set(fh, 'PaperUnits', 'inches');
w = 4; h = 4;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

if ~strcmp(file, '')
    print(fh, '-dpdf', '-painters', strcat(file, '.pdf'));
end

end
