function oes_scaling

t1 = 3/7;
data1 = [1 142; 2 71; 3 8; 4 3];

t2 = 10/7;
data2 = [1 111; 2 109; 3 35; 4 15; 5 4; 6 2];

t3 = 3;
data3 = [1 140; 2 79; 3 38; 4 25; 5 11; 6 4; 7 1; 8 1; 11 1];

t4 = 6;
data4 = [1 73; 2 60; 3 44; 4 26; 5 22; 6 10; 7 3; 8 4; 10 2; 11 4; 13 2; 14 1; 20 1; 28 1];

t5 = 12;
data5 = [1 87; 2 87; 3 50; 4 52; 5 26; 6 24; 7 14; 8 8; 9 9; 10 8; 11 6; 12 7; 13 4; 14 3; 15 3; 16 2; 17 2; 20 2; 22 1; 26 1; 29 1];

t6 = 26;
data6 = [1 40; 2 59; 3 43; 4 18; 5 27; 6 23; 7 16; 8 12; 9 14; 10 10; 11 11; 12 5; 13 8; 14 7; 15 6; 16 3; 17 4; 18 3; 19 3; 20 4; 21 1; 22 2; 23 1; 24 3; 25 4; 27 1; 30 1; 33 1; 34 1; 35 5; 36 1; 37 2; 39 1; 50 2; 56 1; 85 1; 90 1];

t7 = 52;
data7 = [1 12; 2 18; 3 15; 4 14; 5 6; 6 11; 7 8; 8 6; 9 9; 10 9; 11 3; 12 4; 13 4; 14 2; 15 4; 16 4; 17 2; 18 4; 19 2; 20 2; 21 3; 22 4; 23 1; 24 5; 26 2; 27 4; 28 3; 29 1; 30 3; 31 3; 32 3; 33 2; 34 4; 35 2; 36 1; 37 1; 38 1; 40 1; 41 1; 42 1; 43 1; 44 1; 45 1; 47 2; 48 1; 49 1; 54 1; 55 3; 58 1; 60 2; 62 1; 63 1; 70 2; 75 1; 80 1; 85 1; 96 1; 110 1; 130 1; 140 1; 175 1];

ts = [t1 t2 t3 t4 t5 t6 t7];
data = {data1 data2 data3 data4 data5 data6 data7};
av = zeros(size(ts));
tot = zeros(size(ts));
for i = 1:numel(ts)
    tot(i) = sum(data{i}(:,2));
    av(i) = dot(data{i}(:,1), data{i}(:,2) / tot(i));
end

lambda = 0.7676;
r = 0.1926;
rho = 0.5164;
gamma = lambda * rho / (1-rho);
tau = rho / (r * lambda);

fh = figure;
gh = newplot(fh);
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
set(gh, 'NextPlot', 'add');
for i = 1:numel(ts);
    ps = condPb2(abs(exact_pops(lambda, r, lambda * rho / (1-rho), ts(i), floor(10*ts(i))+1)));
    ps = ps(2:end);
    tav = dot(1:numel(ps), ps);
    err = sqrt(ps .* (1-ps) ./ tot(i));
    plot(gh, (1:numel(ps)) ./ tav, log10(ps * tav) + (i-1), '-', 'Color', colours(i,:), 'LineWidth', 1);
    fill_plot(gh, (1:numel(ps)) ./ tav, ...
        log10(max(ps - err, 1e-10) * tav) + (i-1), ...
        log10(max(ps + err, 1e-10) * tav) + (i-1), ...
        colours(i,:), 0.4);

    sdata = sparse(data{i}(:,1), ones(size(data{i}(:,1))), data{i}(:,2), 1000000, 1);
    marker = 1;
    acc = 0;
    lo = zeros(0,0);
    hi = zeros(0,0);
    binned = zeros(0,0);
    for j = 1:numel(ps)
       acc = acc + ps(j);
       if (1-acc) ^ tot(i) < 0.01
           lo(end+1) = marker;
           hi(end+1) = j;
           binned(end+1) = sum(sdata(marker:j));
           marker = j+1;
           acc = 0;
       end
    end

    lh = plot(gh, (lo + hi) / 2 / av(i), log10((binned ./ (hi - lo + 1)) / tot(i) * av(i)) + (i-1), ...
        '^', 'MarkerEdgeColor', colours(i,:), 'MarkerFaceColor', colours(i,:), 'MarkerSize', 4);
    ah = get(lh, 'Annotation'); leh = get(ah, 'LegendInformation'); set(leh, 'IconDisplayStyle', 'off');
    for j = 1:numel(binned)
        lh = plot(gh, [lo(j)-0.5 hi(j)+0.5] / av(i), log10([binned(j) binned(j)] ./ (hi(j) - lo(j) + 1) / tot(i) * av(i)) + (i-1), ...
           '-', 'Color', colours(i,:));
        ah = get(lh, 'Annotation'); leh = get(ah, 'LegendInformation'); set(leh, 'IconDisplayStyle', 'off');
    end


    set(gh, 'XLim', [0 5], 'YLim', [-2 6]);
    pause(0.1);
end
plot(gh, linspace(0,5,50), -linspace(0,5,50).*log10(exp(1)) + 6, '-', 'Color', 'black', 'LineWidth', 2);

xlabel(gh, 'scaled clone size', 'FontSize', 8, 'FontName', 'Times');
set(gh, 'YTick', log10([reshape([1:9]' * 10.^[-2:5], 8*9, 1); 10^6]));
set(gh, 'YTickLabel', '10 ||||||||');
set(gh, 'TickLength', get(gh, 'TickLength') / 2);

set(gh, 'FontName', 'Times', 'FontSize', 8);

legend('3 days', ...
    '10 days', ...
    '3 weeks', ...
    '6 weeks', ...
    '3 months', ...
    '6 months', ...
    '1 year', ...
    'limit', ...
    'Location', 'SouthEast');

set(gh, 'Position', get(gh, 'OuterPosition') - ...
    get(gh, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);

set(fh, 'PaperUnits', 'inches');
w = 6.25; h = 7.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

set(fh, 'Color', 'white');

set(fh, 'Renderer', 'Painters');
pause(0.1);
set(fh, 'Renderer', 'OpenGL');
pause(0.1);
print(fh, '-dpng', '-opengl', '-r300', 'oes-scaling-b');

end
