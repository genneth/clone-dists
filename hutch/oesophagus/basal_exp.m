function basal_exp

oesophagus_data;

for i=1:numel(ts2)
    basal{i}(1+1) = 0;
    tot(i) = sum(basal{i});
    av(i) = dot(0:(numel(basal{i})-1), basal{i}/tot(i));
end

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

fh = figure;
set(fh, 'PaperUnits', 'inches');
w = 7; h = 4.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

ah = newplot(fh);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 10);
xlabel(ah, 'scaled clone size', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
ylabel(ah, 'scaled frequency', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');

for i = 2:numel(ts2)
    plot((2:(numel(basal{i})-1)) ./ av(i), basal{i}((2+1):end) * av(i) / tot(i), 's', ...
        'Color', colours(i-1,:));
end
plot(linspace(0,10,100), exp(-linspace(0,10,100)), 'k-', 'LineWidth', 1);

set(ah, 'XLim', [0 5], 'YLim', [1e-2 1e1]);
set(ah, 'YScale', 'log');

legend(ah, {'10 days', '3 weeks', '6 weeks', '3 months', '6 months', '1 year', 'exponential'});

print(fh, '-dpng', '-painters', '-r300', 'basal-dist');

close(fh);


end