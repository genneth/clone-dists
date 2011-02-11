function oesophagus_atra_comparison_bubble_plot

% control
ts3(1) = 22/7;
data3{1} = sparse([
	 0  0 39 14  7  4  1  0;
	 0 59 24  6  3  0  0  0;
	32 33 17 15  1  1  0  0;
	 5 11 10  2  3  1  0  0;
	 4  6  3  1  1  0  1  0;
	 5  1  1  2  1  0  1  0;
	 1  4  3  2  0  1  0  1;
	 0  1  0  0  0  0  1  0;
	 0  0  0  0  1  0  0  0;
	 0  0  1  0  1  0  0  0;
	 0  0  0  1  0  0  0  0;
	 0  0  0  0  1  0  0  0
	]);
data3{1}(0+1,:) = 0;
tot{1} = sum(sum(data3{1}));

% ATRA
ts3(2) = 3;
data3{2} = sparse([
	0	31	23	21	10	7	1	0	0	1	0	0	0;
	9	21	24	10	10	5	3	0	1	0	0	0	0;
	8	25	16	14	4	1	2	1	1	0	0	0	0;
	4	2	14	5	5	6	0	2	1	0	0	1	0;
	3	2	7	7	4	2	3	3	1	0	0	0	0;
	0	2	2	1	1	3	0	0	1	0	0	0	0;
	0	1	1	2	1	2	2	1	0	1	0	0	0;
	1	0	0	1	0	5	0	0	0	1	0	0	0;
	0	0	0	0	0	1	3	0	0	0	0	1	0;
	0	0	0	0	1	0	2	0	0	0	0	0	0;
	0	0	0	0	0	1	1	0	0	0	0	1	0;
	0	0	0	0	0	0	0	0	0	0	0	1	0;
	0	0	0	0	0	0	1	0	0	0	1	0	0]);
data3{2}(15+1,12+1) = 1;
data3{2}(17+1, 7+1) = 1;
data3{2}(19+1,10+1) = 1;
data3{2}(19+1,20+1) = 1;
data3{2}(0+1,:) = 0;
data3{2}(1+1,0+1) = 0;
tot{2} = sum(sum(data3{2}));

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

kk=6;

fh = figure;
ah = newplot;
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(ah, 'Basal', 'FontName', 'Helvetica', 'FontSize', 9, 'FontWeight', 'bold');
ylabel(ah, 'Suprabasal', 'FontName', 'Helvetica', 'FontSize', 9, 'FontWeight', 'bold');
set(ah, 'DataAspectRatio', [1 1 1]);
set(ah, 'XLim', [0.5 kk+0.5], 'YLim', [-0.5 kk+0.5]);
%set(ah, 'XAxisLocation', 'top', 'YDir', 'reverse', 'Box', 'on');
set(ah, 'Box', 'on');

for i=2:kk+1
    for j=1:kk+1
        if data3{2}(i,j) > 0
            r = sqrt(data3{2}(i,j)/tot{2});
            rectangle('Position', [i-1-r j-1-r 2*r 2*r], 'Curvature', [1 1], ...
                'FaceColor', colours(1,:) * 0.50 + [1 1 1] * 0.50, 'EdgeColor', 'none');
        end
        
        if data3{1}(i,j) > 0
            r = sqrt(data3{1}(i,j)/tot{2});
            rectangle('Position', [i-1-r j-1-r 2*r 2*r], 'Curvature', [1 1], ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1);
        end
        
    end
end

set(fh, 'PaperUnits', 'inches');
w = 4; h = 4;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'oesophaguse-atra-comparison-plot');

close(fh);

% make the legend

fh = figure;
ah = newplot(fh);
set(ah, 'NextPlot', 'add');
set(ah, 'FontName', 'Helvetica', 'FontSize', 9);
set(ah, 'DataAspectRatio', [1 1 1]);
set(ah, 'XLim', [-0.5 2.0], 'YLim', [-0.5 kk+0.5]);
set(ah, 'Visible', 'off');

labels = {'1%', '5%', '10%'};
r = sqrt([1 5 10] / 100);
y = 2 + cumsum(2*r) + [0 0.2 0.4];
for i = 1:numel(r)
    rectangle('Position', [0-r(i) y(i)-r(i) 2*r(i) 2*r(i)], 'Curvature', [1 1], ...
        'FaceColor', colours(1,:) * 0.50 + [1 1 1] * 0.50, 'EdgeColor', 'none');
    rectangle('Position', [1-r(i) y(i)-r(i) 2*r(i) 2*r(i)], 'Curvature', [1 1], ...
        'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1);
    text(1.1+r(i), y(i), labels{i}, ...
        'FontName', 'Helvetica', 'FontSize', 8, 'Color', 'k');
end

text(0, 1.8-sqrt(0.01), 'control', 'Rotation', 90, 'HorizontalAlignment', 'Right', ...
    'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'bold', 'Color', colours(1,:));
text(1, 1.8-sqrt(0.01), 'ATRA', 'Rotation', 90, 'HorizontalAlignment', 'Right', ...
    'FontName', 'Helvetica', 'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k');

set(fh, 'PaperUnits', 'inches');
w = 1; h = 4;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'oesophaguse-atra-comparison-plot-legend');

close(fh);

end
