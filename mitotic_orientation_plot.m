function mitotic_orientation_plot

freq = [6 5 6 18 17 20 30 21 27];
prop = freq ./ sum(freq);
los = deg2rad(0:10:80);
his = deg2rad(10:10:90);

fh = figure;
ah = newplot;

h = polar(ah, 0,0.25);
set(ah, 'NextPlot', 'add');
set(h, 'LineStyle', 'none');

set(ah, 'FontName', 'Helvetica', 'FontSize', 9);
xlabel(ah, 'Proportion', 'FontName', 'Helvetica', 'FontSize', 9, 'FontWeight', 'bold');
ylabel(ah, 'Proportion', 'FontName', 'Helvetica', 'FontSize', 9, 'FontWeight', 'bold');

for i=1:9
    % sector
    X = prop(i) * cos(los(i):0.01:his(i));
    Y = prop(i) * sin(los(i):0.01:his(i));
%     plot([0 X 0], [0 Y 0], '-', ...
%         'Color', [0.152941, 0.113725, 0.490196], ...
%         'LineWidth', 1);
    fill([0 X], [0 Y], ...
        [0.152941, 0.113725, 0.490196], ...
        'EdgeColor', [0 0 0]);
        

    % vertical bar
    X = (prop(i) + [0 sqrt(freq(i))] / sum(freq)) * cos((los(i)+his(i))/2);
    Y = (prop(i) + [0 sqrt(freq(i))] / sum(freq)) * sin((los(i)+his(i))/2);
    plot(X,Y, '-k');
    
    % upper tic
    X = (prop(i) + sqrt(freq(i)) / sum(freq)) * cos((los(i)+his(i))/2 + (-0.02:0.01:0.02));
    Y = (prop(i) + sqrt(freq(i)) / sum(freq)) * sin((los(i)+his(i))/2 + (-0.02:0.01:0.02));
    plot(X,Y, '-k');

%     % lower tic
%     X = (prop(i) - sqrt(freq(i)) / sum(freq)) * cos((los(i)+his(i))/2 + (-0.02:0.01:0.02));
%     Y = (prop(i) - sqrt(freq(i)) / sum(freq)) * sin((los(i)+his(i))/2 + (-0.02:0.01:0.02));
%     plot(X,Y, '-k');
end

axis(ah, [-0.001 inf -0.001 inf]);

set(fh, 'PaperUnits', 'inches');
w = 2.5; h = 2.5;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

print(fh, '-dpdf', '-painters', 'mitotic-orientation-plot');

close(fh);

end
