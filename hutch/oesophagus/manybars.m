% Graph of comparison between theory and experiment, using lots of little
% bar charts. Fun times.
function manybars(t, data, r, gamma, mu, lambda, M, N, output)

% collate things first
databs = zeros(101,1001);
for i=0:100
    for j=0:100
        databs(i+1,j+1) = 0;
        for k=1:numel(data);
            [m,n] = size(data{k});
            if i+1<=m && j+1<=n
                databs(i+1,j+1) = databs(i+1,j+1) + data{k}(i+1,j+1);
            end
        end
    end
end
databs(0+1,:) = 0;
databs(1+1,0+1) = 0;
total = sum(sum(databs));

% theory = clone_dist_bs_shed(r, gamma, mu, t*lambda, 10,10);
% ignore = generating_function2(r, gamma, t*lambda, 0, 0) + theory(1+1,0+1);
% theory = theory ./ (1 - ignore);
% theory(0+1,:)   = 0;
% theory(1+1,0+1) = 0;
% save(output, 'theory');
theory = [];
load(output);

% load the colours
colours = [];
colour_palette;

fh = figure;
set(fh, 'PaperUnits', 'inches');
w = 3; h = 3;
set(fh, 'PaperSize', [w h]);
set(fh, 'PaperPosition', [0 0 w h]);

lim = 30;

margin_x = 0.1; margin_y = 0.1;
for i=1:M
    for j=0:N
        ah = subplot('Position', [
            margin_x+(j)/(N+1)*(1-2*margin_x)
            margin_y+(M-i)/M*(1-2*margin_y)
            (1-2*margin_x)/(N+1)
            (1-2*margin_y)/M]);
        set(ah, 'FontName', 'Helvetica', 'FontSize', 8);
        set(ah, 'Box', 'off');
        set(ah, 'NextPlot', 'add');
        bar(ah, 0, databs(i+1,j+1), ...
            'FaceColor', colours(1,:), ...
            'EdgeColor', colours(1,:), ...
            'LineWidth', 0.5);
        plot(ah, [-0.1 0.1], [total*theory(i+1,j+1) total*theory(i+1,j+1)], ...
            '-', 'LineWidth', 1.0, 'Color', colours(2,:));
        plot(ah, [0 0]', [total - binoinv(0.975, total, 1-theory(i+1,j+1))
            binoinv(0.975, total, theory(i+1,j+1))], ...
            '-', 'LineWidth', 1.0, 'Color', colours(2,:));
        set(ah, 'YLim', [0 lim]);
        set(ah, 'XAxisLocation', 'Top');
        set(ah, 'XTickLabel', {}, 'YTickLabel', {});
        set(ah, 'XTick', [], 'YTick', []);
        if i==1
            xlabel(ah, sprintf('%d', j), 'FontName', 'Helvetica', 'FontSize', 8);
        end
        if j==0
            ylabel(ah, sprintf('%d', i), 'FontName', 'Helvetica', 'FontSize', 8);
        end
        if i==1
            % 2nd axis just to get scale on the right
            ah2 = subplot('Position', get(ah, 'Position'));
            set(ah2, 'YAxisLocation', 'Right');
            set(ah2, 'YTick', [0 lim/2 lim]);
            set(ah2, 'YGrid', 'on', 'YMinorGrid', 'off');
            if j==N
                set(ah2, 'YTickLabel', {0 lim/2 lim});
            else
                set(ah2, 'YTickLabel', {});
            end
        end
        drawnow;
    end
end

print(fh, '-dpdf', '-painters', output);
close(fh);


end