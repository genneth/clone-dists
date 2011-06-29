function atRA_posttreat

atra_posttreat = {}; atra_posttreat_ts = {};
oesophagus_data;
colours = [];
colour_palette;

for i=1:numel(atra_posttreat_ts)

    data = sum(atra_posttreat{i}, 2);
    m = numel(data)-1;
    data(0+1) = 0;
    tot = sum(data);
    theory = inverse_z_transform(@(z) posttreat_basal_gf(3, atra_posttreat_ts(i)-3, z,z), m, 1e-6, 1e-6);
    theory = theory / (1-theory(0+1));
    theory(0+1) = 0;

    fh = figure;
    gh = newplot(fh);
    set(gh, 'NextPlot', 'add');
    set(gh, 'XLim', [0 15]);
    set(gh, 'FontName', 'Helvetica', 'FontSize', 9);
    xlabel(gh, 'no. basal cells', 'FontSize', 9, 'FontName', 'Helvetica', 'FontWeight', 'bold');
    set(gh, 'Box', 'off');
%     set(gh, 'YTick', 0:20:80, 'YTickLabel', {0 20 40 60 80});
    set(gh, 'YGrid', 'on', 'YMinorGrid', 'off');

    plot(gh, 1:m, data(2:end), '^', 'Color', colours(1,:), 'MarkerSize', 2);
    plot(gh, 1:m, tot*theory(2:end), '-', 'Color', colours(2,:), 'LineWidth', 0.5);
    plot(gh, 1:m, binoinv(0.975, tot, theory(2:end)), '-', 'Color', colours(3,:), 'LineWidth', 0.5);
    plot(gh, 1:m, tot - binoinv(0.975, tot, 1-theory(2:end)), '-', 'Color', colours(3,:), 'LineWidth', 0.5);

    set(fh, 'PaperUnits', 'inches');
    w = 4; h = 3;
    set(fh, 'PaperSize', [w h]);
    set(fh, 'PaperPosition', [0 0 w h]);

    print(fh, '-dpdf', '-painters', sprintf('atRA-posttreat-%d-days', atra_posttreat_ts(i)*7));

    close(fh);

end

end