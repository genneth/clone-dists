function oes_atra_bubbles

ts = [3 6 3+1/7 3+9/7 3+1/7 3+9/7];

% pre-treat/atra homoestasis
databs{1} = sparse([
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
databs{1}(15+1,12+1) = 1;
databs{1}(17+1, 7+1) = 1;
databs{1}(19+1,10+1) = 1;
databs{1}(19+1,20+1) = 1;

databs{2} = sparse([
	0	7	12	2	4	2	0	1	0	0	0	0	0	0	0	0;
	3	15	10	7	9	1	0	0	0	0	0	0	0	0	0	0;
	4	26	12	9	4	1	1	0	1	0	0	0	0	0	0	0;
	4	4	17	10	7	2	2	1	0	0	0	0	1	0	0	0;
	1	5	5	8	6	5	3	0	1	0	0	0	1	0	0	0;
	0	1	3	4	2	6	5	1	1	0	1	0	0	0	0	0;
	0	0	2	1	2	3	2	3	0	1	0	1	0	0	0	0;
	0	0	1	3	3	1	0	0	2	0	1	0	0	0	0	0;
	0	0	0	1	0	0	2	0	1	0	0	1	0	0	0	0;
	0	0	0	0	0	0	5	1	0	1	0	0	1	0	0	0;
	0	0	0	0	2	0	0	1	0	3	0	0	1	1	0	0;
	0	0	0	0	0	0	3	1	1	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0;
	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0]);
databs{2}(20+1,20+1) = 1;
databs{2}(21+1,13+1) = 1;
databs{2}(25+1,21+1) = 1;
databs{2}(46+1,40+1) = 1;
databs{2}(51+1,26+1) = 1;

% atra treatement, 3 weeks after induction
databs{3} = sparse([
	0	23	29	17	8	2	2	0	0	0	0;
	14	48	23	7	2	2	0	2	0	0	0;
	15	32	23	15	6	3	0	0	0	0	0;
	8	11	22	15	3	1	1	0	0	0	0;
	0	6	15	10	5	5	0	0	0	0	0;
	1	4	5	3	1	1	0	0	0	0	0;
	0	1	1	3	7	1	0	0	0	0	0;
	0	0	0	0	3	3	0	0	0	0	0;
	0	0	0	4	1	1	0	0	0	0	0;
	0	0	0	0	1	0	0	0	0	0	0;
	0	0	0	0	2	2	0	0	0	1	0;
	0	0	0	1	0	0	0	0	0	0	0;
	0	0	0	0	2	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	1	0	0	0	1;
	0	0	0	0	0	0	1	0	2	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	1	0	0	0	0	0]);

databs{4} = sparse([
	0	22	29	21	11	11	6	1	0	0	0	0	1	0	0	0	0	0	0	0;
	5	24	14	19	11	3	5	1	0	1	0	0	0	0	0	0	0	0	0	0;
	2	19	19	18	10	4	3	1	0	0	0	0	0	0	0	0	0	0	0	0;
	2	3	6	9	8	5	3	0	1	0	0	0	0	0	0	0	0	0	0	0;
	0	1	9	8	4	6	4	2	1	2	0	0	0	0	0	0	0	0	0	1;
	0	0	2	7	2	4	1	2	2	1	0	0	0	0	0	0	0	0	0	0;
	0	0	0	2	3	2	1	0	1	1	0	1	0	0	0	0	0	0	0	0;
	0	0	2	3	1	1	0	1	0	0	1	0	0	0	0	0	0	0	0	0;
	0	0	1	2	2	0	0	2	1	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	1	0	1	0	0	1	1	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1;
	0	0	0	0	0	1	0	1	0	0	0	1	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	1	0	0	1	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0]);

% control
databs{5} = sparse([
	0	69	39	14	7	4	1	0;
	32	59	24	6	3	0	0	0;
	32	33	17	15	1	1	0	0;
	5	11	10	2	3	1	0	0;
	4	6	3	1	1	0	1	0;
	0	1	1	2	1	0	1	0;
	1	4	3	2	0	1	0	1;
	0	1	0	0	0	0	1	0;
	0	0	0	0	1	0	0	0;
	0	0	1	0	1	0	0	0;
	0	0	0	1	0	0	0	0;
	0	0	0	0	1	0	0	0]);

databs{6} = sparse([
	0	29	26	9	3	2	1	0	0	0	0;
	12	35	32	8	9	1	1	0	0	0	0;
	12	26	23	12	8	2	2	1	0	0	0;
	3	8	6	10	7	2	0	0	1	0	0;
	3	5	9	6	7	3	0	0	0	1	0;
	1	0	2	9	6	2	1	0	0	0	0;
	0	2	1	1	0	2	1	0	0	0	0;
	0	1	0	1	1	2	1	3	1	0	0;
	0	0	1	1	1	2	1	1	0	1	0;
	0	0	0	0	2	1	0	1	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	0	0	0	0;
	0	0	0	0	0	0	0	1	0	0	0;
	0	0	0	0	0	0	0	0	0	0	1]);

    function plot_slice(ah, cx, cy, r, theta, colour, linespec, transparency)
        X = cx + r * cos(linspace(theta(1), theta(2), 20));
        Y = cy + r * sin(linspace(theta(1), theta(2), 20));
        axes(ah);
        fill_h = fill(X, Y, linespec);
        set(fill_h, 'EdgeColor', colour, 'FaceColor', colour, ...
            'EdgeAlpha', 1.0, 'FaceAlpha', transparency, 'LineWidth', 1.2);
        %plot(ah, X, Y, linespec, 'Color', colour, 'LineWidth', 1.2);
    end

    function plot_slices(ah, bs, scale, theta, colour, linespec, transparency)
        [m,n,c] = find(bs);
        tot = sum(sum(bs));
        for i = 1:numel(c)
            if c > 0
                    plot_slice(ah, m(i)-1, n(i)-1, sqrt(c(i) / tot) * scale, theta, colour, linespec, transparency);
            end
        end
    end

for dataset = [1 4]

    t = ts(dataset);
    if dataset == 1
        lambda = 1.4;
    elseif dataset == 4
        lambda = (1.4 * 3 + 0.7676 * (t-3)) / t;
    end
    r = 0.1926;
    rho = 0.5164;
    gamma = lambda * rho / (1-rho);

    theorybs = clone_size_distribution_cond(r, gamma, 0.0, [t/2 t] .* lambda, 7, 7);

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

    databs{dataset}(0+1,:)   = 0;
    databs{dataset}(1+1,0+1) = 0;
    plot_slices(gh, databs{dataset}, 1.0, [0 2*pi], colours(1,:), '-', 1.0);
    plot_slices(gh, squeeze(theorybs(:,:,2)), 1.0, [0 2*pi], colours(4,:), '-', 0.0);
    %plot_slices(gh, databs{5}, 1.0, [3*pi/2 2*pi], colours(2,:), 0.95);
    %plot_slices(gh, databs{6}, 1.0, [0      pi/2], colours(3,:), 0.95);

    set(gh, 'DataAspectRatio', [1 1 1]);
    axes(gh);
    %axis tight;
    set(gh, 'XLim', [0.5 7.5], 'YLim', [-0.5 7.5]);

    xlabel(gh, 'Basal count', 'FontName', 'Times', 'FontSize', 10);
    ylabel(gh, 'Suprabasal count', 'FontName', 'Times', 'FontSize', 10);

    set(gh, 'FontName', 'Times', 'FontSize', 9);

    set(fh, 'PaperUnits', 'inches');
    w = 5; h = 5;
    set(fh, 'PaperSize', [w h]);
    set(fh, 'PaperPosition', [0 0 w h]);

    set(fh, 'Color', 'white');

    drawnow

    if dataset == 1
        print(fh, '-dpng', '-opengl', '-r300', 'oes-atra-pre-treat-bubbles');
    elseif dataset == 4
        print(fh, '-dpng', '-opengl', '-r300', 'oes-atra-post-treat-bubbles');
    end

end

end
