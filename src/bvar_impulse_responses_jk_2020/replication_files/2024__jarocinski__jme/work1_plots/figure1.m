% Script that plots the histograms of the variables (Figure 1)
clear all, close all

tab = readtable('../work/data.csv');
Y = tab{:,2:end};
ynames = tab.Properties.VariableNames(2:end);

nbins = [61, 35, 25, 25];
for n = 1:4
    x = Y(:,n);
    
    pd_normal = fitdist(x, 'normal');
    pd_normal.mu = 0;
    pd_t = fitdist(x, 'tlocationscale');
    pd_t.mu = 0;
    
    xmin = min([x;-x]); xmax = max([x;-x]);
    xval = linspace(xmin, xmax)*1.1;
    xval3sig = linspace(-3*pd_normal.sigma, 3*pd_normal.sigma);
    n_pdf = pdf(pd_normal, xval3sig);
    t_pdf = pdf(pd_t, xval);
    
    fh = figure('Units','centimeters','Position',[5 3, 10, 6]);
    hold on
    hh = histogram(x, nbins(n), 'BinLimits',[xmin xmax], 'Normalization','pdf', 'LineWidth',1.2, 'FaceColor',[0 0.4470 0.7410]);
    ln = plot(xval3sig, n_pdf, '-r', 'LineWidth', 3);
    lt = plot(xval, t_pdf, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.2);
    uistack(ln)
    axis tight
    xlim([xmin xmax]*1.1)    
    xlabel('basis points')
    legend([hh ln lt], {'data', sprintf('normal, std=%.0f', pd_normal.sigma), sprintf('t-student, v=%.1f', pd_t.nu)}, 'Box', 'off', 'Location','northwest')
    title(ynames{n})
    exportgraphics(fh, [ynames{n} '.pdf'])
end


