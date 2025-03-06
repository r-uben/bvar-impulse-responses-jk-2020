% Plot shocks over time (Figure 4)
clear all, close all

tab = readtable("../results/U1bp.csv");
names = tab.Properties.VariableNames(2:end);
U = tab{:, 2:end};
[T, N] = size(U);

% displays for each shock the dates with the largest absolute value
Ntop = 10;
for n = 1:N
    [~, idx] = sort(abs(U(:,n)),'descend');
    idx = idx(1:Ntop);
    fprintf('shock%d\n', n);
    disp(tab(idx,:))
end


lab = readtable('labels.txt');
first2 = find(year(lab.date)>2008,1);
offscl = 2;

pos = [5, 2, 25, 10];
fh = figure('Units','centimeters','Position',pos);
plot(tab.date, U, 'LineWidth', 1.5)
grid on
xlim([datetime(1991,1,1), datetime(2009,1,1)])
ymax = max(ylim);
for i = 1:(first2-1)
    try % label only if the event i is in tab
        labi = join(lab(i,:),tab); % error if i is not in tab
        if contains(labi.txt{1}, 'IM')
            xline(labi.date,':','LineWidth',1)
            text(labi.date, ymax+labi.off*offscl+2, labi.txt, 'HorizontalAlignment', 'center');
        else
            xline(labi.date)
            text(labi.date, ymax-.5, labi.txt, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Rotation', 90);
        end
    end
end
box off
legend(names,'Location','best')
exportgraphics(fh, "U1bp_history_1991_2008.pdf")



pos = [5, 2, 25, 10];
fh2 = figure('Units','centimeters','Position',pos);
plot(tab.date, U, 'LineWidth', 1.5)
grid on
xlim([datetime(2009,1,1), datetime(2020,1,1)])
ymax = max(ylim);
for i = first2:height(lab)
    try % label only if the event i is in tab
        labi = join(lab(i,:),tab); % error if i is not in tab
        xline(labi.date)
        text(labi.date, ymax+labi.off*offscl, labi.txt, 'HorizontalAlignment', 'center');
    end
end
box off
legend(names,'Location','southeast')
exportgraphics(fh2, "U1bp_history_2009_2019.pdf")


