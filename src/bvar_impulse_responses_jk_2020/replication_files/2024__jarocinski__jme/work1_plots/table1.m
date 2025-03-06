%% Produce Table 1: Correlations between alternative sets of shocks
clear all, close all

tab1 = readtable("shocks_tocompare/U1s.csv");
tab1.Properties.VariableNames = cellfun(@(s) strrep(s, 'u','tab1_u'), tab1.Properties.VariableNames', 'Uniform', 0);

shockpath = "shocks_tocompare/";
fnames = shockpath + ["U1s_drop3Jan2001.csv", "U1s_dropQE1.csv", "U1s_1991_2004.csv", "U1s_2005_2019.csv"];

for pp = 1:length(fnames)
    fname2 = fnames(pp);
    tab2 = readtable(fname2);
    tab2.Properties.VariableNames = cellfun(@(s) strrep(s, 'u','tab2_u'), tab2.Properties.VariableNames', 'Uniform', 0);    
    tab = innerjoin(tab1, tab2);
    disp(fname2)
    fprintf('N.obs: %d \n', size(tab,1))
    for n = 1:4
        x1 = tab{:,{sprintf('tab1_u%d',n)}};
        x2 = tab{:,{sprintf('tab2_u%d',n)}};
        fprintf('Shock%d: %.3f %.3f\n', n, corr(x1,x2,'Type','Spearman'), corr(x1,x2))
    end
    disp(' ')
end

%% Compare with Swanson (2021)
fname2 = shockpath + "swanson.csv";
tab2 = readtable(fname2);
tab = innerjoin(tab1, tab2);

disp(fname2)
fprintf('N.obs: %d \n', size(tab,1))
for n = 1:3
    x1 = tab{:,n+1};
    x2 = tab{:,n+5};
    fprintf('Shock%d: %.3f %.3f\n', n, corr(x1,x2,'Type','Spearman'), corr(x1,x2))
end
x1 = tab.tab1_u4;
x2 = tab.fg;
fprintf('Shock%d: %.3f %.3f\n', n, corr(x1,x2,'Type','Spearman'), corr(x1,x2))
disp(' ')

%% Compare with Jarocinski, Karadi (2020) FF4 shocks
fname2 = shockpath + "data_fig1_median_d.csv";
tab2 = readtable(fname2);
tab2.date = datetime(tab2.year, tab2.month, tab2.day);
tab = innerjoin(tab1, tab2);

disp(fname2)
fprintf('N.obs: %d \n', size(tab,1))
for n = 1:3
    x1 = tab{:,sprintf('tab1_u%d',n)};
    x2 = tab.MP_med;
    fprintf('Shock%d: %.3f %.3f\n', n, corr(x1,x2,'Type','Spearman'), corr(x1,x2))
end
x1 = tab.tab1_u4;
x2 = tab.CBI_med;
fprintf('Shock%d: %.3f %.3f\n', n, corr(x1,x2,'Type','Spearman'), corr(x1,x2))
disp(' ')


%% Compare with Jarocinski, Karadi (2020) principal component shocks
fname2 = shockpath + "shocks_fed_gss_pr_median_d.csv";
tab2 = readtable(fname2);
tab2.date = datetime(tab2.year, tab2.month, tab2.day);
tab = innerjoin(tab1, tab2);

disp(fname2)
fprintf('N.obs: %d \n', size(tab,1))
for n = 1:3
    x1 = tab{:,sprintf('tab1_u%d',n)};
    x2 = tab.MP_median;
    fprintf('Shock%d: %.3f %.3f\n', n, corr(x1,x2,'Type','Spearman'), corr(x1,x2))
end
x1 = tab.tab1_u4;
x2 = tab.CBI_median;
fprintf('Shock%d: %.3f %.3f\n', n, corr(x1,x2,'Type','Spearman'), corr(x1,x2))
disp(' ')


fh = figure;
plot(tab.date, tab.CBI_median*50, tab.date, tab.tab1_u4);
grid on
legend('CBI (JK 2020, Appendix)', 'u4')


