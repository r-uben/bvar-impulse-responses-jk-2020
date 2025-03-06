% Plots for comparing this estimation with the Student-t
% Run this after main3_plots.m
load("baseline_tocompare/C1s.mat");
C1s_Studentt = C1s;

[fh, varminmax] = plot_resp_two(maxl.C1s, C1s_l1, C1s_u1, C1s_Studentt.l, C1s_Studentt.u,...
    ymaturities, ynames);
exportgraphics(fh, out_path + "C1s_band2.pdf")
