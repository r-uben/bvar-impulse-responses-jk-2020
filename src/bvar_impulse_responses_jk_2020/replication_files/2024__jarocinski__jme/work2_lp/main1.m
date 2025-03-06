% Estimate local projections for the point estimated shocks
clear all, close all

% determine the shocks to use
tabs = readtable("../results/U1s.csv");
shocknames = ["u1","u2","u3","u4"];
tabs{:,shocknames} = tabs{:,shocknames}/100; % shocks in percent

% determine the list of right-hand side variables
fig5 = ["ff","sveny02","sveny10","sp500","bofaml_us_hyld_oas","bkeven05","eurusd"]; % Figure 5
figI1 = ["krippner","krippnerlong","gs3mo"]; % Figure I.1
figI2 = ["bofaml_us_aa_yld","bofaml_us_aa_oas","bofaml_us_bbb_yld","bofaml_us_bbb_oas","nfci_lin"]; % Figure I.2
figI3 = ["bkeven02","bkeven10","bkeven1f4","bkeven5f05"]; % Figure I.3
figI4 = ["broad","afe","eme"]; % Figure I.4
varnames = fig5;
% uncomment for the Appendix figures
%varnames = [figI1 figI2 figI3 figI4];

% load daily
tabd = readtable('../daily/daily.csv');
tabd.date = datetime(tabd.date, 'InputFormat', 'ddMMMyyyy');

% merge
tab = join(tabs, tabd,  "Keys", "date");
clear tabs tabd

% prepare nice names for plot titles
nicenames = varnames;
temp = readtable("nicenames_d.csv", "ReadRowNames", true);
for vv = 1:length(varnames)
    try
        nicenames(vv) = temp.name(varnames(vv));
    catch
        nicenames(vv) = strrep(varnames(vv),'_','-');
    end
end

outdir = "out/";
mkdir(outdir);

%% Estimate
tt = [1, 2, 3, 4, 5, 10, 15, 20, 25]';

fnames = strings(length(shocknames),1);

for vv = 1:length(varnames)
    varname = varnames(vv);
    fhh = [];
    for ss = 1:length(shocknames)
        shockname = shocknames(ss);
        Xnames = shockname;

        rowNames = ["b1","s1","R-sq","N.obs."];
        restab = array2table(nan(length(rowNames), length(tt)), 'RowNames', rowNames, 'VariableNames', "f" + tt' + "l1" + varname);

        % estimate LP for each horizon
        for hh = 1:length(tt)
            yname = "f" + string(tt(hh)) + "l1" + varname;

            mdl = fitlm(tab(:,[Xnames yname]));
            if isnan(mdl.Rsquared.Ordinary)
                coeff = zeros(length(Xnames)+1,1); EHWse = zeros(length(Xnames)+1,1);
            else
                [EHWcov,EHWse,coeff] = hac(mdl, 'type', 'HC', 'display', 'off');
            end

            res_h = [coeff(2), EHWse(2), mdl.Rsquared.Ordinary, mdl.NumObservations];

            restab.(yname) = res_h';
        end

        writetable(restab, outdir + varname + "-" + shockname + ".csv");

        % plot
        toplot = restab{1,:}' + restab{2,:}'*[0, -1, 1, -1.645, 1.645];

        plot_lp_daily
        title(shockname + " $\rightarrow$ " + nicenames(vv), 'FontWeight', 'normal', 'Interpreter','latex')
        fname = outdir + varname + "-" + shockname + ".pdf";
        exportgraphics(fh, fname)
        fnames(ss) = fname;

        fhh = [fhh fh];
    end

    if 0 % align y axes all shocks
        sameaxes('y',fhh);
        for ss = 1:length(fhh)
            exportgraphics(fhh(ss), fnames(ss))
        end
    end
    
end


