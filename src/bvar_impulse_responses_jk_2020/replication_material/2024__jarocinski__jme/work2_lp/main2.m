% Estimate local projections for a distribution of shocks
clear all, close all

% determine the shocks to use
load("../results/chain.mat")
load("../results/data.mat")
load("../results/maxl.mat")
drawsW = cat(3, maxl.W, chain.W);
Y = tab{:,2:5}; % /100 -> shocks in percent
shocknames = ["u1","u2","u3","u4"];

% determine the list of right-hand side variables
fig5 = ["ff","sveny02","sveny10","sp500","bofaml_us_hyld_oas","bkeven05","eurusd"]; % Figure 5
figI1 = ["krippner","krippnerlong","gs3mo"]; % Figure I.1
figI2 = ["bofaml_us_aa_yld","bofaml_us_aa_oas","bofaml_us_bbb_yld","bofaml_us_bbb_oas","nfci_lin"]; % Figure I.2
figI3 = ["bkeven02","bkeven10","bkeven1f4","bkeven5f05"]; % Figure I.3
figI4 = ["broad","afe","eme"]; % Figure I.4
varnames = fig5;
% uncomment for the Appendix figures
varnames = [figI1 figI2 figI3 figI4];

% load and join daily
tabd = readtable('../daily/daily.csv');
tabd.date = datetime(tabd.date, 'InputFormat', 'ddMMMyyyy');
tab = join(tab, tabd, "Keys", "date");
clear tabd

% parameters of the plots
tt = [1, 2, 3, 4, 5, 10, 15, 20, 25]';
qtoplot = [0.5 0.16 0.84 0.05 0.95];
sdtoplot = [0, -1, 1, -1.645, 1.645];

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
ndraws = size(drawsW,3);

%% Estimate

fnames = strings(length(varnames),1);

for vv = 1:length(varnames)
    varname = varnames(vv);
    nicename = nicenames(vv);

    stdincreases = nan(length(tt),length(shocknames));

    for ss = 1:4
        shockname = shocknames(ss);

        % drawing shocks, then coefs|shocks
        bhat_draws = nan(ndraws, length(tt));
        var_draws = nan(ndraws, length(tt));

        tic
        for ii = 1:ndraws

            % draw shock
            UU = Y*drawsW(:,ss,ii);
            UU = UU/diag(std(UU));

            % estimate LP for each horizon
            for hh = 1:length(tt)
                yname = "f" + tt(hh) + "l1" + varname;
                y = tab.(yname);

                % drop missing observations
                UUU = UU(not(isnan(y)),:);
                y = y(not(isnan(y)));
                T = length(y);
                X = [UUU ones(T,1)];
                K = size(X,2);

                % estimate coefs conditional on the drawn shocks
                bhat = X\y;
                e2 = (y-X*bhat).^2;
                XXinv = (X'*X)\eye(K);
                EHWcov = XXinv*X'*diag(e2)*X*XXinv;

                bhat_draws(ii,hh) = bhat(1);
                var_draws(ii,hh) = EHWcov(1,1);
            end
        end
        toc

        % conditional on maximum likelihood point estimates of the shocks
        toplot1 = bhat_draws(1,:)' + sqrt(var_draws(1,:))'*sdtoplot;

        % using law of total variance
        totalvar = mean(var_draws) + var(bhat_draws);
        toplot2 = mean(bhat_draws)' + sqrt(totalvar)'*sdtoplot; %

        % band for the point estimates
        toplot3 = mean(bhat_draws)' + std(bhat_draws)'*sdtoplot; %
        
        toplot = [toplot1(:,1) toplot2(:,2:5)];

        plot_lp_daily

        title(shockname + " $\rightarrow$ " + nicename, 'FontWeight', 'normal', 'Interpreter','latex')
        fname = outdir + varname + "-" + shockname + ".pdf";
        exportgraphics(fh, fname)

        stdinc = sqrt(totalvar)./sqrt(var_draws(1,:));
        stdincreases(:,ss) = stdinc(:);
        disp(max(stdinc(:)))
    end
    writematrix(stdincreases, outdir + varname + "_stdinc.csv")
end


%% optional: report the std increase factors due to uncertainty about U
if 0
allstdinc = [];
for vv = 1:length(varnames)
    disp(varnames(vv))
    stdinc = readmatrix("out/" + varnames(vv) + "_stdinc.csv");
    disp(max(stdinc))
    allstdinc = [allstdinc; stdinc];
end
allsdinc = allstdinc(:);

fh = figure;
histogram(allstdinc, 10)
exportgraphics(fh, "histogram_stdinc.pdf")

mean(allsdinc)
quantile(allsdinc, [0.5; 0.9; 0.95; 1])
end