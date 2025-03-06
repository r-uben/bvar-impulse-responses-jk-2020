% Geweke convergence diagnostics
clear all, close all
load("../results/chain.mat")

[N, ~, ndraws] = size(chain.W);

Wmat = permute(chain.W, [3, 2, 1]);
Wmat = reshape(Wmat, ndraws, N^2);

drawsPar = [chain.v Wmat];
drawsParnames = convertStringsToChars(["v" + (1:4),"W" + (1:16)]);

n_autocov = 10;
p1 = 0.4; p2 = 0.5;
[ndraws, nPar] = size(drawsPar);

show = true;
if 0
    for n = 1:nPar
        figure;
        plot(drawsPar(:,n))
        title(drawsParnames{n})
    end
end

% split the chain into 1st p1 percent and last p2 percent
nobs1 = round(p1*ndraws);
nobs2 = round(p2*ndraws);
draws1 = drawsPar(1:nobs1,:);
draws2 = drawsPar(nobs2:end,:); 
fprintf('first %.2f\n', p1)
res1 = report_mcmc_nw(draws1, n_autocov, drawsParnames, show);
fprintf('from %.2f\n', p2)
res2 = report_mcmc_nw(draws2, n_autocov, drawsParnames, show);

% test equality of the means for each parameter
SS = res1.S/res1.T + res2.S/res2.T;
mean1 = res1.mean';
mean2 = res2.mean';
absdif = abs(mean1 - mean2);
difstd = sqrt(diag(SS));
Z = absdif./difstd;
pvalue = erfc(Z/sqrt(2)); % normcdf, see help for erfc

% round the numbers before displaying
mean1 = round(mean1, 3);
mean2 = round(mean2, 3);
absdif = round(absdif, 4);
difstd = round(difstd, 2, "significant");
Z = round(Z, 2);
pvalue = round(pvalue, 2);
tab = table(mean1, mean2, absdif, difstd, Z, pvalue, RowNames=drawsParnames);

disp(tab)
disp(tab(tab.Z>2, :))

disp([num2str(nPar) ' parameters, ' num2str(sum(Z>2)) ' Z-stats > 2 (share: ' num2str(sum(Z>2)/nPar) ')'])


