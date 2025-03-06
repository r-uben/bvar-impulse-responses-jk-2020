function result = report_mcmc_nw(chain, n_autocov, varnames, show)
% PURPOSE: Compute numerical standard errors of means from a Markov Chain.
% INPUTS:
% chain - T by N matrix, each row is one draw (one realization of the chain)
% n_autocov - number of autocovariances to use in the Newey-West estimator
%             of the variance of the mean of chain elements
% varnames - N by ? character array, names of the NN elemens of the chain
% show - true/false, show plots and print table?
% RETURNS:
% result.mean - 1 by N mean of the chain
% result.acfs - N by n_autocov matrix, rows contain autocorrelations
%               of elements of the chain up to order n_autocov
% result.S - variance of the chain, Newey-West estimator
%            accounting for the autocorrelation of the chain;
%            note that the variance of the chain mean is S/T
% result.T - number of draws in the chain
% result.RNE - Relative Numerical Efficiency

if nargin<4
    show = true;
end

[T, N] = size(chain);

% 1. compute the mean
mm = mean(chain);

% 2. compute the variance of the mean (Newey-West estimator)
% demean chain
chain = chain - repmat(mm,T,1);

% prepare storage for the autocovariances
Gammas = nan(N, N, n_autocov+1);

% compute autocovariances of order 0 (variance) up to n_autocov
for i = 0:n_autocov
    Gammas(:,:,i+1) = chain(i+1:T,:)'*chain(1:T-i,:)/(T-i);
end

% the Newey-West estimator of the variance
S_NW = Gammas(:,:,1);
for i = 1:n_autocov
    S_NW = S_NW + (1 - i/(n_autocov+1))*(Gammas(:,:,i+1) + Gammas(:,:,i+1)');
end

% autocorrelations
acfs = nan(N,n_autocov+1);
for i = 1:(n_autocov+1)
    acfs(:,i) = diag(Gammas(:,:,i))./diag(Gammas(:,:,1));
end

% 3. some reporting: ACF plot and a table printout
if show
    % plot ACFs
    figure()
    plot(repmat((0:n_autocov)',1,N),acfs')
    legend(varnames, 'Location', 'northeast')

    % print out a table with the results
    ParMean = mm';
    NSE_NW = sqrt(diag(S_NW)/T);
    NSE_iid = sqrt(diag(Gammas(:,:,1))/T);
    RNE = NSE_iid./NSE_NW;
    lastcor = acfs(:,end);

    ParMean = round(ParMean, 3);
    NSE_NW = round(NSE_NW, 1, 'significant');
    NSE_iid = round(NSE_iid, 1, 'significant');
    RNE = round(RNE, 3);
    lastcor = round(lastcor, 2);

    disp(table(ParMean, NSE_NW, NSE_iid, RNE, lastcor, 'RowNames', varnames))
end

% 4. write results to be returned
result.mean = mm;
result.acfs = acfs;
result.S = S_NW;
result.T = T;
result.RNE = sqrt(diag(Gammas(:,:,1))./diag(S_NW));
result.sd = sqrt(diag(Gammas(:,:,1)));
end