function [xnew, logacceptprob] = drawmisgamma_modecurvnum(logp, xlast, b_plot)
% PURPOSE: Draw xnew from the distribution given by logp (univariate)
% using the Metropolized Independence Sampler based on Gamma approximation.
% The proposal gamma density has mode and curvature of the target density.
% Uses numerical optimization (derivest) to find the mode and curvature.
% INPUTS:
% logp - function evaluating the log target density
% xlast - previous draw
% b_plot - TRUE/FALSe: plot the approximation?
% OUTPUTS:
% xnew - draw of x (might be equal to the previous draw)
% logacceptprob - log of the acceptance probability used for accepting
% DEPENDS: derivest.m
% REFERENCE: Jiang and Ding (2016) J of Stat. Planning and Inference
% Marek Jarocinski 2021-May, 2022-Mar
if 0
    options = optimoptions('fminunc');
    options.Display = 'off';
    fun = @(x) -logp(x);
    [mm,fval,exitflag,output,grad,hessian] = fminunc(fun, xlast, options);
    ll = -hessian;
else
    options = optimoptions('fmincon');
    options.Display = 'off';
    fun = @(x) -logp(x);
    mm = fmincon(fun,xlast,[],[],[],[],0,[],[],options);
    hessian = derivest(fun,mm,'der',2,'FixedStep',0.1);
    ll = -hessian;
end

if isnan(mm) % exponential approximating density
    aa = 1; bb = ll; mm = 3; % set some mm for scaling the plot
else % find the parameters of the gamma approximating density
    bb = -1/(mm*ll); aa = 1-ll*mm^2;
end

if b_plot % plot the approximation
    acceptprob = @(xprime) exp(logp(xprime) - logp(xlast) + log(gampdf(xlast, aa, bb)) - log(gampdf(xprime, aa, bb)));
    xx = 1e-4:0.01:10;
    figure;
    subplot(2,1,1)
    temp = -max(logp(xx));
    plot(xx, exp(temp+logp(xx)), '-k', 'LineWidth', 2)
    yyaxis right
    plot(xx, gampdf(xx, aa, bb), ':r', 'LineWidth', 2)
    legend('target','approx')
    subplot(2,1,2)
    plot(xx, acceptprob(xx), 'LineWidth', 2);
    xline(xlast)
    yline(1)
    yline(0)
    ylim([0 1.5])
end

xprime = gamrnd(aa, bb);

logacceptprob = logp(xprime) - logp(xlast) + log(gampdf(xlast, aa, bb)) - log(gampdf(xprime, aa, bb));

if log(rand) < logacceptprob
    xnew = xprime;
else
    xnew = xlast;
end
end