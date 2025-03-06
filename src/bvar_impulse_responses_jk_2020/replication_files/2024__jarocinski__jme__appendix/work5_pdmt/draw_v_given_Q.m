function v = draw_v_given_Q(Q, vlast, va, vb)
% PURPOSE: Draw v from p(v|Q) where
% Q is a matrix of i.i.d. Gamma(v/2, 2/v)
% with the prior p(v)=Gamma(va, vb),
% using Metropolized Independence Sampler.
% INPUTS:
% Q - matrix (array) with positive numbers
% v - the previous draw of v
% va, vb - scalars, parameters of the Gamma(va, vb) prior
% OUTPUT:
% v - draw from p(v|Q) (could be the same as vlast)
% DEPENDS: drawmisgamma_modecurv
%
% Marek Jarocinski, 2021-Sep.

N = numel(Q);

% logp(v|.)
d1 = -1/vb - 0.5*sum(Q,'all') -0.5*N*log(2) + 0.5*sum(log(Q),'all');
logp = @(x) (va-1)*log(x) + d1*x + 0.5*N*x.*log(x) - N*gammaln(0.5*x);

% dlogp(v|.)/dv
dlogp = @(x) (va-1)*x.^(-1) + d1 - 0.5*N*psi(0.5*x) + 0.5*N*(log(x)+1);

% mode and curvature at the mode
try
    mm = fzero(dlogp, [1e-4 40]);
catch
    warning('v>40')
    mm = fzero(dlogp, [1e-4 400]);
end
ll = -(va-1)*mm^-2 + 0.5*N/mm - 0.25*N*psi(1,mm/2);

v = drawmisgamma_modecurv(logp, mm, ll, vlast, 0);

end