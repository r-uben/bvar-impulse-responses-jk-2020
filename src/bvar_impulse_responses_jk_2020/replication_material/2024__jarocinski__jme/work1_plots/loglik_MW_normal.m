function ll = loglik_MW_normal(M, W, ~)
% Evaluate log p(M|W,v)
% Model M = U*inv(W), U~i.i.d. normal
U = M*W;
logpU = -0.5*(U.^2);
ll = size(M,1)*log(det(W)) + sum(sum(logpU));